/*!
  @file   ito_plasma_godunov.cpp
  @author Robert Marskar
  @date   June 2020
  @brief  Implementation of ito_plasma_godunov
*/

#include "ito_plasma_godunov.H"
#include "data_ops.H"
#include "units.H"
#include "poisson_multifluid_gmg.H"

#include <ParmParse.H>

#define DEBUG 1

using namespace physics::ito_plasma;

ito_plasma_godunov::ito_plasma_godunov(){
  m_name = "ito_plasma_godunov";
  m_use_old_dt = false;
}

ito_plasma_godunov::ito_plasma_godunov(RefCountedPtr<ito_plasma_physics>& a_physics){
  m_name    = "ito_plasma_godunov";
  m_physics = a_physics;
  m_use_old_dt = false;
}

ito_plasma_godunov::~ito_plasma_godunov(){

}

void ito_plasma_godunov::parse_options() {
  CH_TIME("ito_plasma_godunov::parse_options");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_options" << endl;
  }

  ParmParse pp(m_name.c_str());
  std::string str;
  
  pp.get("verbosity",      m_verbosity);
  pp.get("ppc",            m_ppc);
  pp.get("max_cells_hop",  m_max_cells_hop);
  pp.get("merge_interval", m_merge_interval);
  pp.get("relax_factor",   m_relax_factor);
  pp.get("regrid_super",   m_regrid_superparticles);
  pp.get("algorithm",      str);

  if(str == "euler"){
    m_algorithm = which_algorithm::euler;
  }
  else if(str == "predictor_corrector"){
    m_algorithm = which_algorithm::predictor_corrector;
  }
  else if(str == "semi_implicit"){
    m_algorithm = which_algorithm::semi_implicit;
  }
  else{
    MayDay::Abort("ito_plasma_godunov::parse_options - unknown algorithm requested");
  }
}

void ito_plasma_godunov::allocate_internals(){
  CH_TIME("ito_plasma_godunov::allocate_internals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocate_internals" << endl;
  }

  m_amr->allocate(m_J,            m_phase, SpaceDim);
  m_amr->allocate(m_scratch1,     m_phase, 1);
  m_amr->allocate(m_scratch2,     m_phase, 1);
  m_amr->allocate(m_conduct_cell, m_phase, 1);
  m_amr->allocate(m_conduct_face, m_phase, 1);
  m_amr->allocate(m_conduct_eb,   m_phase, 1);
  
}

void ito_plasma_godunov::compute_dt(Real& a_dt, time_code::which_code& a_timecode){
  CH_TIME("ito_plasma_godunov::compute_dt");
  if(m_verbosity > 5){
    pout() << "ito_plasma_godunov::compute_dt" << endl;
  }
  
  a_dt = m_ito->compute_dt();
  a_dt = a_dt*m_max_cells_hop;
  a_timecode = time_code::cfl;

  if(m_algorithm == which_algorithm::euler || m_algorithm == which_algorithm::predictor_corrector){
    const Real dtRelax = m_relax_factor*m_dt_relax;
    if(dtRelax < a_dt){
      a_dt = dtRelax;
      a_timecode = time_code::relaxation_time;
    }
  }
  else if(m_algorithm == which_algorithm::semi_implicit && m_use_old_dt == true){
    a_dt = m_dt;
    m_use_old_dt = false;
  }

#if 1
  if(procID() == 0) std::cout << "dt = " << a_dt << "\t relax dt = " << m_dt_relax << "\t factor = " << a_dt/m_dt_relax << std::endl;
#endif
}
void ito_plasma_godunov::pre_regrid(const int a_lmin, const int a_old_finest_level){
  CH_TIME("ito_plasma_godunov::pre_regrid");
  if(m_verbosity > 5){
    pout() << "ito_plasma_godunov::pre_regrid" << endl;
  }

  ito_plasma_stepper::pre_regrid(a_lmin, a_old_finest_level);

  if(m_algorithm == which_algorithm::semi_implicit){
    m_use_old_dt = true;
  }
}

Real ito_plasma_godunov::advance(const Real a_dt) {
  CH_TIME("ito_plasma_godunov::advance");
  if(m_verbosity > 5){
    pout() << m_name + "::advance" << endl;
  }

  // Particle algorithms
  if(m_algorithm == which_algorithm::euler){
    this->advance_particles_euler(a_dt);
  }
  else if(m_algorithm == which_algorithm::predictor_corrector){
    this->advance_particles_pc(a_dt);
  }
  else if(m_algorithm == which_algorithm::semi_implicit){
    this->advance_particles_si(a_dt);
  }

  // Compute current and relaxation time.
  this->compute_J(m_J, a_dt);
  m_dt_relax = this->compute_relaxation_time();

  // Move photons
  this->advance_photons(a_dt);

  // Sort the particles and photons per cell so we can call reaction algorithms
  m_ito->sort_particles_by_cell();
  this->sort_bulk_photons_by_cell();
  this->sort_source_photons_by_cell();

  // Chemistry kernel.
  this->advance_reaction_network(a_dt);

  // Make superparticles
  if((m_step+1) % m_merge_interval == 0 && m_merge_interval > 0){
    m_ito->make_superparticles(m_ppc);
  }

  // Sort particles per patch. 
  m_ito->sort_particles_by_patch();
  this->sort_bulk_photons_by_patch();
  this->sort_source_photons_by_patch();

  // Clear other data holders for now. BC comes later
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->clear(solver_it()->get_eb_particles());
    solver_it()->clear(solver_it()->get_domain_particles());
  }

  // Prepare next step
  this->compute_ito_velocities();
  this->compute_ito_diffusion();
  
  return a_dt;
}

void ito_plasma_godunov::advance_particles_euler(const Real a_dt){
  CH_TIME("ito_plasma_godunov::advance_particles_euler");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_particles_euler" << endl;
  }

  this->set_old_positions();
  this->advect_particles_euler(a_dt);
  this->diffuse_particles_euler(a_dt);
  
  m_ito->remap();
  m_ito->deposit_particles();
  this->intersect_particles(a_dt);
  this->solve_poisson();
}

void ito_plasma_godunov::advance_particles_pc(const Real a_dt){
  CH_TIME("ito_plasma_godunov::advance_particles_pc");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_particles_pc" << endl;
  }

  // ---------- PREDICTOR BEGIN -------
  this->set_old_positions();
  this->advect_particles_rk2(a_dt);
  this->diffuse_particles_euler(a_dt);
  
  m_ito->remap();
  m_ito->deposit_particles();
  this->solve_poisson();
  // ---------- PREDICTOR END ---------

  // ---------- CORRECTOR BEGIN ---------
  this->rewind_particles();
  m_ito->remap();
  
  this->compute_ito_velocities();
  this->compute_ito_diffusion();

  this->advect_particles_rk2(a_dt);
  this->diffuse_particles_euler(a_dt);
  m_ito->remap();
  this->intersect_particles(a_dt);
  
  m_ito->deposit_particles();
  this->solve_poisson();
  // ---------- CORRECTOR END ---------
}

void ito_plasma_godunov::advance_particles_si(const Real a_dt){
  CH_TIME("ito_plasma_godunov::advance_particles_si");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_particles_si" << endl;
  }

  this->set_old_positions();

  // Compute conductivity and setup poissonp
  this->compute_conductivity();
  this->setup_semi_implicit_poisson(a_dt);

  // Diffuse particles
  this->diffuse_particles_euler(a_dt);

  // Remap and deposit
  m_ito->remap();
  m_ito->deposit_particles();

  // Now compute the electric field
  this->solve_poisson();

  // Compute new ito velocities and advect the particles
  this->compute_ito_velocities();
  this->advect_particles_euler(a_dt);

  // Remap, intersect, and redeposit
  m_ito->remap();
  this->intersect_particles(a_dt);
  m_ito->deposit_particles();
}

void ito_plasma_godunov::advect_particles_euler(const Real a_dt){
  CH_TIME("ito_plasma_godunov::advect_particles_euler");
  if(m_verbosity > 5){
    pout() << m_name + "::advect_particles_euler" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    if(solver->is_mobile()){
      
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	const DisjointBoxLayout& dbl          = m_amr->get_grids()[lvl];
	ParticleData<ito_particle>& particles = solver->get_particles()[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	  List<ito_particle>& particleList = particles[dit()].listItems();
	  ListIterator<ito_particle> lit(particleList);

	  // First step
	  for (lit.rewind(); lit; ++lit){
	    ito_particle& p = particleList[lit];
	    p.oldPosition() = p.position();
	    p.position() += p.velocity()*a_dt;
	  }
	}
      }
    }
  }
}

void ito_plasma_godunov::advect_particles_rk2(const Real a_dt){
  CH_TIME("ito_plasma_godunov::advect_particles_rk2");
  if(m_verbosity > 5){
    pout() << m_name + "::advect_particles_rk2" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    if(solver->is_mobile()){
      
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	const DisjointBoxLayout& dbl          = m_amr->get_grids()[lvl];
	ParticleData<ito_particle>& particles = solver->get_particles()[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	  List<ito_particle>& particleList = particles[dit()].listItems();
	  ListIterator<ito_particle> lit(particleList);

	  // First step
	  for (lit.rewind(); lit; ++lit){
	    ito_particle& p = particleList[lit];
	    p.oldPosition() = p.position();
	    p.position() += 0.5*p.velocity()*a_dt;
	  }
	  // Interpolate velocities
	  solver->interpolate_velocities(lvl, dit());

	  // Second step
	  for (lit.rewind(); lit; ++lit){
	    ito_particle& p = particleList[lit];
	    p.position() = p.oldPosition() + p.velocity()*a_dt;
	  }
	}
      }
    }
  }
}



void ito_plasma_godunov::diffuse_particles_euler(const Real a_dt){
  CH_TIME("ito_plasma_godunov::diffuse_particles_euler");
  if(m_verbosity > 5){
    pout() << m_name + "::diffuse_particles_euler" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    
    if(solver->is_diffusive()){
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	const DisjointBoxLayout& dbl          = m_amr->get_grids()[lvl];
	ParticleData<ito_particle>& particles = solver->get_particles()[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	  List<ito_particle>& particleList = particles[dit()].listItems();
	  ListIterator<ito_particle> lit(particleList);

	  // Diffusion hop.
	  for (lit.rewind(); lit.ok(); ++lit){
	    ito_particle& p = particleList[lit];
	    const RealVect ran = solver->random_gaussian();
	    const RealVect hop = ran*sqrt(2.0*p.diffusion()*a_dt);
	    p.position() += hop;
	  }
	}
      }
    }
  }
}

void ito_plasma_godunov::set_old_positions(){
  CH_TIME("ito_plasma_godunov::set_old_positions()");
  if(m_verbosity > 5){
    pout() << m_name + "::set_old_positions()" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    
    if(solver->is_diffusive()){
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	const DisjointBoxLayout& dbl          = m_amr->get_grids()[lvl];
	ParticleData<ito_particle>& particles = solver->get_particles()[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	  List<ito_particle>& particleList = particles[dit()].listItems();


	  for (ListIterator<ito_particle> lit(particleList); lit.ok(); ++lit){
	    ito_particle& p = particleList[lit];
	    p.oldPosition() = p.position();
	  }
	}
      }
    }
  }
}

void ito_plasma_godunov::rewind_particles(){
  CH_TIME("ito_plasma_godunov::rewind_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::rewind_particles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    
    if(solver->is_diffusive()){
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	const DisjointBoxLayout& dbl          = m_amr->get_grids()[lvl];
	ParticleData<ito_particle>& particles = solver->get_particles()[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	  List<ito_particle>& particleList = particles[dit()].listItems();
	  ListIterator<ito_particle> lit(particleList);

	  // Diffusion hop.
	  for (lit.rewind(); lit; ++lit){
	    ito_particle& p = particleList[lit];
	    p.position() = p.oldPosition();
	  }
	}
      }
    }
  }
}

void ito_plasma_godunov::intersect_particles(const Real a_dt){
  CH_TIME("ito_plasma_godunov::intersect_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::intersect_particles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();

    solver->intersect_particles();
  }
}

void ito_plasma_godunov::compute_conductivity(){
  CH_TIME("ito_plasma_godunov::compute_conductivity");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_conductivity" << endl;
  }

  // Get handle to E gas-side E
  EBAMRCellData Egas;
  m_amr->allocate_ptr(Egas);
  m_amr->alias(Egas, m_phase, m_poisson->get_E());

  // Compute |E| and reset conductivity
  data_ops::vector_length(m_scratch1, Egas);
  data_ops::set_value(m_conduct_cell, 0.0);
  
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<ito_solver>& solver   = solver_it();
    const RefCountedPtr<ito_species>& species = solver->get_species();

    const int q = species->get_charge();
    
    if(solver->is_mobile() &&  q != 0){

      const EBAMRCellData& velo  = solver->get_velo_cell();
      const EBAMRCellData& state = solver->get_state();

      data_ops::vector_length(m_scratch2, velo);
      data_ops::divide_scalar(m_scratch2, m_scratch1);
      data_ops::multiply(m_scratch2, state);
      data_ops::incr(m_conduct_cell, m_scratch2, Abs(q));
    }
  }

  m_amr->average_down(m_conduct_cell, m_phase);
  m_amr->interp_ghost_pwl(m_conduct_cell, m_phase);

  // This code does averaging from cell to face. 
  data_ops::average_cell_to_face_allcomps(m_conduct_face, m_conduct_cell, m_amr->get_domains());

  // This code computes the conductivity on the EB
  const irreg_amr_stencil<eb_centroid_interp>& ebsten = m_amr->get_eb_centroid_interp_stencils(m_phase);
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    ebsten.apply(m_conduct_eb, m_conduct_cell, lvl);
  }
}

void ito_plasma_godunov::setup_semi_implicit_poisson(const Real a_dt){
  CH_TIME("ito_plasma_godunov::setup_semi_implicit_poisson");
  if(m_verbosity > 5){
    pout() << m_name + "::setup_semi_implicit_poisson" << endl;
  }

  poisson_multifluid_gmg* poisson = (poisson_multifluid_gmg*) (&(*m_poisson));

  // Set coefficients as usual
  poisson->set_coefficients();

  // Get bco and increment with mobilities
  MFAMRFluxData& bco   = poisson->get_bco();
  MFAMRIVData& bco_irr = poisson->get_bco_irreg();
  
  EBAMRFluxData bco_gas;
  EBAMRIVData   bco_irr_gas;
  
  m_amr->allocate_ptr(bco_gas);
  m_amr->allocate_ptr(bco_irr_gas);
  
  m_amr->alias(bco_gas,     phase::gas, bco);
  m_amr->alias(bco_irr_gas, phase::gas, bco_irr);

  // Increment with conductivity
  data_ops::incr(bco_gas,     m_conduct_face, units::s_Qe*a_dt/(units::s_eps0));
  data_ops::incr(bco_irr_gas, m_conduct_eb,   units::s_Qe*a_dt/(units::s_eps0));

  // Set up the multigrid solver
  poisson->setup_operator_factory();
  poisson->setup_solver();
  poisson->set_needs_setup(false);
}