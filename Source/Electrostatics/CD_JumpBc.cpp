/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_JumpBc.cpp
  @brief  Implementation of CD_JumpBc.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>

// Our includes
#include <CD_JumpBc.H>
#include <CD_LeastSquares.H>
#include <CD_NamespaceHeader.H>
  
bool JumpBc::s_quadrant_based = true;
int  JumpBc::s_lsq_radius     = 1;

JumpBc::JumpBc(){
  CH_TIME("JumpBc::JumpBc(weak)");
}

JumpBc::JumpBc(const MFLevelGrid&            a_mflg,
	       const LevelData<MFBaseIVFAB>& a_bco,
	       const Real&                   a_dx,
	       const int                     a_order,
	       const LayoutData<IntVectSet>* a_cfivs){
  CH_TIME("JumpBc::JumpBc(full)");

  this->define(a_mflg, a_bco, a_dx, a_order, a_cfivs);
}

JumpBc::~JumpBc(){
  CH_TIME("JumpBc::~JumpBc");
}

void JumpBc::getFirstOrderStencil(Real&             a_weight,
				  VoFStencil&       a_stencil,
				  const VolIndex&   a_vof,
				  const EBISBox&    a_ebisbox,
				  const IntVectSet& a_cfivs){
  const RealVect normal   = a_ebisbox.normal(a_vof);
  const RealVect centroid = a_ebisbox.bndryCentroid(a_vof);

  if(s_quadrant_based){
    EBArith::getLeastSquaresGradSten(a_stencil, a_weight, a_vof, a_ebisbox, m_dx*RealVect::Unit, m_domain, 0);
  }
  else{
    EBArith::getLeastSquaresGradStenAllVoFsRad(a_stencil,
					       a_weight,
					       normal,
					       centroid,
					       a_vof,
					       a_ebisbox,
					       m_dx*RealVect::Unit,
					       m_domain,
					       0,
					       s_lsq_radius);
  }

  // Oh shit, couldn't find a stencil. Try including more cells
  if(s_quadrant_based && a_stencil.size() == 0){
    EBArith::getLeastSquaresGradStenAllVoFsRad(a_stencil,
					       a_weight,
					       normal,
					       centroid,
					       a_vof,
					       a_ebisbox,
					       m_dx*RealVect::Unit,
					       m_domain,
					       0,
					       s_lsq_radius);
  }

  // Damn, still haven't found a stencil. Approximate using one side only
  if(a_stencil.size() == 0){
#if DEBUG_JUMP
    pout() << "JumpBc::get_first_order sten - no sten on domain = "
	   << m_domain << "\t vof = "
	   << a_vof.gridIndex() << endl;
    MayDay::Warning("JumpBc::getFirstOrderStencil - could not find a stencil.");
#endif

    // Make an approximation to the cell-centered gradient
    for (int dir = 0; dir < SpaceDim; dir++){
      VoFStencil sten; 
      EBArith::getFirstDerivStencilWidthOne(sten, a_vof, a_ebisbox, dir, m_dx, (IntVectSet*) &a_cfivs, 0);

      if(sten.size() > 0){
	sten      *= normal[dir];
	a_stencil += sten;
      }
    }
  }
}

bool JumpBc::getSecondOrderStencil(Real&             a_weight,
				   VoFStencil&       a_stencil,
				   const VolIndex&   a_vof,
				   const EBISBox&    a_ebisbox,
				   const IntVectSet& a_cfivs){
  a_stencil.clear();
  bool drop_order = false;

  Vector<VoFStencil> point_stencils;
  Vector<Real> distance_along_lines;
  
  EBArith::johanStencil(drop_order, point_stencils, distance_along_lines, a_vof, a_ebisbox, m_dx*RealVect::Unit, a_cfivs);
  if(drop_order){
    return true;
  }

  const Real& x1   = distance_along_lines[0];
  const Real& x2   = distance_along_lines[1];
  const Real denom = x2*x2*x1 - x1*x1*x2;

  VoFStencil& phi1Sten = point_stencils[0];
  VoFStencil& phi2Sten = point_stencils[1];
  
  phi1Sten *= -x2*x2/denom;
  phi2Sten *=  x1*x1/denom;

  a_weight   = -x1*x1/denom + x2*x2/denom;
  a_stencil +=  phi1Sten;
  a_stencil +=  phi2Sten;

  return false;
}

void JumpBc::define(const MFLevelGrid&            a_mflg,
		    const LevelData<MFBaseIVFAB>& a_bco,
		    const Real&                   a_dx,
		    const int                     a_order,
		    const LayoutData<IntVectSet>* a_cfivs){
  CH_TIME("JumpBc::define");
  m_mflg   = a_mflg;
  m_dx     = a_dx;
  m_domain = m_mflg.getDomain();
  m_multifluidIndexSpace   = m_mflg.getMfIndexSpace();
  m_grids  = m_mflg.getGrids();
  m_order  = a_order;
  m_cfivs  = a_cfivs;

  m_ivs.define(m_grids);
  m_bco.define(m_grids);
  m_weights.define(m_grids);
  m_stencils.define(m_grids);
  m_inhomo.define(m_grids);
  m_homog.define(m_grids);
  m_vofit_gas.define(m_grids);
  m_vofit_sol.define(m_grids);
  m_avgWeights.define(m_grids);
  m_avgStencils.define(m_grids);
  m_avgBco.define(m_grids);
  m_avgJump.define(m_grids);
  m_avgFactor.define(m_grids);

  int num = 0;

  // If we only have one phase, there are no jump cells and we don't need to do anything. 
  if(m_multifluidIndexSpace->numPhases() > 1){
    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit){
      MFInterfaceFAB<Real>& bco            = m_bco[dit()];
      MFInterfaceFAB<Real>& weights        = m_weights[dit()];
      MFInterfaceFAB<Real>& inhomo         = m_inhomo[dit()];
      MFInterfaceFAB<Real>& homog          = m_homog[dit()];
      MFInterfaceFAB<Real>& avgWeights     = m_avgWeights[dit()];
      MFInterfaceFAB<Real>& avgBco         = m_avgBco[dit()];
      MFInterfaceFAB<Real>& avgJump        = m_avgJump[dit()];
      MFInterfaceFAB<Real>& avgFactor      = m_avgFactor[dit()];
      MFInterfaceFAB<VoFStencil>& stens    = m_stencils[dit()];
      MFInterfaceFAB<VoFStencil>& avgStens = m_avgStencils[dit()];

      bco.define(m_mflg,        dit());
      homog.define(m_mflg,      dit());
      stens.define(m_mflg,      dit());
      inhomo.define(m_mflg,     dit());
      avgBco.define(m_mflg,     dit());
      avgJump.define(m_mflg,    dit());
      avgFactor.define(m_mflg,  dit());
      weights.define(m_mflg,    dit());
      avgStens.define(m_mflg,   dit());
      avgWeights.define(m_mflg, dit());

      m_ivs[dit()] = bco.getIVS();
    }

    this->defineVofIterator();
    this->setBcoefficient(a_bco);
    this->buildStencils();
  }

  m_defined = true;
}

void JumpBc::defineVofIterator(){
  CH_TIME("JumpBc::defineVofIterator");

  const int phase1 = 0;
  const int phase2 = 1;
  
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit){
    VoFIterator& vofit_g = m_vofit_gas[dit()];
    VoFIterator& vofit_s = m_vofit_sol[dit()];

    const IntVectSet& ivs = m_bco[dit()].getIVS();
    const EBGraph& graph1 = m_bco[dit()].getIVFAB(phase1).getEBGraph();
    const EBGraph& graph2 = m_bco[dit()].getIVFAB(phase2).getEBGraph();

    vofit_g.define(ivs, graph1);
    vofit_s.define(ivs, graph2);
  }
}

void JumpBc::setBcoefficient(const LevelData<MFBaseIVFAB>& a_bco){
  CH_TIME("JumpBc::setBcoefficient");

  const int comp = 0;
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit){
    for (int iphase = 0; iphase < m_multifluidIndexSpace->numPhases(); iphase ++){

      BaseIVFAB<Real>& bco             = m_bco[dit()].getIVFAB(iphase);
      const BaseIVFAB<Real>& bco_irreg = a_bco[dit()].getIVFAB(iphase);

      for (VoFIterator vofit(bco.getIVS(), bco.getEBGraph()); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	bco(vof, comp) = bco_irreg(vof, comp);
      }
    }
  }
}

void JumpBc::buildStencils(){
  CH_TIME("JumpBc::buildStencils");

  const int comp = 0;

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit){
    for (int iphase = 0; iphase < m_multifluidIndexSpace->numPhases(); iphase ++){

      BaseIVFAB<Real>& bco               = m_bco[dit()].getIVFAB(iphase);
      BaseIVFAB<Real>& avgBco            = m_avgBco[dit()].getIVFAB(iphase);
      BaseIVFAB<Real>& weights           = m_weights[dit()].getIVFAB(iphase);
      BaseIVFAB<Real>& avgWeights        = m_avgWeights[dit()].getIVFAB(iphase);
      BaseIVFAB<VoFStencil>& stencils    = m_stencils[dit()].getIVFAB(iphase);
      BaseIVFAB<VoFStencil>& avgStencils = m_avgStencils[dit()].getIVFAB(iphase);

      avgBco.setVal(0.0);
      avgWeights.setVal(0.0);

      const EBLevelGrid& eblg = m_mflg.getEBLevelGrid(iphase);
      const EBISLayout ebisl  = eblg.getEBISL();

      const EBISBox& ebisbox  = ebisl[dit()];
      const IntVectSet& cfivs = (*m_cfivs)[dit()];

      // Get stencils just like we would do for regular BCs. 
      for (VoFIterator vofit(weights.getIVS(), weights.getEBGraph()); vofit.ok(); ++vofit){
	const VolIndex& vof     = vofit();
	Real& cur_weight        = weights(vof, comp);
	VoFStencil& cur_stencil = stencils(vof, comp);

	bool drop_order = false;
	
	if(m_order == 2){
	  drop_order = this->getSecondOrderStencil(cur_weight, cur_stencil, vof, ebisbox, cfivs);
	}
	else if(m_order == 1){
	  this->getFirstOrderStencil(cur_weight, cur_stencil, vof, ebisbox, cfivs);
	}

	if(m_order == 2 && drop_order){
	  this->getFirstOrderStencil(cur_weight, cur_stencil, vof, ebisbox, cfivs);
	}
      }

      // Now build the average stencils, weights, and coefficients
      const IntVectSet& ivs = avgWeights.getIVS();
      for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit){
	const IntVect iv = ivsit();
	
	const VolIndex vof0(iv, 0);
	const Vector<VolIndex> vofs = ebisbox.getVoFs(iv);

	Real& curBco           = avgBco(vof0, comp);
	Real& curWeight        = avgWeights(vof0, comp);
	VoFStencil& curStencil = avgStencils(vof0, comp);

	curBco    = 0.0;
	curWeight = 0.0;
	curStencil.clear();
	
	for (int ivof = 0; ivof < vofs.size(); ivof++){
	  const VolIndex& vof = vofs[ivof];

	  curBco     += bco(vof, comp);
	  curWeight  += weights(vof, comp);
	  curStencil += stencils(vof, comp);
	}

	const Real invNum = 1./vofs.size();

	curBco     *= invNum;
	curWeight  *= invNum;
	curStencil *= invNum;


#if DEBUG_JUMP
	if(std::isnan(curWeight)) MayDay::Abort("jumpc_bc::buildStencils - got NaN weight");
	if(std::isnan(curBco))    MayDay::Abort("jumpc_bc::buildStencils - got NaN bco");
#endif
      }
    }

    

    // Build the average factor on each side.
    const int phase1 = 0;
    const int phase2 = 1;

    BaseIVFAB<Real>& factor1 = m_avgFactor[dit()].getIVFAB(phase1);
    BaseIVFAB<Real>& factor2 = m_avgFactor[dit()].getIVFAB(phase2);

    const BaseIVFAB<Real>& avgBco1 = m_avgBco[dit()].getIVFAB(phase1);
    const BaseIVFAB<Real>& avgBco2 = m_avgBco[dit()].getIVFAB(phase2);

    const BaseIVFAB<Real>& avgWeights1 = m_avgWeights[dit()].getIVFAB(phase1);
    const BaseIVFAB<Real>& avgWeights2 = m_avgWeights[dit()].getIVFAB(phase2);

    factor1.setVal(0.0);
    factor2.setVal(0.0);
    
    for (IVSIterator ivsit(m_ivs[dit()]); ivsit.ok(); ++ivsit){
      const IntVect iv = ivsit();

      const VolIndex vof0(iv, 0);
      
      const Real& bco1 = avgBco1(vof0, 0);
      const Real& bco2 = avgBco2(vof0, 0);

      const Real& w1 = avgWeights1(vof0, 0);
      const Real& w2 = avgWeights2(vof0, 0);

      const Real factor = 1.0/(bco1*w1 + bco2*w2);

#if DEBUG_JUMP
      if(std::isnan(factor)) MayDay::Abort("JumpBc::buildStencils -- factor is NaN");
#endif

      factor1(vof0, 0) = factor;
      factor2(vof0, 0) = factor;
    }
  }
}

void JumpBc::matchBc(LevelData<BaseIVFAB<Real> >&       a_phibc,
		     const LevelData<BaseIVFAB<Real> >& a_jump,
		     const LevelData<MFCellFAB>&        a_phi,
		     const bool                         a_homogeneous){
  CH_TIME("JumpBc::matchBc(1)");

  DataIterator dit = a_phibc.dataIterator();
  this->matchBc(a_phibc, a_jump, a_phi, dit, a_homogeneous);
}

void JumpBc::matchBc(LevelData<BaseIVFAB<Real> >&       a_phibc,
		     const LevelData<BaseIVFAB<Real> >& a_jump,
		     const LevelData<MFCellFAB>&        a_phi,
		     DataIterator&                      a_dit,
		     const bool                         a_homogeneous){
  CH_TIME("JumpBc::matchBc(1)");

  for (a_dit.reset(); a_dit.ok(); ++a_dit){
    this->matchBc(a_phibc[a_dit()], a_jump[a_dit()], a_phi[a_dit()], a_dit());
  }
}



void JumpBc::computeDnPhi(Vector<LevelData<BaseIVFAB<Real> > >&       a_dphidn,
			  const Vector<LevelData<BaseIVFAB<Real> > >& a_phibc,
			  const LevelData<MFCellFAB>&                 a_phi){
  for (int iphase = 0; iphase < m_multifluidIndexSpace->numPhases(); iphase++){
    this->computeDnPhi(a_dphidn[iphase], a_phibc[iphase], a_phi, iphase);
  }
}

void JumpBc::computeDnPhi(LevelData<BaseIVFAB<Real> >&       a_dphidn,
			  const LevelData<BaseIVFAB<Real> >& a_phibc,
			  const LevelData<MFCellFAB>&        a_phi,
			  const int                          a_phase){

  for (DataIterator dit = a_dphidn.dataIterator(); dit.ok(); ++dit){
    this->computeDnPhi(a_dphidn[dit()], a_phibc[dit()], a_phi[dit()].getPhase(a_phase), dit(), a_phase);
  }
}

void JumpBc::computeDnPhi(LevelData<BaseIVFAB<Real> >&       a_dphidn,
			  const LevelData<BaseIVFAB<Real> >& a_phibc,
			  const LevelData<EBCellFAB>&        a_phi,
			  const int                          a_phase){
  
  for (DataIterator dit = a_dphidn.dataIterator(); dit.ok(); ++dit){
    this->computeDnPhi(a_dphidn[dit()], a_phibc[dit()], a_phi[dit()], dit(), a_phase);
  }
}

void JumpBc::computeDnPhi(BaseIVFAB<Real>&       a_dphidn,
			  const BaseIVFAB<Real>& a_phibc,
			  const EBCellFAB&       a_phi,
			  const DataIndex&       a_dit,
			  const int              a_phase){
  const int comp = 0.;

  const BaseIVFAB<Real>& wb             = m_weights[a_dit].getIVFAB(a_phase);
  const BaseIVFAB<VoFStencil>& stencils = m_stencils[a_dit].getIVFAB(a_phase);

  for (VoFIterator vofit(m_ivs[a_dit], a_dphidn.getEBGraph()); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();

    // Apply stencil wb*phi_B + sum(wi*phi_i)
    Real dphidn = wb(vof, comp)*a_phibc(vof,comp);

    const VoFStencil& sten = stencils(vof, comp);
    for (int i = 0; i < sten.size(); i++){
      const VolIndex& ivof = sten.vof(i);
      const Real& weight   = sten.weight(i);
      dphidn += weight*a_phi(ivof, comp);
    }

    a_dphidn(vof, comp) = dphidn;
  }
}

LayoutData<MFInterfaceFAB<VoFStencil> >& JumpBc::getStencils(){
  return m_stencils;
}

LayoutData<MFInterfaceFAB<Real> >& JumpBc::getWeights(){
  return m_weights;
}

LayoutData<MFInterfaceFAB<Real> >& JumpBc::getAvgBco(){
  return m_bco;
}

LayoutData<MFInterfaceFAB<VoFStencil> >& JumpBc::getAvgStencils(){
  return m_avgStencils;
}

LayoutData<MFInterfaceFAB<Real> >& JumpBc::getAvgWeights(){
  return m_avgWeights;
}

LayoutData<MFInterfaceFAB<Real> >& JumpBc::get_avgBco(){
  return m_avgBco;
}

LayoutData<MFInterfaceFAB<Real> >& JumpBc::getInhomogeneous(){
  return m_inhomo;
}

LayoutData<MFInterfaceFAB<Real> >& JumpBc::getHomogeneous(){
  return m_homog;
}

#include <CD_NamespaceFooter.H>