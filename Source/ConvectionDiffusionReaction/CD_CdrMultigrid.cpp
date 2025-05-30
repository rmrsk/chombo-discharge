/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrMultigrid.cpp
  @brief  Implementation of CD_CdrMultigrid.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>
#include <EBAMRIO.H>

// Our includes
#include <CD_CdrMultigrid.H>
#include <CD_DataOps.H>
#include <CD_EBHelmholtzNeumannDomainBCFactory.H>
#include <CD_EBHelmholtzDirichletDomainBCFactory.H>
#include <CD_EBHelmholtzNeumannEBBCFactory.H>
#include <CD_NamespaceHeader.H>

CdrMultigrid::CdrMultigrid() : CdrSolver()
{
  CH_TIME("CdrMultigrid::CdrMultigrid()");

  // Default settings
  m_name               = "CdrMultigrid";
  m_className          = "CdrMultigrid";
  m_hasMultigridSolver = false;
}

CdrMultigrid::~CdrMultigrid()
{}

void
CdrMultigrid::registerOperators()
{
  CH_TIME("CdrMultigrid::registerOperators()");
  if (m_verbosity > 5) {
    pout() << m_name + "::registerOperators()" << endl;
  }

  // Need to register everything that the base class registered, plus the amr interpolator
  CdrSolver::registerOperators();

  m_amr->registerOperator(s_eb_multigrid, m_realm, m_phase);
  m_amr->registerOperator(s_eb_flux_reg, m_realm, m_phase);
}

void
CdrMultigrid::allocate()
{
  CH_TIME("CdrMultigrid::allocate()");
  if (m_verbosity > 5) {
    pout() << m_name + "::allocate()" << endl;
  }

  CdrSolver::allocate();
}

void
CdrMultigrid::preRegrid(const int a_lbase, const int a_oldFinestLevel)
{
  CH_TIME("CdrMultigrid::preRegrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::preRegrid" << endl;
  }

  CdrSolver::preRegrid(a_lbase, a_oldFinestLevel);

  m_hasMultigridSolver = false;
}

void
CdrMultigrid::resetAlphaAndBeta(const Real a_alpha, const Real a_beta)
{
  CH_TIME("CdrMultigrid::resetAlphaAndBeta(Real, Real");
  if (m_verbosity > 5) {
    pout() << m_name + "::resetAlphaAndBeta(Real, Real)" << endl;
  }

  if (m_isDiffusive) {
    Vector<MGLevelOp<LevelData<EBCellFAB>>*> multigridOperators = m_multigridSolver->getAllOperators();

    for (int i = 0; i < multigridOperators.size(); i++) {
      TGAHelmOp<LevelData<EBCellFAB>>* helmholtzOperator = (TGAHelmOp<LevelData<EBCellFAB>>*)multigridOperators[i];

      helmholtzOperator->setAlphaAndBeta(a_alpha, a_beta);
    }
  }
}

void
CdrMultigrid::setMultigridSolverCoefficients()
{
  CH_TIME("CdrMultigrid::setMultigridSolverCoefficients()");
  if (m_verbosity > 5) {
    pout() << "CdrMultigrid::setMultigridSolverCoefficients" << endl;
  }

  if (!m_hasMultigridSolver) {
    MayDay::Error("CdrMultigrid::setMultigridSolverCoefficients -- must set up solver first!");
  }

  // Get the AMR operators and update the coefficients.
  Vector<AMRLevelOp<LevelData<EBCellFAB>>*>& operatorsAMR = m_multigridSolver->getAMROperators();

  CH_assert(operatorsAMR.size() == 1 + m_amr->getFinestLevel());

  // Set coefficients for the AMR levels.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    CH_assert(!(operatorsAMR[lvl] == nullptr));

    EBHelmholtzOp& op = static_cast<EBHelmholtzOp&>(*operatorsAMR[lvl]);

    op.setAcoAndBco(m_helmAcoef[lvl], m_faceCenteredDiffusionCoefficient[lvl], m_ebCenteredDiffusionCoefficient[lvl]);
  }

  // Get the deeper multigrid levels and coarsen onto that data as well. Strictly speaking, we don't
  // have to do this but it facilitates multigrid convergence and is therefore good practice. The operator
  // factory has routines for the coefficients that belong to the multigrid levels. The factory does not
  // have access to the operator, so we fetch those using AMRMultiGrid and call setAcoAndBco from there.
  // access to the operator.
  m_helmholtzOpFactory->coarsenCoefficientsMG();
  Vector<Vector<MGLevelOp<LevelData<EBCellFAB>>*>> operatorsMG = m_multigridSolver->getOperatorsMG();

  for (int amrLevel = 0; amrLevel < operatorsMG.size(); amrLevel++) {
    for (int mgLevel = 0; mgLevel < operatorsMG[amrLevel].size(); mgLevel++) {
      CH_assert(!(operatorsMG[amrLevel][mgLevel] == nullptr));

      EBHelmholtzOp& op = static_cast<EBHelmholtzOp&>(*operatorsMG[amrLevel][mgLevel]);

      op.setAcoAndBco(op.getAcoef(), op.getBcoef(), op.getBcoefIrreg());
    }
  }
}

void
CdrMultigrid::computeKappaLphi(EBAMRCellData& a_kappaLphi, const EBAMRCellData& a_phi)
{
  CH_TIME("CdrMultigrid::computeKappaLphi(EBAMRCellData, EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeKappaLphi(EBAMRCellData, EBAMRCellData)" << endl;
  }

  const int coarsestLevel = 0;
  const int finestLevel   = m_amr->getFinestLevel();

  Vector<LevelData<EBCellFAB>*> LphiPtr;
  Vector<LevelData<EBCellFAB>*> phiPtr;

  m_amr->alias(LphiPtr, a_kappaLphi);
  m_amr->alias(phiPtr, a_phi);

  // Need to reset to make sure we are computing kappa*div(D*grad(phi)) and not something like kappa*(phi - div(D*grad(phi))).
  this->resetAlphaAndBeta(0.0, 1.0);

  m_multigridSolver->computeAMROperator(LphiPtr, phiPtr, finestLevel, coarsestLevel, false);
}

void
CdrMultigrid::advanceEuler(EBAMRCellData&       a_newPhi,
                           const EBAMRCellData& a_oldPhi,
                           const EBAMRCellData& a_source,
                           const Real           a_dt)
{
  CH_TIME("CdrMultigrid::advanceEuler(EBAMRCellData, EBAMRCellData, EBAMRCellData, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::advanceEuler(EBAMRCellData, EBAMRCellData, EBAMRCellData, Real)" << endl;
  }

  if (m_isDiffusive) {
    // Allocate some data = 0 which we use for computing the residual.
    EBAMRCellData zero;
    EBAMRCellData scratch;

    m_amr->allocate(zero, m_realm, m_phase, 1);
    m_amr->allocate(scratch, m_realm, m_phase, 1);

    DataOps::setValue(zero, 0.0);

    // TLDR: Recall that the elliptic operator solves
    //
    //          kappa*L(phi) = kappa*rho
    //
    //       rather than L(phi) = rho. So this means that our right-hand side needs to be kappa-weighted before we pass this into multigrid. We assume that the user
    //       has provided a source-term that comes in is already weighted, but that the old solution comes in unweighted.

    // Set up multigrid again because the diffusion coefficients might have changed underneath us.
    if (!m_hasMultigridSolver) {
      this->setupDiffusionSolver();
    }
    else {
      this->setMultigridSolverCoefficients();
    }

    // Make the right-hand side for the Euler equation. Since we are solving
    //
    //     kappa * phi^(k+1) - kappa*dt*L(phi^(k+1)) = kappa * phi^k + kappa*dt*S
    //
    // we need to put the right-hand side somewhere. We use scratch for holding kappa*phi^k + kappa*dt*S. Note that we assume that S comes in weighted.
    DataOps::copy(scratch, a_oldPhi);
    DataOps::kappaScale(scratch);
    DataOps::incr(scratch, a_source, a_dt);

    // As above, the alpha and beta-coefficients for the Helmholtz operator need to be 1 and -a_dt. The kappas on the left-hand side
    // in the above equation are absorbed into the Helmholtz operator so we don't need to worry about those.
    this->resetAlphaAndBeta(1.0, -a_dt);

    // Aliasing because Chombo did not always use smart pointers.
    Vector<LevelData<EBCellFAB>*> newPhi;
    Vector<LevelData<EBCellFAB>*> eulerRHS;
    Vector<LevelData<EBCellFAB>*> resid;
    Vector<LevelData<EBCellFAB>*> zer;

    m_amr->alias(newPhi, a_newPhi);
    m_amr->alias(eulerRHS, scratch);
    m_amr->alias(resid, m_residual);
    m_amr->alias(zer, zero);

    const int coarsestLevel = 0;
    const int finestLevel   = m_amr->getFinestLevel();

    // Figure out how far away we are form a "converged" solution.
    const Real zeroResid = m_multigridSolver->computeAMRResidual(zer, eulerRHS, finestLevel, coarsestLevel);

    // Set the convergence metric.
    m_multigridSolver->m_convergenceMetric = zeroResid;

    // Always from previous solution.
    DataOps::copy(a_newPhi, a_oldPhi);

    // Do the multigrid solve.
    m_multigridSolver->solveNoInitResid(newPhi, resid, eulerRHS, finestLevel, coarsestLevel, false);
  }
  else {
    DataOps::copy(a_newPhi, a_oldPhi);
  }
}

void
CdrMultigrid::advanceCrankNicholson(EBAMRCellData&       a_newPhi,
                                    const EBAMRCellData& a_oldPhi,
                                    const EBAMRCellData& a_source,
                                    const Real           a_dt)
{
  CH_TIME("CdrMultigrid::advanceCrankNicholson(EBAMRCellData, EBAMRCellData, EBAMRCellData, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::advanceCrankNicholson(EBAMRCellData, EBAMRCellData, EBAMRCellData, Real)" << endl;
  }

  if (m_isDiffusive) {
    // Allocate some data = 0 which we use for computing the residual.
    EBAMRCellData zero;
    EBAMRCellData scratch;

    m_amr->allocate(zero, m_realm, m_phase, 1);
    m_amr->allocate(scratch, m_realm, m_phase, 1);

    DataOps::setValue(zero, 0.0);

    const int coarsestLevel = 0;
    const int finestLevel   = m_amr->getFinestLevel();

    // TLDR: Recall that the elliptic operator solves
    //
    //          kappa*L(phi) = kappa*rho
    //
    //       rather than L(phi) = rho. So this means that our right-hand side needs to be kappa-weighted before we pass this into multigrid. We assume that the user
    //       has provided a source-term that comes in is already weighted, but that the old solution comes in unweighted.

    // Set up multigrid again because the diffusion coefficients might have changed underneath us.
    if (!m_hasMultigridSolver) {
      this->setupDiffusionSolver();
    }
    else {
      this->setMultigridSolverCoefficients();
    }

    // Make the right-hand side for the Euler equation. Since we are solving
    //
    //     kappa * phi^(k+1) - 0.5 * kappa*dt*L(phi^(k+1)) = kappa * phi^k + 0.5 * kappa*dt*L(phi^k) + kappa*dt*S^(k+1/2)
    //
    // we need to put the right-hand side somewhere. We use scratch for holding the right-hand side. Note that S^(k+1/2) should come in weighted and the old solution
    // come in unweighted.

    // First, put kappa*phi^k in scratch.
    DataOps::copy(scratch, a_oldPhi);
    DataOps::kappaScale(scratch);

    // Add the source term, which by assumption comes in weighted by the volume fraction.
    DataOps::incr(scratch, a_source, a_dt);

    // Compute kappa*L(phi^k) and add it to the scratch data holder.
    EBAMRCellData kappaLphi;
    EBAMRCellData phiScratch;

    m_amr->allocate(kappaLphi, m_realm, m_phase, 1);
    m_amr->allocate(phiScratch, m_realm, m_phase, 1);

    DataOps::copy(phiScratch, a_oldPhi);

    this->computeKappaLphi(kappaLphi, phiScratch);

    // After this we have scratch = kappa * phi^k + 0.5 * kappa*dt*L(phi^k) + kappa*dt*S^(k+1/2)
    DataOps::incr(scratch, kappaLphi, 0.5 * a_dt);

    // From the equation above, the alpha and beta-coefficients for the Helmholtz operator need to be 1 and -0.5*a_dt. The kappas on the left-hand side
    // in the above equation are absorbed into the Helmholtz operator so we don't need to worry about those.
    this->resetAlphaAndBeta(1.0, -0.5 * a_dt);

    // Aliasing because Chombo did not always use smart pointers.
    Vector<LevelData<EBCellFAB>*> newPhi;
    Vector<LevelData<EBCellFAB>*> eulerRHS;
    Vector<LevelData<EBCellFAB>*> resid;
    Vector<LevelData<EBCellFAB>*> zer;

    m_amr->alias(newPhi, a_newPhi);
    m_amr->alias(eulerRHS, scratch);
    m_amr->alias(resid, m_residual);
    m_amr->alias(zer, zero);

    // Figure out how far away we are form a "converged" solution.
    const Real zeroResid = m_multigridSolver->computeAMRResidual(zer, eulerRHS, finestLevel, coarsestLevel);

    // Set the convergence metric.
    m_multigridSolver->m_convergenceMetric = zeroResid;

    // Always from previous solution.
    DataOps::copy(a_newPhi, a_oldPhi);

    // Do the multigrid solve.
    m_multigridSolver->solveNoInitResid(newPhi, resid, eulerRHS, finestLevel, coarsestLevel, false);
  }
  else {
    DataOps::copy(a_newPhi, a_oldPhi);
  }
}

void
CdrMultigrid::setupDiffusionSolver()
{
  CH_TIME("CdrMultigrid::setupDiffusionSolver()");
  if (m_verbosity > 5) {
    pout() << m_name + "::setupDiffusionSolver()" << endl;
  }

  // This is storage which is needed if we are doing an implicit diffusion solve. I know that not all
  // diffusion solves are implicit, but this is really the easiest way of
  if (m_isDiffusive) {
    m_amr->allocate(m_helmAcoef, m_realm, m_phase, m_nComp);
    m_amr->allocate(m_residual, m_realm, m_phase, m_nComp);

    DataOps::setValue(m_helmAcoef, 1.0);
    DataOps::setValue(m_residual, 0.0);

    // This sets up the multigrid Helmholtz solver.
    this->setupHelmholtzFactory();
    this->setupMultigrid();

    m_hasMultigridSolver = true;
  }
}

void
CdrMultigrid::setupHelmholtzFactory()
{
  CH_TIME("CdrMultigrid::setupHelmholtzFactory()");
  if (m_verbosity > 5) {
    pout() << m_name + "::setupHelmholtzFactory()" << endl;
  }

  const Vector<RefCountedPtr<EBLevelGrid>>&             levelGrids   = m_amr->getEBLevelGrid(m_realm, m_phase);
  const Vector<RefCountedPtr<EBCoarAve>>&               coarAve      = m_amr->getCoarseAverage(m_realm, m_phase);
  const Vector<RefCountedPtr<EBReflux>>&                fluxReg      = m_amr->getFluxRegister(m_realm, m_phase);
  const Vector<RefCountedPtr<EBMultigridInterpolator>>& interpolator = m_amr->getMultigridInterpolator(m_realm,
                                                                                                       m_phase);

  // Coarsest domain used for multigrid. The user specifies the minimum number of cells in any
  // coordinate direction, and we coarsen until we have a domain which satisfies that constraint.
  ProblemDomain bottomDomain = m_amr->getDomains()[0];
  while (bottomDomain.domainBox().shortside() >= 2 * m_minCellsBottom) {
    bottomDomain.coarsen(2);
  }

  // Number of ghost cells in data holders
  const IntVect ghostPhi = m_amr->getNumberOfGhostCells() * IntVect::Unit;
  const IntVect ghostRhs = m_amr->getNumberOfGhostCells() * IntVect::Unit;

  // Handle to boundary condition factories for domain and EB in an EBHelmholtzOp context.
  auto domainBcFactory = RefCountedPtr<EBHelmholtzDomainBCFactory>(new EBHelmholtzNeumannDomainBCFactory(0.0));
  auto ebbcFactory     = RefCountedPtr<EBHelmholtzEBBCFactory>(new EBHelmholtzNeumannEBBCFactory(0.0));

  // Temp alpha and beta. Diffusion solvers will reset these later.
  const Real alpha = 1.0;
  const Real beta  = 1.0;

  // We need to have the Helmholtz A-coefficient=1 because we will end up solving equations like
  //
  //    phi^(k+1) - dt*Div(D*Grad(phi^(k+1)) = phi^k + dt*rho
  //
  // and in that case we have alpha * A = 1 and beta = -dt. EBHelmholtzOpFactory shouldn't be doing anything
  // with this data but diffusion solvers might change the alpha and beta under us. We set the A-coefficient to one.
  DataOps::setValue(m_helmAcoef, 1.0);

  // Set up the operator
  m_helmholtzOpFactory = RefCountedPtr<EBHelmholtzOpFactory>(
    new EBHelmholtzOpFactory(Location::Cell::Center,
                             alpha,
                             beta,
                             m_amr->getProbLo(),
                             levelGrids,
                             m_amr->getValidCells(m_realm),
                             interpolator,
                             fluxReg,
                             coarAve,
                             m_amr->getRefinementRatios(),
                             m_amr->getDx(),
                             m_helmAcoef.getData(),
                             m_faceCenteredDiffusionCoefficient.getData(),
                             m_ebCenteredDiffusionCoefficient.getData(),
                             domainBcFactory,
                             ebbcFactory,
                             ghostPhi,
                             ghostRhs,
                             m_smoother,
                             bottomDomain,
                             m_amr->getMaxBoxSize()));
}

void
CdrMultigrid::setupMultigrid()
{
  CH_TIME("CdrMultigrid::setupMultigrid()");
  if (m_verbosity > 5) {
    pout() << m_name + "::setupMultigrid()" << endl;
  }

  // Select the bottom solver
  LinearSolver<LevelData<EBCellFAB>>* botsolver = NULL;
  switch (m_bottomSolverType) {
  case BottomSolverType::Simple: {
    botsolver = &m_simpleSolver;

    break;
  }
  case BottomSolverType::BiCGStab: {
    botsolver = &m_bicgstab;

    break;
  }
  case BottomSolverType::GMRES: {
    botsolver = &m_gmres;

    m_gmres.m_verbosity = 0; // Shut up.
  }
  default: {
    MayDay::Error("CdrMultigrid::setupMultigrid() - logic bust in bottom solver setup");

    break;
  }
  }

  // Make m_multigridType into an int for multigrid
  int gmgType;
  switch (m_multigridType) {
  case MultigridType::VCycle: {
    gmgType = 1;

    break;
  }
  case MultigridType::WCycle: {
    gmgType = 2;

    break;
  }
  default: {
    MayDay::Error("CdrMultigrid::setupMultigrid() -- logic bust in multigrid type selection");

    break;
  }
  }

  const int           finestLevel    = m_amr->getFinestLevel();
  const ProblemDomain coarsestDomain = m_amr->getDomains()[0];

  // Define AMRMultiGrid
  m_multigridSolver = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB>>>(new AMRMultiGrid<LevelData<EBCellFAB>>());
  m_multigridSolver->define(coarsestDomain, *m_helmholtzOpFactory, botsolver, 1 + finestLevel);
  m_multigridSolver->setSolverParameters(m_multigridPreSmooth,
                                         m_multigridPostSmooth,
                                         m_multigridBottomSmooth,
                                         gmgType,
                                         m_multigridMaxIterations,
                                         m_multigridExitTolerance,
                                         m_multigridExitHang,
                                         1.E-99); // Residue set through other means

  m_multigridSolver->m_imin      = m_multigridMinIterations;
  m_multigridSolver->m_verbosity = m_multigridVerbosity;

  // Dummies for init
  EBAMRCellData dummy1;
  EBAMRCellData dummy2;

  m_amr->allocate(dummy1, m_realm, m_phase, m_nComp);
  m_amr->allocate(dummy2, m_realm, m_phase, m_nComp);

  DataOps::setValue(dummy1, 0.0);
  DataOps::setValue(dummy2, 0.0);

  // Aliasing
  Vector<LevelData<EBCellFAB>*> phi;
  Vector<LevelData<EBCellFAB>*> rhs;

  m_amr->alias(phi, dummy1);
  m_amr->alias(rhs, dummy2);

  // Init solver. This instantiates all the operators in AMRMultiGrid so we can just call "solve"
  m_multigridSolver->init(phi, rhs, finestLevel, 0);
}

void
CdrMultigrid::computeDivJ(EBAMRCellData& a_divJ,
                          EBAMRCellData& a_phi,
                          const Real     a_extrapDt,
                          const bool     a_conservativeOnly,
                          const bool     a_ebFlux,
                          const bool     a_domainFlux)
{
  CH_TIME("CdrMultigrid::computeDivJ(EBAMRCellData, EBAMRCelLData, Real, bool, bool, bool)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDivJ(EBAMRCellData, EBAMRCelLData, Real, bool, bool, bool)" << endl;
  }

  if (m_isMobile || m_isDiffusive) {
    EBAMRFluxData scratchFlux;
    m_amr->allocate(scratchFlux, m_realm, m_phase, m_nComp);

    // Fill ghost cells
    m_amr->interpGhostPwl(a_phi, m_realm, m_phase);

    // We will let scratchFlux hold the total flux = advection + diffusion fluxes
    DataOps::setValue(scratchFlux, 0.0);

    if (m_isMobile && !m_isDiffusive) {
      m_amr->interpGhostPwl(m_cellVelocity, m_realm, m_phase);

      // Update face velocity and advect to faces.
      this->averageVelocityToFaces();
      this->advectToFaces(m_faceStates, a_phi, a_extrapDt);
      this->computeAdvectionFlux(scratchFlux, m_faceVelocity, m_faceStates, a_domainFlux);
    }
    else if (!m_isMobile && m_isDiffusive) {
      this->computeDiffusionFlux(scratchFlux, a_phi, a_domainFlux); // Domain flux needs to come in through here.
      DataOps::scale(scratchFlux, -1.0);
    }
    else if (m_isMobile && m_isDiffusive) { //
      m_amr->interpGhostPwl(m_cellVelocity, m_realm, m_phase);
      this->averageVelocityToFaces();
      this->advectToFaces(m_faceStates, a_phi, a_extrapDt);

      this->computeAdvectionDiffusionFlux(scratchFlux,
                                          a_phi,
                                          m_faceStates,
                                          m_faceVelocity,
                                          m_faceCenteredDiffusionCoefficient,
                                          a_domainFlux);
    }

    // General divergence computation. Also inject charge. Domain fluxes came in above but eb fluxes come in here.
    EBAMRIVData* ebflux;
    if (a_ebFlux) {
      ebflux = &m_ebFlux;
    }
    else {
      ebflux = &m_ebZero;
    }
    this->computeDivG(a_divJ, scratchFlux, *ebflux, a_conservativeOnly);
  }
  else {
    DataOps::setValue(a_divJ, 0.0);
  }

  return;
}

void
CdrMultigrid::computeDivF(EBAMRCellData& a_divF,
                          EBAMRCellData& a_phi,
                          const Real     a_extrapDt,
                          const bool     a_conservativeOnly,
                          const bool     a_ebFlux,
                          const bool     a_domainFlux)
{
  CH_TIME("CdrMultigrid::computeDivF(EBAMRCellData, EBAMRCellData, Real, bool, bool, bool)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDivF(EBAMRCellData, EBAMRCellData, Real, bool, bool, bool)" << endl;
  }

  if (m_isMobile) {
    EBAMRFluxData scratchFlux;
    m_amr->allocate(scratchFlux, m_realm, m_phase, m_nComp);

    // Fill ghost cells
    m_amr->interpGhostPwl(a_phi, m_realm, m_phase);
    m_amr->interpGhostPwl(m_cellVelocity, m_realm, m_phase);

    // Cell-centered velocities become face-centered velocities.
    this->averageVelocityToFaces();

    // Face extrapolation to cell-centered faces
    this->advectToFaces(m_faceStates, a_phi, a_extrapDt);

    // Compute face-centered fluxes
    this->computeAdvectionFlux(scratchFlux, m_faceVelocity, m_faceStates, a_domainFlux);

    EBAMRIVData* ebflux;
    if (a_ebFlux) {
      ebflux = &m_ebFlux;
    }
    else {
      ebflux = &m_ebZero;
    }

    // Compute div(F) -- this includes interpolation to centroids and redistribution.
    this->computeDivG(a_divF, scratchFlux, *ebflux, a_conservativeOnly);
  }
  else {
    DataOps::setValue(a_divF, 0.0);
  }
}

void
CdrMultigrid::computeDivD(EBAMRCellData& a_divD,
                          EBAMRCellData& a_phi,
                          const bool     a_conservativeOnly,
                          const bool     a_ebFlux,
                          const bool     a_domainFlux)
{
  CH_TIME("CdrMultigrid::computeDivD(EBAMRCellData, EBAMRCellData, bool, bool)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDivD(EBAMRCellData, EBAMRCellData, bool, bool)" << endl;
  }

  if (m_isDiffusive) {
    EBAMRFluxData scratchFlux;
    m_amr->allocate(scratchFlux, m_realm, m_phase, m_nComp);

    // Fill ghost cells
    m_amr->interpGhostPwl(a_phi, m_realm, m_phase);

    this->computeDiffusionFlux(scratchFlux, a_phi, a_domainFlux); // Compute the face-centered diffusion flux

    EBAMRIVData* ebflux;
    if (a_ebFlux) {
      ebflux = &m_ebFlux;
    }
    else {
      ebflux = &m_ebZero;
    }
    this->computeDivG(a_divD,
                      scratchFlux,
                      *ebflux,
                      a_conservativeOnly); // General face-centered flux to divergence magic.

    m_amr->conservativeAverage(a_divD, m_realm, m_phase);
    m_amr->interpGhost(a_divD, m_realm, m_phase);
  }
  else {
    DataOps::setValue(a_divD, 0.0);
  }
}

void
CdrMultigrid::parseMultigridSettings()
{
  CH_TIME("CdrMultigrid::parseMultigridSettings()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseMultigridSettings()" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;

  pp.get("gmg_verbosity", m_multigridVerbosity);
  pp.get("gmg_pre_smooth", m_multigridPreSmooth);
  pp.get("gmg_post_smooth", m_multigridPostSmooth);
  pp.get("gmg_bott_smooth", m_multigridBottomSmooth);
  pp.get("gmg_max_iter", m_multigridMaxIterations);
  pp.get("gmg_min_iter", m_multigridMinIterations);
  pp.get("gmg_exit_tol", m_multigridExitTolerance);
  pp.get("gmg_exit_hang", m_multigridExitHang);
  pp.get("gmg_min_cells", m_minCellsBottom);

  // Fetch the desired bottom solver from the input script. We look for things like CdrMultigrid.gmg_bottom_solver = bicgstab or '= simple <number>'
  // where <number> is the number of relaxation for the smoothing solver.
  const int num = pp.countval("gmg_bottom_solver");
  if (num == 1) {
    pp.get("gmg_bottom_solver", str);
    if (str == "bicgstab") {
      m_bottomSolverType = BottomSolverType::BiCGStab;
    }
    else if (str == "gmres") {
      m_bottomSolverType = BottomSolverType::GMRES;
    }
    else {
      MayDay::Error(
        "CdrMultigrid::parseMultigridSettings - logic bust, you've specified one parameter and I expected either 'bicgstab' or 'gmres'");
    }
  }
  else if (num == 2) {
    int numSmooth;
    pp.get("gmg_bottom_solver", str, 0);
    pp.get("gmg_bottom_solver", numSmooth, 1);
    if (str == "simple") {
      m_bottomSolverType = BottomSolverType::Simple;
      m_simpleSolver.setNumSmooths(numSmooth);
    }
    else {
      MayDay::Error(
        "CdrMultigrid::parseMultigridSettings - logic bust, you've specified two parameters and I expected 'simple <number>'");
    }
  }
  else {
    MayDay::Error(
      "CdrMultigrid::parseMultigridSettings - logic bust in bottom solver. You must specify ' = bicgstab', ' = gmres', or ' = simple <number>'");
  }

  // Relaxation type
  pp.get("gmg_smoother", str);
  if (str == "jacobi") {
    m_smoother = EBHelmholtzOp::Smoother::PointJacobi;
  }
  else if (str == "red_black") {
    m_smoother = EBHelmholtzOp::Smoother::GauSaiRedBlack;
  }
  else if (str == "multi_color") {
    m_smoother = EBHelmholtzOp::Smoother::GauSaiMultiColor;
  }
  else {
    MayDay::Error("CdrMultigrid::parseMultigridSettings - unknown relaxation method requested");
  }

  // Cycle type
  pp.get("gmg_cycle", str);
  if (str == "vcycle") {
    m_multigridType = MultigridType::VCycle;
  }
  else {
    MayDay::Error("CdrMultigrid::parseMultigridSettings - unknown cycle type requested");
  }

  // No lower than 2.
  if (m_minCellsBottom < 2) {
    m_minCellsBottom = 2;
  }
}

#include <CD_NamespaceFooter.H>
