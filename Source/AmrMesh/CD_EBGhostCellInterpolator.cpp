/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBGhostCellInterpolator.cpp
  @brief  Implementation of CD_EBGhostCellInterpolator.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <IntVectSet.H>
#include <EBCellFactory.H>
#include <EBAlias.H>
#include <NeighborIterator.H>

// Our includes
#include <CD_EBGhostCellInterpolator.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

EBGhostCellInterpolator::EBGhostCellInterpolator() noexcept
{
  CH_TIME("EBGhostCellInterpolator::EBGhostCellInterpolator(weak)");

  m_isDefined = false;
}

EBGhostCellInterpolator::EBGhostCellInterpolator(const EBLevelGrid& a_eblgFine,
                                                 const EBLevelGrid& a_eblgCoFi,
                                                 const EBLevelGrid& a_eblgCoar,
                                                 const IntVect&     a_ghostVector,
                                                 const int          a_refRat,
                                                 const int          a_ghostCF) noexcept
{
  CH_TIME("EBGhostCellInterpolator::EBGhostCellInterpolator(strong)");

  this->define(a_eblgFine, a_eblgCoFi, a_eblgCoar, a_ghostVector, a_refRat, a_ghostCF);
}

EBGhostCellInterpolator::~EBGhostCellInterpolator() noexcept
{
  CH_TIME("EBGhostCellInterpolator::~EBGhostCellInterpolator");
}

void
EBGhostCellInterpolator::define(const EBLevelGrid& a_eblgFine,
                                const EBLevelGrid& a_eblgCoFi,
                                const EBLevelGrid& a_eblgCoar,
                                const IntVect&     a_ghostVector,
                                const int          a_refRat,
                                const int          a_ghostCF) noexcept
{
  CH_TIME("EBGhostCellInterpolator::define()");

  CH_assert(a_ghostCF >= 0);
  CH_assert(a_refRat >= 2);

  m_eblgFine    = a_eblgFine;
  m_eblgCoFi    = a_eblgCoFi;
  m_eblgCoar    = a_eblgCoar;
  m_ghostVector = a_ghostVector;
  m_refRat      = a_refRat;
  m_ghostCF     = a_ghostCF;

  this->defineGhostRegions();

  m_isDefined = true;
}

void
EBGhostCellInterpolator::defineGhostRegions() noexcept
{
  CH_TIME("EBGhostCellInterpolator::defineGhostRegions");

  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();
  const DisjointBoxLayout& dblCoar = m_eblgCoFi.getDBL();

  const EBISLayout& ebislFine = m_eblgFine.getEBISL();
  const EBISLayout& ebislCoar = m_eblgCoFi.getEBISL();

  const ProblemDomain& domainFine = m_eblgFine.getDomain();
  const ProblemDomain& domainCoar = m_eblgCoFi.getDomain();

  m_regularGhostRegions.define(dblFine);
  m_fineIrregCells.define(dblFine);
  m_coarIrregCells.define(dblCoar);

  for (DataIterator dit(dblFine); dit.ok(); ++dit) {
    const Box fineCellBox = dblFine[dit()];
    const Box coarCellBox = dblCoar[dit()];

    const EBISBox& fineEBISBox = ebislFine[dit()];
    const EBISBox& coarEBISBox = ebislCoar[dit()];

    const EBGraph& fineGraph = fineEBISBox.getEBGraph();
    const EBGraph& coarGraph = coarEBISBox.getEBGraph();

    // Regular ghost cells.
    for (int dir = 0; dir < SpaceDim; dir++) {
      for (SideIterator sit; sit.ok(); ++sit) {

        // Compute the layer of ghost cells on the lo/hi side for this patch. Note that
        // we need to grow the ghosted box in directions orthogonal to 'dir' because we
        // also do the corner ghost cells.
        IntVectSet cfivs = IntVectSet(adjCellBox(fineCellBox, dir, sit(), m_ghostCF));
        for (int otherDir = 0; otherDir < SpaceDim; otherDir++) {
          if (otherDir != dir) {
            cfivs.grow(otherDir, m_ghostCF);
          }
        }
        cfivs &= domainFine;
        cfivs.recalcMinBox();

        const Box fineGhostBox = cfivs.minBox();
        const Box coarGhostBox = coarsen(fineGhostBox, m_refRat);

        const IntVectSet coarIrregIVS = coarEBISBox.getIrregIVS(coarGhostBox);
        const IntVectSet fineIrregIVS = refine(coarIrregIVS, m_refRat);

        // Define our objects.
        m_regularGhostRegions[dit()].emplace(std::make_pair(dir, sit()), fineGhostBox);
        m_fineIrregCells[dit()].emplace(std::make_pair(dir, sit()), VoFIterator(fineIrregIVS, fineGraph));
        m_coarIrregCells[dit()].emplace(std::make_pair(dir, sit()), VoFIterator(coarIrregIVS, coarGraph));
      }
    }
  }
}

void
EBGhostCellInterpolator::interpolate(LevelData<EBCellFAB>&       a_phiFine,
                                     const LevelData<EBCellFAB>& a_phiCoar,
                                     const Interval              a_variables,
                                     const Type                  a_interpType) const noexcept
{
  CH_TIME("EBGhostCellInterpolator::interpolate(LD<EBCellFAB>");

  CH_assert(a_phiFine.nComp() > a_variables.end());
  CH_assert(a_phiCoar.nComp() > a_variables.end());

  // Define buffer that we need. Need two ghost cells for doing the coarse-side slopes.
  const int numCoarGhostCells = 2;

  const DisjointBoxLayout& coFiGrids = m_eblgCoFi.getDBL();
  const EBISLayout&        coFiEBISL = m_eblgCoFi.getEBISL();

  LevelData<EBCellFAB> grownCoarData(coFiGrids, 1, numCoarGhostCells * IntVect::Unit, EBCellFactory(coFiEBISL));

  for (int icomp = a_variables.begin(); icomp <= a_variables.end(); icomp++) {
    const Interval srcInterv = Interval(icomp, icomp);
    const Interval dstInterv = Interval(0, 0);
    a_phiCoar.copyTo(srcInterv, grownCoarData, dstInterv);
    grownCoarData.exchange();

    // Fill invalid regions.
    for (DataIterator dit(m_eblgFine.getDBL()); dit.ok(); ++dit) {
      this->interpolateRegular(a_phiFine[dit()].getFArrayBox(),
                               grownCoarData[dit()].getFArrayBox(),
                               dit(),
                               icomp,
                               0,
                               a_interpType);
    }
  }

  // Fill valid regions.
  a_phiFine.exchange(a_variables);
}

void
EBGhostCellInterpolator::interpolateRegular(FArrayBox&       a_phiFine,
                                            const FArrayBox& a_phiCoar,
                                            const DataIndex& a_dit,
                                            const int        a_fineVar,
                                            const int        a_coarVar,
                                            const Type       a_interpType) const noexcept
{
  CH_TIMERS("EBGhostCellInterpolator::interpolateRegular(EBCellFAB)");
  CH_TIMER("EBGhostCellInterpolator::interpolateRegular(EBCellFAB)::set_fine_to_coar", t1);
  CH_TIMER("EBGhostCellInterpolator::interpolateRegular(EBCellFAB)::compute_slopes", t2);
  CH_TIMER("EBGhostCellInterpolator::interpolateRegular(EBCellFAB)::apply_slopes", t3);

  CH_assert(a_phiFine.nComp() > a_fineVar);
  CH_assert(a_phiCoar.nComp() > a_coarVar);

  const ProblemDomain& domainCoar = m_eblgCoFi.getDomain();

  // Storage for grid slopes
  FArrayBox slopes(a_phiCoar.box(), 1);

  // Interpolate on all coarse-fine interface sides.
  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      const Box interpBox = m_regularGhostRegions[a_dit].at(std::make_pair(dir, sit()));

      // Kernel for setting phiFine = phiCoar in the ghost cells.
      auto regSetFineToCoar = [&](const IntVect fineIV) -> void {
        const IntVect coarIV = coarsen(fineIV, m_refRat);

        a_phiFine(fineIV, a_fineVar) = a_phiCoar(coarIV, a_coarVar);
      };

      // Set phiFine = phiCoar
      CH_START(t1);
      BoxLoops::loop(interpBox, regSetFineToCoar);
      CH_STOP(t1);

      // Add contributions from slopes.
      const bool doSlopes = a_interpType != EBGhostCellInterpolator::Type::PWC &&
                            a_interpType != EBGhostCellInterpolator::Type::ConservativePWC;

      if (doSlopes) {
        for (int slopeDir = 0; slopeDir < SpaceDim; slopeDir++) {
          const IntVect s = BASISV(slopeDir);

          // Figure out how to compute the slopes.
          std::function<void(const IntVect& iv)> regularSlopeKernel;
          switch (a_interpType) {
          case EBGhostCellInterpolator::Type::MinMod: {
            regularSlopeKernel = [&](const IntVect& iv) -> void {
              Real dwl = 0.0;
              Real dwr = 0.0;

              if (domainCoar.contains(iv - s)) {
                dwl = a_phiCoar(iv, a_coarVar) - a_phiCoar(iv - s, a_coarVar);
              }
              if (domainCoar.contains(iv + s)) {
                dwr = a_phiCoar(iv + s, a_coarVar) - a_phiCoar(iv, a_coarVar);
              }

              slopes(iv, 0) = this->minmod(dwl, dwr);
            };

            break;
          }
          case EBGhostCellInterpolator::Type::MonotonizedCentral: {
            regularSlopeKernel = [&](const IntVect& iv) -> void {
              Real dwl = 0.0;
              Real dwr = 0.0;

              if (domainCoar.contains(iv - s)) {
                dwl = a_phiCoar(iv, a_coarVar) - a_phiCoar(iv - s, a_coarVar);
              }
              if (domainCoar.contains(iv + s)) {
                dwr = a_phiCoar(iv + s, a_coarVar) - a_phiCoar(iv, a_coarVar);
              }

              slopes(iv, 0) = this->monotonizedCentral(dwl, dwr);
            };

            break;
          }
          case EBGhostCellInterpolator::Type::Superbee: {
            regularSlopeKernel = [&](const IntVect& iv) -> void {
              Real dwl = 0.0;
              Real dwr = 0.0;

              if (domainCoar.contains(iv - s)) {
                dwl = a_phiCoar(iv, a_coarVar) - a_phiCoar(iv - s, a_coarVar);
              }
              if (domainCoar.contains(iv + s)) {
                dwr = a_phiCoar(iv + s, a_coarVar) - a_phiCoar(iv, a_coarVar);
              }

              slopes(iv, 0) = this->superbee(dwl, dwr);
            };

            break;
          }
          default: {
            MayDay::Error("EBGhostCellInterpolator::interpolate(EBCellFAB) - logic bust");
          }
          }

          // Kernel for adding in slope limiter.
          auto addRegularSlopeContribution = [&](const IntVect& fineIV) -> void {
            const IntVect  coarIV = coarsen(fineIV, m_refRat);
            const RealVect delta = (RealVect(fineIV) - m_refRat * RealVect(coarIV) + 0.5 * (1.0 - m_refRat)) / m_refRat;

            a_phiFine(fineIV, a_fineVar) += slopes(coarIV, 0) * delta[slopeDir];
          };

          // Compute slopes and add contributions into phiFine.
          const Box coarsenedInterpBox = coarsen(interpBox, m_refRat);

          CH_START(t2);
          BoxLoops::loop(coarsenedInterpBox, regularSlopeKernel);
          CH_STOP(t2);

          CH_START(t3);
          BoxLoops::loop(interpBox, addRegularSlopeContribution);
          CH_STOP(t3);
        }
      }
    }
  }
}

Real
EBGhostCellInterpolator::minmod(const Real& dwl, const Real& dwr) const noexcept
{
  Real slope = 0.0;

  if (dwl * dwr > 0.0) {
    slope = std::abs(dwl) < std::abs(dwr) ? dwl : dwr;
  }

  return slope;
}

Real
EBGhostCellInterpolator::superbee(const Real& dwl, const Real& dwr) const noexcept
{
  Real slope = 0.0;

  if (dwl * dwr > 0.0) {
    const Real s1 = this->minmod(dwl, 2 * dwr);
    const Real s2 = this->minmod(dwr, 2 * dwl);

    if (s1 * s2 > 0.0) {
      slope = std::abs(s1) > std::abs(s2) ? s1 : s2;
    }
  }

  return slope;
}

Real
EBGhostCellInterpolator::monotonizedCentral(const Real& dwl, const Real& dwr) const noexcept
{
  Real slope = 0.0;

  if (dwl * dwr > 0.0) {
    const Real dwc = dwl + dwr;
    const Real sgn = Real((dwc > 0.0) - (dwc < 0.0));

    slope = sgn * std::min(0.5 * std::abs(dwc), 2.0 * std::min(std::abs(dwl), std::abs(dwr)));
  }

  return slope;
}

#include <CD_NamespaceFooter.H>
