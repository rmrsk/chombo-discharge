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

  m_regularGhostCells.define(dblFine);
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

        // Compute the layer of ghost cells on the lo/hi side for this pathc.
        IntVectSet cfivs = IntVectSet(adjCellBox(fineCellBox, dir, sit(), m_ghostCF));

        NeighborIterator nit(dblFine);
        for (nit.begin(dit()); nit.ok(); ++nit) {
          cfivs -= dblFine[nit()];
        }

        cfivs &= domainFine;
        cfivs.recalcMinBox();

        const Box ghostBox = cfivs.minBox();

        // Coarsen the ghosted box and figure out if it contains cut-cells. We need to be careful when computing the coarse-grid slopes
        // in these cells.
        const Box        coarsenedGhostBox = coarsen(ghostBox, m_refRat);
        const IntVectSet coarIrregIVS      = coarEBISBox.getIrregIVS(coarsenedGhostBox);
        const IntVectSet refCoarIrregIVS   = refine(coarIrregIVS, m_refRat);

        // Define our objects.
        m_regularGhostCells[dit()].emplace(std::make_pair(dir, sit()), ghostBox);
        m_fineIrregCells[dit()].emplace(std::make_pair(dir, sit()), VoFIterator(refCoarIrregIVS, fineGraph));
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

  for (int fineVar = a_variables.begin(); fineVar <= a_variables.end(); fineVar++) {
    const int coarVar = 0;

    const Interval srcInterv = Interval(fineVar, fineVar);
    const Interval dstInterv = Interval(coarVar, coarVar);
    a_phiCoar.copyTo(srcInterv, grownCoarData, dstInterv);
    grownCoarData.exchange();

    const DisjointBoxLayout& dbl = m_eblgFine.getDBL();

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      this->interpolate(a_phiFine[dit()], grownCoarData[dit()], dit(), fineVar, coarVar, a_interpType);
    }
  }

  a_phiFine.exchange(a_variables);
}

void
EBGhostCellInterpolator::interpolate(EBCellFAB&       a_phiFine,
                                     const EBCellFAB& a_phiCoar,
                                     const DataIndex& a_dit,
                                     const int        a_fineVar,
                                     const int        a_coarVar,
                                     const Type       a_interpType) const noexcept
{
  CH_TIMERS("EBGhostCellInterpolator::interpolate(EBCellFAB)");
  CH_TIMER("EBGhostCellInterpolator::interpolate(EBCellFAB)::set_fine_to_coar", t1);
  CH_TIMER("EBGhostCellInterpolator::interpolate(EBCellFAB)::compute_slopes", t2);
  CH_TIMER("EBGhostCellInterpolator::interpolate(EBCellFAB)::apply_slopes", t3);

  CH_assert(a_phiFine.nComp() > a_fineVar);
  CH_assert(a_phiCoar.nComp() > a_coarVar);

  const EBISLayout&    ebisl  = m_eblgCoFi.getEBISL();
  const ProblemDomain& domain = m_eblgCoFi.getDomain();

  const EBISBox& ebisBox = ebisl[a_dit];

  EBCellFAB slopes(ebisBox, a_phiCoar.getRegion(), 1);

  FArrayBox&       phiFineReg = a_phiFine.getFArrayBox();
  FArrayBox&       slopesReg  = slopes.getFArrayBox();
  const FArrayBox& phiCoarReg = a_phiCoar.getFArrayBox();

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      const Box interpBox = m_regularGhostCells[a_dit].at(std::make_pair(dir, sit()));
      const Box coarseBox = coarsen(interpBox, m_refRat);

      // Set phiFine = phiCoar
      auto regSetFineToCoar = [&](const IntVect fineIV) -> void {
        const IntVect coarIV = coarsen(fineIV, m_refRat);

        phiFineReg(fineIV, a_fineVar) = phiCoarReg(coarIV, a_coarVar);
      };

      CH_START(t1);
      BoxLoops::loop(interpBox, regSetFineToCoar);
      CH_STOP(t1);

      // Add contributions from slopes.
      for (int slopeDir = 0; slopeDir < SpaceDim; slopeDir++) {
        const IntVect s = BASISV(slopeDir);

        std::function<void(const IntVect& iv)>   regularSlopeKernel;
        std::function<void(const VolIndex& vof)> irregularSlopeKernel;

        switch (a_interpType) {
        case EBGhostCellInterpolator::Type::PWC: {
          regularSlopeKernel = [&](const IntVect& iv) -> void {
            slopesReg(iv, 0) = 0.0;
          };

          irregularSlopeKernel = [&](const VolIndex& vof) -> void {

          };

          break;
        }
        case EBGhostCellInterpolator::Type::MinMod: {
          regularSlopeKernel = [&](const IntVect& iv) -> void {
            Real dwl = 0.0;
            Real dwr = 0.0;

            if (domain.contains(iv - s)) {
              dwl = phiCoarReg(iv, a_coarVar) - phiCoarReg(iv - s, a_coarVar);
            }
            if (domain.contains(iv + s)) {
              dwr = phiCoarReg(iv + s, a_coarVar) - phiCoarReg(iv, a_coarVar);
            }

            slopesReg(iv, 0) = this->minmod(dwl, dwr);
          };

          irregularSlopeKernel = [&](const VolIndex& vof) -> void {
            slopes(vof, 0) = 0.0;
          };

          break;
        }
        case EBGhostCellInterpolator::Type::MonotonizedCentral: {
          regularSlopeKernel = [&](const IntVect& iv) -> void {
            Real dwl = 0.0;
            Real dwr = 0.0;

            if (domain.contains(iv - s)) {
              dwl = phiCoarReg(iv, a_coarVar) - phiCoarReg(iv - s, a_coarVar);
            }
            if (domain.contains(iv + s)) {
              dwr = phiCoarReg(iv + s, a_coarVar) - phiCoarReg(iv, a_coarVar);
            }

            slopesReg(iv, 0) = this->monotonizedCentral(dwl, dwr);
          };

          irregularSlopeKernel = [&](const VolIndex& vof) -> void {
            slopes(vof, 0) = 0.0;
          };

          break;
        }
        case EBGhostCellInterpolator::Type::Superbee: {
          regularSlopeKernel = [&](const IntVect& iv) -> void {
            Real dwl = 0.0;
            Real dwr = 0.0;

            if (domain.contains(iv - s)) {
              dwl = phiCoarReg(iv, a_coarVar) - phiCoarReg(iv - s, a_coarVar);
            }
            if (domain.contains(iv + s)) {
              dwr = phiCoarReg(iv + s, a_coarVar) - phiCoarReg(iv, a_coarVar);
            }

            slopesReg(iv, 0) = this->superbee(dwl, dwr);
          };

          irregularSlopeKernel = [&](const VolIndex& vof) -> void {
            slopes(vof, 0) = 0.0;
          };

          break;
        }
        default: {
          MayDay::Error("EBGhostCellInterpolator::interpolate(EBCellFAB) - logic bust");
        }
        }

        auto addRegularSlopeContribution = [&](const IntVect& fineIV) -> void {
          const IntVect  coarIV = coarsen(fineIV, m_refRat);
          const RealVect delta  = (RealVect(fineIV) - m_refRat * RealVect(coarIV) + 0.5 * (1.0 - m_refRat)) / m_refRat;

          phiFineReg(fineIV, a_fineVar) += slopesReg(coarIV, 0) * delta[slopeDir];
        };

        // CH_START(t2);
        // BoxLoops::loop(coarseBox, regularSlopeKernel);
        // CH_STOP(t2);

        // CH_START(t3);
        // BoxLoops::loop(interpBox, addRegularSlopeContribution);
        // CH_STOP(t3);
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
