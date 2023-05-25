/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBRedistribution.cpp
  @brief  Implementation of CD_EBRedistribution.cpp
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <EBCellFactory.H>
#include <NeighborIterator.H>

// Our includes
#include <CD_EBRedistribution.H>
#include <CD_VofUtils.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

EBRedistribution::EBRedistribution() noexcept
{
  CH_TIME("EBRedistribution::EBRedistribution(weak)");

  m_isDefined = false;
}

EBRedistribution::EBRedistribution(const EBLevelGrid& a_eblgCoar,
                                   const EBLevelGrid& a_eblgRefinedCoar,
                                   const EBLevelGrid& a_eblg,
                                   const EBLevelGrid& a_eblgCoarsenedFine,
                                   const EBLevelGrid& a_eblgFine,
                                   const int          a_refToCoar,
                                   const int          a_refToFine,
                                   const int          a_redistributionRadius) noexcept
{
  CH_TIME("EBRedistribution::EBRedistribution(full)");

  this->define(a_eblgCoar,
               a_eblgRefinedCoar,
               a_eblg,
               a_eblgCoarsenedFine,
               a_eblgFine,
               a_refToCoar,
               a_refToFine,
               a_redistributionRadius);
}

EBRedistribution::~EBRedistribution() noexcept {
  CH_TIME("EBRedistribution::~EBRedistribution");
}

void
EBRedistribution::define(const EBLevelGrid& a_eblgCoar,
                         const EBLevelGrid& a_eblgRefinedCoar,
                         const EBLevelGrid& a_eblg,
                         const EBLevelGrid& a_eblgCoarsenedFine,
                         const EBLevelGrid& a_eblgFine,
                         const int          a_refToCoar,
                         const int          a_refToFine,
                         const int          a_redistributionRadius) noexcept
{
  CH_TIME("EBRedistribution::define");

  CH_assert(a_eblg.isDefined());
  CH_assert(a_redistributionRadius >= 1);

  m_refToCoar = -1;
  m_refToFine = -1;

  m_hasCoar = false;
  m_hasFine = false;

  m_eblg         = a_eblg;
  m_redistRadius = a_redistributionRadius;

  if (a_eblgCoar.isDefined()) {
    CH_assert(a_refToCoar >= 2);
    CH_assert(a_refToCoar % 2 == 0);
    CH_assert(a_eblgRefinedCoar.isDefined());

    m_eblgRefinedCoar = a_eblgRefinedCoar;
    m_eblgCoar        = a_eblgCoar;

    m_refToCoar = a_refToCoar;
    m_hasCoar   = true;
  }

  if (a_eblgFine.isDefined()) {
    CH_assert(a_refToFine >= 2);
    CH_assert(a_refToFine % 2 == 0);
    CH_assert(a_eblgCoarsenedFine.isDefined());

    m_eblgCoarsenedFine = a_eblgCoarsenedFine;
    m_eblgFine          = a_eblgFine;

    m_refToFine = a_refToFine;
    m_hasFine   = true;
  }

  this->defineStencils();
  this->defineBuffers();

  m_isDefined = true;
}

void
EBRedistribution::defineStencils() noexcept
{
  CH_TIMERS("EBRedistribution::defineStencils");

  // TLDR: This is a bit involved since we need to define stencils on the valid cut-cells on this level. If there's an EBCF interface we
  //       need to know about it because we don't want to:
  //       1) Redistribute from this level into ghost cells on the other side of the coarse-fine interface. That mass should go on the coarse level.
  //       2) Redistribute from this level into regions covered by the finer grid. That mass should go on the fine level instead.

  const DisjointBoxLayout& dbl   = m_eblg.getDBL();
  const EBISLayout&        ebisl = m_eblg.getEBISL();

  // These are maps of the valid cells on this level, and cells that lie on the interface. The interfaceCells data is used to figure out
  // which cells we will redistribute to when we redistribute from a cut-cell and across the coarse-fine interface into a coarse-grid cell. The
  // validCells is used for 1) restricting which cells we redistribute from, and 2) which fine-grid cells redistribute to.
  LevelData<BaseFab<bool>> interfaceCellsLD;
  LevelData<BaseFab<bool>> validCellsLD;

  this->defineValidCells(validCellsLD);
  this->defineInterfaceCells(interfaceCellsLD);

  m_vofit.define(dbl);
  m_redistStencilsCoar.define(dbl);
  m_redistStencilsLevel.define(dbl);
  m_redistStencilsFine.define(dbl);

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box            box            = dbl[dit()];
    const EBISBox&       ebisbox        = ebisl[dit()];
    const EBGraph&       ebgraph        = ebisbox.getEBGraph();
    const IntVectSet     irregIVS       = ebisbox.getIrregIVS(box);
    const BaseFab<bool>& validCells     = validCellsLD[dit()];
    const BaseFab<bool>& interfaceCells = interfaceCellsLD[dit()];

    IntVectSet redistCells;
    for (IVSIterator ivsIt(irregIVS); ivsIt.ok(); ++ivsIt) {
      const IntVect& iv = ivsIt();
      if (validCells(iv)) {
        redistCells |= iv;
      }
    }

    VoFIterator&           vofit         = m_vofit[dit()];
    BaseIVFAB<VoFStencil>& stencilsCoar  = m_redistStencilsCoar[dit()];
    BaseIVFAB<VoFStencil>& stencilsLevel = m_redistStencilsLevel[dit()];
    BaseIVFAB<VoFStencil>& stencilsFine  = m_redistStencilsFine[dit()];

    vofit.define(redistCells, ebgraph);
    stencilsCoar.define(redistCells, ebgraph, 1);
    stencilsLevel.define(redistCells, ebgraph, 1);
    stencilsFine.define(redistCells, ebgraph, 1);

    for (vofit.reset(); vofit.ok(); ++vofit) {
      const VolIndex& vof = vofit();

      // Get the VoFs in the neighborhood of this VoF. When we redistribute, also permit self-redistribution back into this cell.
      const bool             includeSelf = true;
      const Vector<VolIndex> neighborVofs =
        VofUtils::getVofsInRadius(vof, ebisbox, m_redistRadius, VofUtils::Connectivity::SimplyConnected, includeSelf);
    }
  }
}

void
EBRedistribution::defineValidCells(LevelData<BaseFab<bool>>& a_validCells) const noexcept
{
  CH_TIME("EBRedistribution::defineValidCells");

  const DisjointBoxLayout& dbl = m_eblg.getDBL();

  a_validCells.define(dbl, 1, IntVect::Zero);
  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    a_validCells[dit()].setVal(true);
  }

  // If there's a finer level we need to figure out which cells on the fine level overlap with this level. Then we set those cells
  // to 'false'. If there's no finer level then all cells on this level are valid cells.
  if (m_hasFine) {
    const DisjointBoxLayout& dblCoFi = m_eblgCoarsenedFine.getDBL();

    // Create some data = 0 on the coarse grid and = 1 on the fine grid.
    LevelData<FArrayBox> data(dbl, 1, IntVect::Zero);
    LevelData<FArrayBox> dataCoFi(dblCoFi, 1, IntVect::Zero);

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      data[dit()].setVal(0.0);
    }
    for (DataIterator dit(dblCoFi); dit.ok(); ++dit) {
      dataCoFi[dit()].setVal(1.0);
    }

    dataCoFi.copyTo(data);

    // Go through the coarse grid and set cells to false wherever we find data > 0.0
    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box cellBox = dbl[dit()];

      BaseFab<bool>&   validCells = a_validCells[dit()];
      const FArrayBox& mask       = data[dit()];

      auto kernel = [&](const IntVect& iv) -> void {
        if (mask(iv) > 0.0) {
          validCells(iv) = false;
        }
      };

      BoxLoops::loop(cellBox, kernel);
    }
  }
}

void
EBRedistribution::defineInterfaceCells(LevelData<BaseFab<bool>>& a_interfaceCells) const noexcept
{
  CH_TIME("EBRedistribution::defineInterfaceCells");

  const DisjointBoxLayout& dbl = m_eblg.getDBL();
  a_interfaceCells.define(dbl, 1, m_redistRadius * IntVect::Unit);

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    a_interfaceCells[dit()].setVal(false);
  }

  // If there's a coarse grid then some of the cells we would redistribute to overlap with the valid region
  // on the coarse grid.
  if (m_hasCoar) {
    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box grownBox = grow(dbl[dit()], m_redistRadius);

      DenseIntVectSet divs(grownBox, true);

      NeighborIterator nit(dbl);
      for (nit.begin(dit()); nit.ok(); ++nit) {
        divs -= dbl[nit()];
      }
      divs -= dbl[dit()];

      // Flag the cells on the other side of the interface.
      BaseFab<bool>& interfaceCells = a_interfaceCells[dit()];
      for (DenseIntVectSetIterator ivsIt(divs); ivsIt.ok(); ++ivsIt) {
        interfaceCells(ivsIt()) = true;
      }
    }
  }
}

void
EBRedistribution::defineBuffers() noexcept
{
  CH_TIME("EBRedistribution::defineBuffers");
}

void
EBRedistribution::redistributeAMR(LevelData<EBCellFAB>*             a_phiCoar,
                                  LevelData<EBCellFAB>*             a_phi,
                                  LevelData<EBCellFAB>*             a_phiFine,
                                  const LevelData<BaseIVFAB<Real>>& a_deltaM,
                                  const Real                        a_scaleCoar,
                                  const Real                        a_scale,
                                  const Real                        a_scaleFine,
                                  const Interval&                   a_variables) const noexcept
{
  CH_TIME("EBRedistribution::redistributeAMR");

  CH_assert(m_isDefined);
  CH_assert(a_phi != nullptr);

  if (m_hasCoar) {
    CH_assert(a_phiCoar != nullptr);

    this->redistributeCoar(*a_phiCoar, a_deltaM, a_scaleCoar, a_variables);
  }

  this->redistributeLevel(*a_phi, a_deltaM, a_scale, a_variables);

  if (m_hasFine) {
    CH_assert(a_phiFine != nullptr);

    this->redistributeFine(*a_phiFine, a_deltaM, a_scaleFine, a_variables);
  }
}

void
EBRedistribution::redistributeCoar(LevelData<EBCellFAB>&             a_phiCoar,
                                   const LevelData<BaseIVFAB<Real>>& a_deltaM,
                                   const Real&                       a_scaleCoar,
                                   const Interval&                   a_variables) const noexcept
{
  CH_TIME("EBRedistributino::redistributeCoar");

  CH_assert(m_isDefined);
  CH_assert(a_phiCoar.isDefined());
  CH_assert(a_deltaM.isDefined());
  CH_assert(a_phiCoar.nComp() > a_variables.end());
  CH_assert(a_deltaM.nComp() > a_variables.end());
}
void
EBRedistribution::redistributeLevel(LevelData<EBCellFAB>&             a_phi,
                                    const LevelData<BaseIVFAB<Real>>& a_deltaM,
                                    const Real&                       a_scale,
                                    const Interval&                   a_variables) const noexcept
{
  CH_TIME("EBRedistributino::redistributeLevel");

  CH_assert(m_isDefined);
  CH_assert(a_phi.isDefined());
  CH_assert(a_deltaM.isDefined());
  CH_assert(a_phi.nComp() > a_variables.end());
  CH_assert(a_deltaM.nComp() > a_variables.end());
}

void
EBRedistribution::redistributeFine(LevelData<EBCellFAB>&             a_phiFine,
                                   const LevelData<BaseIVFAB<Real>>& a_deltaM,
                                   const Real&                       a_scaleFine,
                                   const Interval&                   a_variables) const noexcept
{
  CH_TIME("EBRedistributino::redistributeFine");

  CH_assert(m_isDefined);
  CH_assert(a_phiFine.isDefined());
  CH_assert(a_deltaM.isDefined());
  CH_assert(a_phiFine.nComp() > a_variables.end());
  CH_assert(a_deltaM.nComp() > a_variables.end());
}

#include <CD_NamespaceFooter.H>
