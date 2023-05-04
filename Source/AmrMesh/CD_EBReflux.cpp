/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBReflux.cpp
  @brief  Implementation of CD_EBReflux.H
  @author Robert Marskar
*/

// Chombo includes
#include <NeighborIterator.H>
#include <CH_Timer.H>

// Our includes
#include <CD_EBReflux.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

EBReflux::EBReflux() noexcept
{
  CH_TIME("EBReflux::EBReflux(weak)");

  m_isDefined = false;
}

EBReflux::EBReflux(const EBLevelGrid& a_eblg,
                   const EBLevelGrid& a_eblgFine,
                   const EBLevelGrid& a_eblgCoFi,
                   const int          a_refRat) noexcept
{
  CH_TIME("EBReflux::EBReflux(full)");

  this->define(a_eblg, a_eblgFine, a_eblgCoFi, a_refRat);
}

EBReflux::~EBReflux() noexcept { CH_TIME("EBReflux::~EBReflux"); }

void
EBReflux::define(const EBLevelGrid& a_eblg,
                 const EBLevelGrid& a_eblgFine,
                 const EBLevelGrid& a_eblgCoFi,
                 const int          a_refRat) noexcept
{
  CH_TIME("EBReflux::define");

  m_eblg     = a_eblg;
  m_eblgFine = a_eblgFine;
  m_eblgCoFi = a_eblgCoFi;

  this->defineRegionsCF();

  m_isDefined = true;
}

void
EBReflux::defineRegionsCF() noexcept
{
  CH_TIMERS("EBReflux::defineRegionsCF");
  CH_TIMER("EBReflux::defineRegionsCF::layout_define", t1);
  CH_TIMER("EBReflux::defineRegionsCF::regular_regions", t2);
  CH_TIMER("EBReflux::defineRegionsCF::irregular_regions", t3);

  const DisjointBoxLayout& dbl     = m_eblg.getDBL();
  const DisjointBoxLayout& dblCoFi = m_eblgCoFi.getDBL();

  const ProblemDomain& domain     = m_eblg.getDomain();
  const ProblemDomain& domainCoFi = m_eblgCoFi.getDomain();

  const EBISLayout& ebisl     = m_eblg.getEBISL();
  const EBISLayout& ebislCoFi = m_eblgCoFi.getEBISL();

  CH_START(t1);
  m_regularCoarseFineRegions.define(dbl);
  m_irregularCoarseFineRegions.define(dbl);
  CH_STOP(t1);

  // Define the "regular" coarse-fine interface cells on the coarse side.
  CH_START(t2);
  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box boxCoar = dbl[dit()];

    auto& regularCoarseFineRegions = m_regularCoarseFineRegions[dit()];

    // I don't like this loop but hopefully it _shouldn't_ become a performance bottleneck. The issue is that
    // it will iterate through all fine-grid boxes and there can be many of them. I _could_ do this in conjunction
    // with how we define the irregular cells on the coarse-fine interface (using the masks), but I really want
    // to describe the regions regions using a collection boxes rather than an IntVectSet/VoFIterator.
    for (LayoutIterator lit = dblCoFi.layoutIterator(); lit.ok(); ++lit) {
      const Box boxCoFi = dblCoFi[lit()];

      // Check if the coarsened fine box intersects this box in any direction.
      for (int dir = 0; dir < SpaceDim; dir++) {
        for (SideIterator sit; sit.ok(); ++sit) {
          const Box sideBoxCoFi = adjCellBox(boxCoFi, dir, sit(), 1);

          if (boxCoar.intersectsNotEmpty(sideBoxCoFi)) {
            regularCoarseFineRegions[std::make_pair(dir, sit())].emplace_back(boxCoar & sideBoxCoFi);
          }
        }
      }
    }
  }
  CH_STOP(t2);

  // Define the irregular coarse-fine interface. This is way more involved since we need to figure out which cells
  // on the coarse-grid side interfaces with cells on the fine-grid side.
  CH_START(t3);
  LevelData<FArrayBox> coarMask(dbl, 1, IntVect::Zero);
  LevelData<FArrayBox> coFiMask(dblCoFi, 1, IntVect::Unit);

  Copier copier;
  copier.ghostDefine(dblCoFi, dbl, domain, IntVect::Unit);

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        coarMask[dit()].setVal(0.0);
      }

      for (DataIterator dit(dblCoFi); dit.ok(); ++dit) {
        coFiMask[dit()].setVal(0.0);

        const Box boxCoFi     = dblCoFi[dit()];
        const Box sideBoxCoFi = adjCellBox(boxCoFi, dir, sit(), 1);

        DenseIntVectSet cfivs = DenseIntVectSet(sideBoxCoFi, true);

        NeighborIterator nit(dblCoFi);
        for (nit.begin(dit()); nit.ok(); ++nit) {
          cfivs -= dblCoFi[nit()];
        }

        for (DenseIntVectSetIterator ivsIt(cfivs); ivsIt.ok(); ++ivsIt) {
          coFiMask[dit()](ivsIt(), 0) = 1.0;
        }
      }

      // Copy the data to the coarse grid.
      const Interval interv(0, 0);
      coFiMask.copyTo(interv, coarMask, interv, copier, LDaddOp<FArrayBox>());

      // Now define the VoFIterator on this side.
      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        const Box        cellBox = dbl[dit()];
        const FArrayBox& mask    = coarMask[dit()];
        const EBISBox&   ebisBox = ebisl[dit()];
        const EBGraph&   ebGraph = ebisBox.getEBGraph();

        IntVectSet irregCells;

        auto findIrregCells = [&](const IntVect& iv) -> void {
          if (mask(iv, 0) > 0.0 && ebisBox.isIrregular(iv)) {
            irregCells |= iv;
          }
        };

        BoxLoops::loop(cellBox, findIrregCells);

        // Define appropriate iterators.
        auto& irregularCoarseFineRegions = m_irregularCoarseFineRegions[dit()];
        irregularCoarseFineRegions[std::make_pair(dir, sit())].define(irregCells, ebGraph);
      }
    }
  }
  CH_STOP(t3);
}

void
EBReflux::reflux(LevelData<EBCellFAB>&       a_Lphi,
                 const LevelData<EBFluxFAB>& a_coarFlux,
                 const LevelData<EBFluxFAB>& a_fineFlux,
                 const Interval              a_variables) const noexcept
{
  CH_TIME("EBReflux::reflux");
}

#include <CD_NamespaceFooter.H>
