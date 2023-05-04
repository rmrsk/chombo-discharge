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
#include <EBCellFactory.H>
#include <EBFluxFactory.H>
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
  m_refRat   = a_refRat;

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

  Copier copier;
  copier.ghostDefine(dblCoFi, dbl, domain, IntVect::Unit);
  CH_STOP(t1);

  // Define the irregular coarse-fine interface. This is way more involved since we need to figure out which cells
  // on the coarse-grid side interfaces with cells on the fine-grid side.
  CH_START(t2);
  LevelData<FArrayBox> coarMask(dbl, 1, IntVect::Zero);
  LevelData<FArrayBox> coFiMask(dblCoFi, 1, IntVect::Unit);

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        coarMask[dit()].setVal(0.0);
      }

      for (DataIterator dit(dblCoFi); dit.ok(); ++dit) {
        coFiMask[dit()].setVal(0.0);

        const Box boxCoFi     = dblCoFi[dit()];
        const Box sideBoxCoFi = adjCellBox(boxCoFi, dir, sit(), 1);

        if (domain.contains(sideBoxCoFi)) {
          DenseIntVectSet cfivs = DenseIntVectSet(sideBoxCoFi, true);

          NeighborIterator nit(dblCoFi);
          for (nit.begin(dit()); nit.ok(); ++nit) {
            cfivs -= dblCoFi[nit()];
          }

          for (DenseIntVectSetIterator ivsIt(cfivs); ivsIt.ok(); ++ivsIt) {
            coFiMask[dit()](ivsIt(), 0) = 1.0;
          }
        }
      }

      // Copy the data to the coarse grid.
      const Interval interv(0, 0);
      coFiMask.copyTo(interv, coarMask, interv, copier, LDaddOp<FArrayBox>());

      // Define regular cells
      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        const FArrayBox& mask = coarMask[dit()];
        DenseIntVectSet  cfivs(dbl[dit()], false);

        for (BoxIterator bit(dbl[dit()]); bit.ok(); ++bit) {
          if (mask(bit(), 0) > 0.0) {
            cfivs |= bit();
          }
        }

        cfivs.recalcMinBox();
	m_regularCoarseFineRegions[dit()][std::make_pair(dir, sit())] = cfivs;
      }

      // Define irregular cells.
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
  CH_STOP(t2);
}

void
EBReflux::reflux(LevelData<EBCellFAB>&       a_Lphi,
                 const LevelData<EBFluxFAB>& a_flux,
                 const LevelData<EBFluxFAB>& a_fineFlux,
                 const Interval              a_variables,
                 const Real                  a_scaleCoarFlux,
                 const Real                  a_scaleFineFlux) const noexcept
{
  CH_TIMERS("EBReflux::reflux");
  CH_TIMER("EBReflux::reflux::define_buffers", t1);

  CH_assert(m_isDefined);
  CH_assert(a_Lphi.nComp() > a_variables.end());
  CH_assert(a_flux.nComp() > a_variables.end());
  CH_assert(a_fineFlux.nComp() > a_variables.end());

  const DisjointBoxLayout& dbl     = m_eblg.getDBL();
  const DisjointBoxLayout& dblCoFi = m_eblgCoFi.getDBL();

  const EBISLayout& ebisl     = m_eblg.getEBISL();
  const EBISLayout& ebislCoFi = m_eblgCoFi.getEBISL();

  CH_START(t1);
  LevelData<EBFluxFAB> fluxCoFi(dblCoFi, 1, IntVect::Unit, EBFluxFactory(ebislCoFi));
  LevelData<EBFluxFAB> fluxCoar(dbl, 1, IntVect::Unit, EBFluxFactory(ebisl));
  CH_STOP(t1);

  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {

    // Coarsen fluxes and copy to fluxCoar
    this->coarsenFluxesCF(fluxCoFi, a_fineFlux, 0, ivar);
    fluxCoFi.copyTo(fluxCoar);
    fluxCoar.exchange();

    // Reflux the coarse level.
    this->refluxIntoCoarse(a_Lphi, a_flux, fluxCoar, ivar, 0, ivar, a_scaleCoarFlux, a_scaleFineFlux);
  }
}

void
EBReflux::coarsenFluxesCF(LevelData<EBFluxFAB>&       a_coarFluxes,
                          const LevelData<EBFluxFAB>& a_fineFluxes,
                          const int                   a_coarVar,
                          const int                   a_fineVar) const noexcept
{
  CH_TIMERS("EBReflux::coarsenFluxes");
  CH_TIMER("EBReflux::coarsenFluxes::regular_faces", t1);
  CH_TIMER("EBReflux::coarsenFluxes::irregular_faces", t2);

  MayDay::Warning(
    "EBReflux::coarsenFluxesCF - since we only reflux on the CF, we can optimize this by only coarsening fluxes on the CF");

  CH_assert(m_isDefined);
  CH_assert(a_coarVar >= 0);
  CH_assert(a_fineVar >= 0);
  CH_assert(a_coarFluxes.nComp() > a_coarVar);
  CH_assert(a_fineFluxes.nComp() > a_fineVar);

  const DisjointBoxLayout& dblCoar = m_eblgCoFi.getDBL();
  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();

  const EBISLayout& ebislCoar = m_eblgCoFi.getEBISL();
  const EBISLayout& ebislFine = m_eblgFine.getEBISL();

  const Real dxCoar         = 1.0;
  const Real dxFine         = dxCoar / m_refRat;
  const Real invFinePerCoar = 1.0 / std::pow(m_refRat, SpaceDim - 1);

  for (DataIterator dit(dblCoar); dit.ok(); ++dit) {
    for (int dir = 0; dir < SpaceDim; dir++) {
      const Box coarCellBox = dblCoar[dit()];
      const Box fineCellBox = dblFine[dit()];
      const Box coarFaceBox = surroundingNodes(coarCellBox, dir);

      const EBISBox& coarEBISBox = ebislCoar[dit()];
      const EBISBox& fineEBISBox = ebislFine[dit()];

      const EBGraph& coarGraph = coarEBISBox.getEBGraph();
      const EBGraph& fineGraph = fineEBISBox.getEBGraph();

      EBFaceFAB&       coarFlux = a_coarFluxes[dit()][dir];
      const EBFaceFAB& fineFlux = a_fineFluxes[dit()][dir];

      FArrayBox&       coarFluxReg = coarFlux.getFArrayBox();
      const FArrayBox& fineFluxReg = fineFlux.getFArrayBox();

      // Some hooks for deciding how to do the averaging.
      const int xDoLoop = (dir == 0) ? 0 : 1;
      const int yDoLoop = (dir == 1) ? 0 : 1;
#if CH_SPACEDIM == 3
      const int zDoLoop = (dir == 2) ? 0 : 1;
#endif

      // Kernel for regular cells.
      auto regularKernel = [&](const IntVect& iv) -> void {
        coarFluxReg(iv, a_coarVar) = 0.0;

#if CH_SPACEDIM == 3
        for (int k = 0; k <= (m_refRat - 1) * zDoLoop; k++) {
#endif
          for (int j = 0; j <= (m_refRat - 1) * yDoLoop; j++) {
            for (int i = 0; i <= (m_refRat - 1) * xDoLoop; i++) {
              const IntVect ivFine = iv * m_refRat + IntVect(D_DECL(i, j, k));

              coarFluxReg(iv, a_coarVar) += fineFluxReg(ivFine, a_fineVar);
            }
          }
#if CH_SPACEDIM == 3
        }
#endif

        coarFluxReg(iv, a_coarVar) *= invFinePerCoar;
      };

      auto irregularKernel = [&](const FaceIndex& face) -> void {
        coarFlux(face, a_coarVar) = 0.0;

        const Real areaCoar = coarEBISBox.areaFrac(face);

        if (areaCoar > 0.0) {
          const Vector<FaceIndex>& fineFaces = ebislCoar.refine(face, m_refRat, dit());

          for (int i = 0; i < fineFaces.size(); i++) {
            const FaceIndex& fineFace = fineFaces[i];
            //            coarFlux(face, a_coarVar) += fineEBISBox.areaFrac(fineFace) * fineFlux(fineFace, a_fineVar);
            coarFlux(face, a_coarVar) += fineFlux(fineFace, a_fineVar);
          }

          coarFlux(face, a_coarVar) *= invFinePerCoar;
        }
      };

      const IntVectSet coarIrregIVS = coarEBISBox.getIrregIVS(coarCellBox);
      FaceIterator     faceIt(coarIrregIVS, coarGraph, dir, FaceStop::SurroundingWithBoundary);

      CH_START(t1);
      BoxLoops::loop(coarFaceBox, regularKernel);
      CH_STOP(t1);

      CH_START(t2);
      BoxLoops::loop(faceIt, irregularKernel);
      CH_STOP(t2);
    }
  }
}

void
EBReflux::refluxIntoCoarse(LevelData<EBCellFAB>&       a_Lphi,
                           const LevelData<EBFluxFAB>& a_oldFluxes,
                           const LevelData<EBFluxFAB>& a_newFluxes,
                           const int                   a_phiVar,
                           const int                   a_oldFluxVar,
                           const int                   a_newFluxVar,
                           const Real                  a_scaleCoarFlux,
                           const Real                  a_scaleFineFlux) const noexcept
{
  CH_TIMERS("EBReflux::refluxIntoCoarse");
  CH_TIMER("EBReflux::refluxIntoCoarse::regular_cells", t1);
  CH_TIMER("EBReflux::refluxIntoCoarse::irregular_cells", t2);

  CH_assert(m_isDefined);
  CH_assert(a_Lphi.nComp() > a_phivar);
  CH_assert(a_oldFluxes.nComp() > a_oldFluxVar);
  CH_assert(a_newFluxes.nComp() > a_newFluxVar);

  const DisjointBoxLayout& dbl   = m_eblg.getDBL();
  const EBISLayout&        ebisl = m_eblg.getEBISL();

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box      cellBox = dbl[dit()];
    const EBISBox& ebisBox = ebisl[dit()];

    EBCellFAB& Lphi    = a_Lphi[dit()];
    FArrayBox& LphiReg = Lphi.getFArrayBox();

    for (int dir = 0; dir < SpaceDim; dir++) {
      const EBFaceFAB& oldFlux = a_oldFluxes[dit()][dir];
      const EBFaceFAB& newFlux = a_newFluxes[dit()][dir];

      const FArrayBox& oldFluxReg = oldFlux.getFArrayBox();
      const FArrayBox& newFluxReg = newFlux.getFArrayBox();

      for (SideIterator sit; sit.ok(); ++sit) {
        const int     iHiLo = sign(flip(sit()));
        const IntVect shift = (sit() == Side::Lo) ? BASISV(dir) : IntVect::Zero;

        auto regularKernel = [&](const IntVect& iv) -> void {
          if (ebisBox.isRegular(iv)) {
            LphiReg(iv, a_phiVar) -= iHiLo * a_scaleCoarFlux * oldFluxReg(iv + shift, a_oldFluxVar);
            LphiReg(iv, a_phiVar) += iHiLo * a_scaleFineFlux * newFluxReg(iv + shift, a_newFluxVar);
          }
        };

        auto irregularKernel = [&](const VolIndex& vof) -> void {
          const Vector<FaceIndex>& faces = ebisBox.getFaces(vof, dir, flip(sit()));

          for (int iface = 0; iface < faces.size(); iface++) {
            const FaceIndex& face     = faces[iface];
            const Real       faceArea = ebisBox.areaFrac(face);

            Lphi(vof, a_phiVar) -= iHiLo * faceArea * a_scaleCoarFlux * oldFlux(face, a_oldFluxVar);
            Lphi(vof, a_phiVar) += iHiLo * faceArea * a_scaleFineFlux * newFlux(face, a_oldFluxVar);
          }
        };

        const DenseIntVectSet& regularCFIVS   = m_regularCoarseFineRegions[dit()].at(std::make_pair(dir, sit()));
        VoFIterator&           irregularCFIVS = m_irregularCoarseFineRegions[dit()].at(std::make_pair(dir, sit()));

        CH_START(t1);
        BoxLoops::loop(regularCFIVS, regularKernel);
        CH_STOP(t1);

        CH_START(t1);
        BoxLoops::loop(irregularCFIVS, irregularKernel);
        CH_STOP(t1);
      }
    }
  }
}

#include <CD_NamespaceFooter.H>
