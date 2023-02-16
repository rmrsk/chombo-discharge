/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBAMRSurfaceDeposition.cpp
  @brief  Implementation of CD_EBAMRSurfaceDeposition.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>
#include <BaseIVFactory.H>
#include <NeighborIterator.H>
#include <CH_Timer.H>

// Our includes
#include <CD_IrregAddOp.H>
#include <CD_EBAMRSurfaceDeposition.H>
#include <CD_NamespaceHeader.H>

EBAMRSurfaceDeposition::EBAMRSurfaceDeposition() noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::EBAMRSurfaceDeposition(weak)");

  m_debug   = false;
  m_verbose = false;
}

EBAMRSurfaceDeposition::EBAMRSurfaceDeposition(const Vector<RefCountedPtr<EBLevelGrid>>& a_ebGrids,
                                               const Vector<RefCountedPtr<EBLevelGrid>>& a_ebGridsCoarsenedFine,
                                               const Vector<RefCountedPtr<EBLevelGrid>>& a_ebGridsRefinedCoar,
                                               const Vector<int>&                        a_refRat,
                                               const Vector<Real>&                       a_dx,
                                               const RealVect&                           a_probLo,
                                               const int                                 a_finestLevel,
                                               const int                                 a_radius) noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::EBAMRSurfaceDeposition(full)");

  this->define(a_ebGrids,
               a_ebGridsCoarsenedFine,
               a_ebGridsRefinedCoar,
               a_refRat,
               a_dx,
               a_probLo,
               a_finestLevel,
               a_radius);
}

EBAMRSurfaceDeposition::~EBAMRSurfaceDeposition() noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::~EBAMRSurfaceDeposition");
}

void
EBAMRSurfaceDeposition::define(const Vector<RefCountedPtr<EBLevelGrid>>& a_ebGrids,
                               const Vector<RefCountedPtr<EBLevelGrid>>& a_ebGridsCoarsenedFine,
                               const Vector<RefCountedPtr<EBLevelGrid>>& a_ebGridsRefinedCoar,
                               const Vector<int>&                        a_refRat,
                               const Vector<Real>&                       a_dx,
                               const RealVect&                           a_probLo,
                               const int                                 a_finestLevel,
                               const int                                 a_radius) noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::define");
  if (m_verbose) {
    pout() << "EBAMRSurfaceDeposition::define" << endl;
  }

  m_debug                = false;
  m_verbose              = true;
  m_ebGrids              = a_ebGrids;
  m_ebGridsCoarsenedFine = a_ebGridsCoarsenedFine;
  m_ebGridsRefinedCoar   = a_ebGridsRefinedCoar;
  m_refRat               = a_refRat;
  m_dx                   = a_dx;
  m_probLo               = a_probLo;
  m_finestLevel          = a_finestLevel;

  CH_assert(m_finestLevel >= 0);
  CH_assert(m_radius >= 0);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const bool hasFine = lvl < m_finestLevel;
    const bool hasCoar = lvl > 0;

    CH_assert(!(a_ebGrids[lvl].isNull()));

    if (lvl < m_finestLevel) {
      CH_assert(a_refRat[lvl] >= 2);
      CH_assert(a_refRat[lvl] % 2 == 0);
    }

    if (hasCoar) {
      CH_assert(!(a_ebGridsRefinedCoar[lvl].isNull()));
    }
    if (hasFine) {
      CH_assert(!(a_ebGridsCoarsenedFine[lvl].isNull()));
    }
  }

  // Put in debug mode or not.
  ParmParse pp("EBAMRSurfaceDeposition");
  pp.query("debug", m_debug);
  pp.query("verbose", m_verbose);

  this->defineBuffers();
  this->defineDataMotion();
  this->defineDepositionStencils();

  m_isDefined = true;
}

inline void
EBAMRSurfaceDeposition::defineBuffers() noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::defineBuffers");
  if (m_verbose) {
    pout() << "EBAMRSurfaceDeposition::defineBuffers" << endl;
  }

  constexpr int nComp  = 1;
  constexpr int nGhost = 0;

  m_refinedCoarData.resize(1 + m_finestLevel);

  for (int lvl = 1; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& refinedCoarDBL   = m_ebGridsRefinedCoar[lvl]->getDBL();
    const EBISLayout&        refinedCoarEBISL = m_ebGridsRefinedCoar[lvl]->getEBISL();

    LayoutData<IntVectSet> irregCells(refinedCoarDBL);

    for (DataIterator dit(refinedCoarDBL); dit.ok(); ++dit) {
      const Box      box     = refinedCoarDBL[dit()];
      const EBISBox& ebisbox = refinedCoarEBISL[dit()];

      irregCells[dit()] = ebisbox.getIrregIVS(box);
    }

    m_refinedCoarData[lvl] = RefCountedPtr<LevelData<BaseIVFAB<Real>>>(
      new LevelData<BaseIVFAB<Real>>(refinedCoarDBL,
                                     nComp,
                                     nGhost * IntVect::Unit,
                                     BaseIVFactory<Real>(refinedCoarEBISL, irregCells)));
  }
}

inline void
EBAMRSurfaceDeposition::defineDataMotion() noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::defineDataMotion");
  if (m_verbose) {
    pout() << "EBAMRSurfaceDeposition::defineDataMotion" << endl;
  }

  m_copierLevel.resize(1 + m_finestLevel);
  m_copierRefinedCoarToFineNoGhosts.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const bool hasCoar = lvl > 0;
    const bool hasFine = lvl < m_finestLevel;

    const EBLevelGrid&       eblg   = *m_ebGrids[lvl];
    const ProblemDomain&     domain = eblg.getDomain();
    const DisjointBoxLayout& dbl    = eblg.getDBL();

    // Define Copier. Note that the below code is basically the same as ghostDefine().
    const bool doExchange = true;

    // Define Copier as going from valid -> valid+ghost.
    m_copierLevel[lvl].define(dbl, dbl, domain, m_radius * IntVect::Unit, doExchange);

    // Define Copier as going from valid+ghost -> valid.
    m_copierLevel[lvl].reverse();

    if (hasCoar) {
      const DisjointBoxLayout& refinedCoarDBL = m_ebGridsRefinedCoar[lvl]->getDBL();

      m_copierRefinedCoarToFineNoGhosts[lvl].define(refinedCoarDBL, dbl, domain);
    }
  }
}

inline void
EBAMRSurfaceDeposition::defineDepositionStencils() noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::defineDepositionStencils");
  if (m_verbose) {
    pout() << "EBAMRSurfaceDeposition::defineDepositionStencils" << endl;
  }

  constexpr int nComp = 1;

  m_depositionStencils.resize(1 + m_finestLevel);

  // Define over valid cut-cells (i.e., cut-cells not covered by a finer grid)
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = m_ebGrids[lvl]->getDBL();
    const EBISLayout&        ebisl = m_ebGrids[lvl]->getEBISL();

    m_depositionStencils[lvl] = RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil>>>(
      new LayoutData<BaseIVFAB<VoFStencil>>(dbl));

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box        box     = dbl[dit()];
      const EBISBox&   ebisbox = ebisl[dit()];
      const EBGraph&   ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs     = ebisbox.getIrregIVS(box);

      (*m_depositionStencils[lvl])[dit()].define(ivs, ebgraph, 1);
    }
  }

  // Compute deposition weights for current cells.
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl    = m_ebGrids[lvl]->getDBL();
    const EBISLayout&        ebisl  = m_ebGrids[lvl]->getEBISL();
    const ProblemDomain&     domain = m_ebGrids[lvl]->getDomain();

    const Real dx = m_dx[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box      box     = dbl[dit()];
      const EBISBox& ebisBox = ebisl[dit()];
      const EBGraph& ebGraph = ebisBox.getEBGraph();

      BaseIVFAB<VoFStencil>& stencils = (*m_depositionStencils[lvl])[dit()];

      // Build stencils.
      const IntVectSet& ivs     = stencils.getIVS();
      const EBGraph&    ebgraph = stencils.getEBGraph();

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit) {
        const VolIndex curVoF = vofit();
        const IntVect  curIV  = curVoF.gridIndex();

        VoFStencil& stencil = stencils(curVoF, 0);
        stencil.clear();

        Real             totalArea = 0.0;
        Vector<VolIndex> stencilVoFs;

        // These are the cells that we are going to deposit into.
        const IntVectSet neighborhood = ebisBox.getIrregIVS(grow(Box(curIV, curIV), m_radius));

        for (IVSIterator ivsIt(neighborhood); ivsIt.ok(); ++ivsIt) {
          const IntVect iv = ivsIt();

          const Vector<VolIndex> vofs = ebisBox.getVoFs(iv);

          for (int i = 0; i < vofs.size(); i++) {
            totalArea += ebisBox.bndryArea(vofs[i]) * std::pow(dx, SpaceDim - 1);
            stencilVoFs.push_back(vofs[i]);
          }
        }

        if (totalArea > std::numeric_limits<Real>::min()) {
          totalArea = 1.0 / totalArea;
        }
        else {
          totalArea = 0.0;
        }

        for (int i = 0; i < stencilVoFs.size(); i++) {
          if (domain.contains(stencilVoFs[i].gridIndex())) {
            stencil.add(stencilVoFs[i], totalArea);
          }
        }
      };
    }
  }
}

inline void
EBAMRSurfaceDeposition::addInvalidCoarseDataToFineData(EBAMRIVData& a_meshData) const noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::addInvalidCoarseDataToFineData");
  if (m_verbose) {
    pout() << "EBAMRSurfaceDeposition::addInvalidCoarseDataToFineData" << endl;
  }

  // TLDR: We want to ADD an interpolation of the coarse-grid data to the fine grid.
  for (int lvl = 1; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dblCoar = m_ebGrids[lvl - 1]->getDBL();

    const Real dxCoar   = m_dx[lvl - 1];
    const Real dxFine   = m_dx[lvl];
    const Real dxFactor = std::pow(dxCoar / dxFine, SpaceDim - 1);

    const EBISLayout& ebislFine = m_ebGridsRefinedCoar[lvl]->getEBISL();
    const EBISLayout& ebislCoar = m_ebGrids[lvl - 1]->getEBISL();

    for (DataIterator dit(dblCoar); dit.ok(); ++dit) {
      const Box coarBox = dblCoar[dit()];

      const EBISBox& ebisBoxCoar = ebislCoar[dit()];
      const EBISBox& ebisBoxFine = ebislFine[dit()];

      const EBGraph&   ebGraphCoar = ebisBoxCoar.getEBGraph();
      const IntVectSet irregCoar   = ebisBoxCoar.getIrregIVS(coarBox);

      BaseIVFAB<Real>&       fineData = (*m_refinedCoarData[lvl])[dit()];
      const BaseIVFAB<Real>& coarData = (*a_meshData[lvl - 1])[dit()];

      fineData.setVal(0.0);

      // Piecewise constant conservative interpolation of coarse-grid data to the fine grid.
      for (VoFIterator vofit(irregCoar, ebGraphCoar); vofit.ok(); ++vofit) {
        const VolIndex& coarVoF = vofit();

        const Vector<VolIndex> fineVoFs = ebislCoar.refine(coarVoF, m_refRat[lvl - 1], dit());

        Real fineArea = 0.0;
        for (int i = 0; i < fineVoFs.size(); i++) {
          fineArea += ebisBoxFine.bndryArea(fineVoFs[i]);
        }

        if (fineArea > std::numeric_limits<Real>::min()) {
          for (int i = 0; i < fineVoFs.size(); i++) {
            fineData(fineVoFs[i], 0) = 0.0;
          }
        }
        else {
          const Real coarArea = ebisBoxCoar.bndryArea(coarVoF);
          const Real fineVal  = coarArea * coarData(coarVoF, 0) * dxFactor / fineArea;

          for (int i = 0; i < fineVoFs.size(); i++) {
            fineData(fineVoFs[i], 0) = fineVal;
          }
        }
      }
    }

    // Add the interpolated data to the fine-grid data.
    const Interval variables = Interval(0, 0);
    m_refinedCoarData[lvl]->copyTo(variables,
                                   *a_meshData[lvl],
                                   variables,
                                   m_copierRefinedCoarToFineNoGhosts[lvl],
                                   IrregAddOp());
  }
}

inline void
EBAMRSurfaceDeposition::addFineGhostDataToValidCoarData(EBAMRIVData& a_meshData) const noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::addFineGhostDataToValidCoarData");
  if (m_verbose) {
    pout() << "EBAMRSurfaceDeposition::addFineGhostDataToValidCoarData" << endl;
  }
}

#include <CD_NamespaceFooter.H>
