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

EBAMRSurfaceDeposition::EBAMRSurfaceDeposition(const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgs,
                                               const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgsCoFi,
                                               const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgsFiCo,
                                               const Vector<int>&                        a_refRat,
                                               const Vector<Real>&                       a_dx,
                                               const RealVect&                           a_probLo,
                                               const int                                 a_finestLevel,
                                               const int                                 a_radius) noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::EBAMRSurfaceDeposition(full)");

  this->define(a_eblgs, a_eblgsCoFi, a_eblgsFiCo, a_refRat, a_dx, a_probLo, a_finestLevel, a_radius);
}

EBAMRSurfaceDeposition::~EBAMRSurfaceDeposition() noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::~EBAMRSurfaceDeposition");
}

inline void
EBAMRSurfaceDeposition::define(const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgs,
                               const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgsCoFi,
                               const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgsFiCo,
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

  m_debug       = false;
  m_verbose     = false;
  m_eblgs       = a_eblgs;
  m_eblgsCoFi   = a_eblgsCoFi;
  m_eblgsFiCo   = a_eblgsFiCo;
  m_refRat      = a_refRat;
  m_dx          = a_dx;
  m_probLo      = a_probLo;
  m_finestLevel = a_finestLevel;

  CH_assert(m_finestLevel >= 0);
  CH_assert(m_radius >= 0);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const bool hasFine = lvl < m_finestLevel;
    const bool hasCoar = lvl > m_finestLevel;

    CH_assert(!(a_eblgs[lvl].isNull()));

    if (lvl < m_finestLevel) {
      CH_assert(a_refRat[lvl] >= 2);
      CH_assert(a_refRat[lvl] % 2 == 0);
    }

    if (hasFine) {
      CH_assert(!(a_eblgsCoFi[lvl].isNull()));
    }
    if (hasCoar) {
      CH_assert(!(a_eblgsFiCo[lvl].isNull()));
    }
  }

  // Put in debug mode or not.
  ParmParse pp("EBAMRSurfaceDeposition");
  pp.query("debug", m_debug);
  pp.query("verbose", m_verbose);

  this->defineLevelMotion();
  this->defineDepositionStencils();

  m_isDefined = true;
}

inline void
EBAMRSurfaceDeposition::defineLevelMotion() noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::defineLevelMotion");
  if (m_verbose) {
    pout() << "EBAMRSurfaceDeposition::defineLevelMotion" << endl;
  }

  m_levelCopiers.resize(1 + m_finestLevel);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const EBLevelGrid&       eblg   = *m_eblgs[lvl];
    const ProblemDomain&     domain = eblg.getDomain();
    const DisjointBoxLayout& dbl    = eblg.getDBL();

    // Define Copier. Note that the below code is basically the same as ghostDefine().
    const bool doExchange = true;

    // Define Copier as going from valid -> valid+ghost.
    m_levelCopiers[lvl].define(dbl, dbl, domain, m_radius * IntVect::Unit, doExchange);

    // Define Copier as going from valid+ghost -> valid.
    m_levelCopiers[lvl].reverse();
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
    const DisjointBoxLayout& dbl   = m_eblgs[lvl]->getDBL();
    const EBISLayout&        ebisl = m_eblgs[lvl]->getEBISL();

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
    const DisjointBoxLayout& dbl   = m_eblgs[lvl]->getDBL();
    const EBISLayout&        ebisl = m_eblgs[lvl]->getEBISL();
    const Real               dx    = m_dx[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box      box     = dbl[dit()];
      const EBISBox& ebisBox = ebisl[dit()];
      const EBGraph& ebGraph = ebisBox.getEBGraph();

      BaseIVFAB<VoFStencil>& stencils = (*m_depositionStencils[lvl])[dit()];

      // Build a local representation of the coarse-fine interface around this patch.
      DenseIntVectSet  cfivs(grow(box, m_radius), 1);
      NeighborIterator nit(dbl);
      for (nit.begin(dit()); nit.ok(); ++nit) {
        cfivs -= dbl[nit()];
      }
      cfivs -= box;

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
          stencil.add(stencilVoFs[i], totalArea);
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
