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
#include <CH_Timer.H>

// Our includes
#include <CD_EBReflux.H>
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
  CH_TIME("EBReflux::defineRegionsCF");

  const DisjointBoxLayout& dbl     = m_eblg.getDBL();
  const DisjointBoxLayout& dblCoFi = m_eblgCoFi.getDBL();

  const EBISLayout& ebisl     = m_eblg.getEBISL();
  const EBISLayout& ebislCoFi = m_eblgCoFi.getEBISL();

  m_regularCoarseFineRegions.define(dbl);
  m_irregularCoarseFineRegions.define(dbl);

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box boxCoar = dbl[dit()];

    // Define map.
    auto& regularCoarseFineRegions   = m_regularCoarseFineRegions[dit()];
    auto& irregularCoarseFineRegions = m_irregularCoarseFineRegions[dit()];

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
