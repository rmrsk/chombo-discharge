/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBRedistribution.cpp
  @brief  Implementation of CD_EBRedistribution.cpp
  @author Robert Marskar
*/

#ifndef CD_EBRedistribution_H
#define CD_EBRedistribution_H

// Chombo includes
#include <CH_Timer.H>
#include <EBCellFactory.H>

// Our includes
#include <CD_EBRedistribution.H>
#include <CD_NamespaceHeader.H>

EBRedistribution::EBRedistribution() noexcept {
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

  if (a_eblgCoar.isDefined()) {
    CH_assert(a_refToCoar >= 2);
    CH_assert(a_refToCoar % 2 == 0);
    CH_assert(a_eblgRefinedCoar.isDefined());

    m_hasCoar = true;
  }

  if (a_eblgFine.isDefined()) {
    CH_assert(a_refToFine >= 2);
    CH_assert(a_refToFine % 2 == 0);
    CH_assert(a_eblgCoarsenedFine.isDefined());

    m_hasFine = true;
  }

  m_isDefined = true;
}

#include <CD_NamespaceFooter.H>
