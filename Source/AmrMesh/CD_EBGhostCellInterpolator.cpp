/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBGhostCellInterpolator.cpp
  @brief  Implementation of CD_EBGhostCellInterpolator.H
  @author Robert Marskar
*/

#include <CH_Timer.H>

// Our includes
#include <CD_EBGhostCellInterpolator.H>
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

  MayDay::Error("EBGhostCellInterpolator::define - not yet implemented");
}

void
EBGhostCellInterpolator::interpolate(LevelData<EBCellFAB>&       a_phiFine,
                                     const LevelData<EBCellFAB>& a_phiCoar,
                                     const Interval              a_variables,
                                     const Type                  a_interpType) const noexcept
{
  CH_TIME("EBGhostCellInterpolator::interpolate");

  MayDay::Error("EBGhostCellInterpolator::interpolate - not yet implemented");
}

void
EBGhostCellInterpolator::defineBuffer() const noexcept
{
  CH_TIME("EBGhostCellInterpolator::defineBuffer");
}

void
EBGhostCellInterpolator::undefineBuffer() const noexcept
{
  CH_TIME("EBGhostCellInterpolator::undefineBuffer");
}

#include <CD_NamespaceFooter.H>
