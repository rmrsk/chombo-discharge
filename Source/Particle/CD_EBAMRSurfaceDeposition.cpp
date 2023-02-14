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
#include <CH_Timer.H>

// Our includes
#include <CD_EBAMRSurfaceDeposition.H>
#include <CD_NamespaceHeader.H>

EBAMRSurfaceDeposition::EBAMRSurfaceDeposition() noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::EBAMRSurfaceDeposition(weak)");

  m_debug   = false;
  m_verbose = false;
}

EBAMRSurfaceDeposition::EBAMRSurfaceDeposition(const Vector<RefCountedPtr<EBLevelGrid>>&              a_eblgs,
                                               const Vector<RefCountedPtr<EBLevelGrid>>&              a_eblgsCoFi,
                                               const Vector<RefCountedPtr<EBLevelGrid>>&              a_eblgsFiCo,
                                               const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_validCells,
                                               const Vector<int>&                                     a_refRat,
                                               const Vector<Real>&                                    a_dx,
                                               const RealVect&                                        a_probLo,
                                               const IntVect&                                         a_ghost,
                                               const int                                              a_finestLevel,
                                               const int                                              a_radius) noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::EBAMRSurfaceDeposition(full)");

  this->define(a_eblgs,
               a_eblgsCoFi,
               a_eblgsFiCo,
               a_validCells,
               a_refRat,
               a_dx,
               a_probLo,
               a_ghost,
               a_finestLevel,
               a_radius);
}

EBAMRSurfaceDeposition::~EBAMRSurfaceDeposition() noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::~EBAMRSurfaceDeposition");
}

inline void
EBAMRSurfaceDeposition::define(const Vector<RefCountedPtr<EBLevelGrid>>&              a_eblgs,
                               const Vector<RefCountedPtr<EBLevelGrid>>&              a_eblgsCoFi,
                               const Vector<RefCountedPtr<EBLevelGrid>>&              a_eblgsFiCo,
                               const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_validCells,
                               const Vector<int>&                                     a_refRat,
                               const Vector<Real>&                                    a_dx,
                               const RealVect&                                        a_probLo,
                               const IntVect&                                         a_ghost,
                               const int                                              a_finestLevel,
                               const int                                              a_radius) noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::define");

  m_debug       = false;
  m_verbose     = false;
  m_eblgs       = a_eblgs;
  m_eblgsCoFi   = a_eblgsCoFi;
  m_eblgsFiCo   = a_eblgsFiCo;
  m_validCells  = a_validCells;
  m_refRat      = a_refRat;
  m_dx          = a_dx;
  m_probLo      = a_probLo;
  m_ghost       = a_ghost;
  m_finestLevel = a_finestLevel;

  // Put in debug mode or not.
  ParmParse pp("EBAMRSurfaceDeposition");
  pp.query("debug", m_debug);
  pp.query("verbose", m_verbose);

  m_isDefined = true;
}

#include <CD_NamespaceFooter.H>
