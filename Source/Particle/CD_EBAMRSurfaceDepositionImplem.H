/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBAMRSurfaceDepositionImplem.H
  @brief  Implementation of CD_EBAMRSurfaceDeposition.H
  @author Robert Marskar
*/

#ifndef CD_EBAMRSurfaceDepositionImplem_H
#define CD_EBAMRSurfaceDepositionImplem_H

// Chombo includes
#include <ParmParse.H>
#include <CH_Timer.H>

// Our includes
#include <CD_EBAMRSurfaceDeposition.H>
#include <CD_NamespaceHeader.H>

EBAMRSurfaceDeposition::EBAMRSurfaceDeposition() noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::EBAMRSurfaceDeposition(weak)");

  m_debug = false;
}

EBAMRSurfaceDeposition::EBAMRSurfaceDeposition(const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgs,
                                               const Vector<int>&                        a_refRat,
                                               const Vector<Real>&                       a_dx,
                                               const RealVect&                           a_probLo,
                                               const IntVect&                            a_ghost,
                                               const int                                 a_finestLevel) noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::EBAMRSurfaceDeposition(full)");

  this->define(a_eblgs, a_refRat, a_dx, a_probLo, a_ghost, a_finestLevel);
}

EBAMRSurfaceDeposition::~EBAMRSurfaceDeposition() noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::~EBAMRSurfaceDeposition");
}

inline void
EBAMRSurfaceDeposition::define(const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgs,
                               const Vector<int>&                        a_refRat,
                               const Vector<Real>&                       a_dx,
                               const RealVect&                           a_probLo,
                               const IntVect&                            a_ghost,
                               const int                                 a_finestLevel) noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::define");

  m_debug       = false;
  m_eblgs       = a_eblgs;
  m_refRat      = a_refRat;
  m_dx          = a_dx;
  m_probLo      = a_probLo;
  m_ghost       = a_ghost;
  m_finestLevel = a_finestLevel;

  // Put in debug mode or not.
  ParmParse pp("EBAMRSurfaceDeposition");
  pp.query("debug", m_debug);
}

template <class P, const Real& (P::*particleScalarField)() const>
void
EBAMRSurfaceDeposition::deposit(EBAMRIVData& a_meshData, const ParticleContainer<P>& a_particles) const noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::deposit");
}

template <class P, Real (P::*particleScalarField)()>
void
EBAMRSurfaceDeposition::deposit(EBAMRIVData& a_meshData, const ParticleContainer<P>& a_particles) const noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::deposit");
}

#include <CD_NamespaceFooter.H>

#endif
