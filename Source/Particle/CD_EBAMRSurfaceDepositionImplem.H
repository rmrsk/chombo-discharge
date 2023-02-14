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
#include <CH_Timer.H>

// Our includes
#include <CD_EBAMRSurfaceDeposition.H>
#include <CD_NamespaceHeader.H>

template <class P, const Real& (P::*particleScalarField)() const>
void
EBAMRSurfaceDeposition::deposit(EBAMRIVData& a_meshData, const ParticleContainer<P>& a_particles) const noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::deposit");
  if (m_verbose) {
    pout() << "EBAMRSurfaceDeposition::deposit<P, const Real& P::*func const>" << endl;
  }
}

template <class P, Real (P::*particleScalarField)()>
void
EBAMRSurfaceDeposition::deposit(EBAMRIVData& a_meshData, const ParticleContainer<P>& a_particles) const noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::deposit");
  if (m_verbose) {
    pout() << "EBAMRSurfaceDeposition::deposit<P, Real P::*func>" << endl;
  }
}

#include <CD_NamespaceFooter.H>

#endif
