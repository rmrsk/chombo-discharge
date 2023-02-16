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
#include <CD_IrregAddOp.H>
#include <CD_DataOps.H>
#include <CD_ParticleOps.H>
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

  CH_assert(a_meshData.getRealm() == a_particles.getRealm());

  DataOps::setValue(a_meshData, 0.0);

  // Deposit on this level
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    CH_assert(!(a_meshData[lvl].isNull()));
    CH_assert(a_meshData[lvl]->nComp() == 1);
    CH_assert(a_meshData[lvl]->ghostVect() >= m_radius * IntVect::Unit);

    const DisjointBoxLayout& dbl   = m_ebGrids[lvl]->getDBL();
    const EBISLayout&        ebisl = m_ebGrids[lvl]->getEBISL();

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box                    cellBox  = dbl[dit()];
      const EBISBox&               ebisbox  = ebisl[dit()];
      const BaseIVFAB<VoFStencil>& stencils = (*m_depositionStencils[lvl])[dit()];

      BaseIVFAB<Real>& meshData  = (*a_meshData[lvl])[dit()];
      const List<P>&   particles = a_particles[lvl][dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        const P& p = lit();

        const IntVect iv = ParticleOps::getParticleCellIndex(p.position(), m_probLo, m_dx[lvl]);

        // If these fail the particle is not in a valid cell.
        CH_assert(ebisbox.isIrregular(iv));

        const VoFStencil& stencil = stencils(VolIndex(iv, 0), 0);

        for (int i = 0; i < stencil.size(); i++) {
          const VolIndex& stencilVoF       = stencil.vof(i);
          const Real&     stencilWeight    = stencil.weight(i);
          const Real      depositionWeight = stencilWeight * (p.*particleScalarField)();

          meshData(stencilVoF, 0) += depositionWeight;
        }
      }
    }

    // Exchange ghost data. This adds the ghosted data into to the valid region.
    a_meshData[lvl]->exchange(Interval(0, 0), m_copierLevel[lvl], IrregAddOp());
  }

  // Monkey with coarse-fine interfaces.
  this->addInvalidCoarseDataToFineData(a_meshData);
  this->addFineGhostDataToValidCoarData(a_meshData);
}

template <class P, Real (P::*particleScalarField)()>
void
EBAMRSurfaceDeposition::deposit(EBAMRIVData& a_meshData, const ParticleContainer<P>& a_particles) const noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::deposit");
  if (m_verbose) {
    pout() << "EBAMRSurfaceDeposition::deposit<P, Real P::*func>" << endl;
  }

  MayDay::Error("EBAMRSurfaceDeposition::deposit -- this version is not supported (yet)");
}

#include <CD_NamespaceFooter.H>

#endif
