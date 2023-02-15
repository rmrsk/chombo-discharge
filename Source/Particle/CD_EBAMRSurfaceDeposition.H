/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBAMRSurfaceDeposition.H
  @brief  Declaration of a class for handling surface deposition of particles with EB and AMR
  @author Robert Marskar
*/

#ifndef CD_EBAMRSurfaceDeposition_H
#define CD_EBAMRSurfaceDeposition_H

// Chombo includes
#include <RefCountedPtr.H>
#include <Copier.H>
#include <EBLevelGrid.H>

// Our includes
#include <CD_EBAMRData.H>
#include <CD_ParticleContainer.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief class for handling surface deposition of particles with EB and AMR.
*/
class EBAMRSurfaceDeposition
{
public:
  /*!
    @brief Default constructor. Leaves the object in undefined state and requires the user to call the define function
  */
  EBAMRSurfaceDeposition() noexcept;

  /*!
    @brief Full constructor. Calls the define function.
    @param[in] a_eblgs            Grids
    @param[in] a_eblgsCoFi        Coarsened fine grids
    @param[in] a_eblgsFiCo        Refined coarse grids
    @param[in] a_refRat           Refinement factory between levels
    @param[in] a_dx               Grid resolutions.
    @πaram[in] a_probLo           Lower-left corner of physical domain.
    @param[in] a_finestLevel      Finest grid level.
    @param[in] a_radius           Deposition radius
  */
  EBAMRSurfaceDeposition(const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgs,
                         const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgsCoFi,
                         const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgsFiCo,
                         const Vector<int>&                        a_refRat,
                         const Vector<Real>&                       a_dx,
                         const RealVect&                           a_probLo,
                         const int                                 a_finestLevel,
                         const int                                 a_radius) noexcept;

  /*!
    @brief Destructor
  */
  virtual ~EBAMRSurfaceDeposition() noexcept;

  /*!
    @brief Define function. Puts objects in usable state. 
    @param[in] a_eblgs            Grids
    @param[in] a_eblgsCoFi        Coarsened fine grids
    @param[in] a_eblgsFiCo        Refined coarse grids
    @param[in] a_refRat           Refinement factory between levels
    @param[in] a_dx               Grid resolutions.
    @πaram[in] a_probLo           Lower-left corner of physical domain.
    @param[in] a_finestLevel      Finest grid level.
    @param[in] a_radius           Deposition radius
  */
  void
  define(const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgs,
         const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgsCoFi,
         const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgsFiCo,
         const Vector<int>&                        a_refRat,
         const Vector<Real>&                       a_dx,
         const RealVect&                           a_probLo,
         const int                                 a_finestLevel,
         const int                                 a_radius) noexcept;

  /*!
    @brief Deposit function. Deposits particle on surface.
    @details P is the particle type
    @note This will throw an error if a particle is located on a grid cell that is NOT an irregular cell
  */
  template <class P, const Real& (P::*particleScalarField)() const>
  void
  deposit(EBAMRIVData& a_meshData, const ParticleContainer<P>& a_particles) const noexcept;

  /*!
    @brief Deposit function. Deposits particle on surface.
    @details P is the particle type
    @note This will throw an error if a particle is located on a grid cell that is NOT an irregular cell
  */
  template <class P, Real (P::*particleScalarField)()>
  void
  deposit(EBAMRIVData& a_meshData, const ParticleContainer<P>& a_particles) const noexcept;

protected:
  /*!
    @brief Debug or not
  */
  bool m_debug;

  /*!
    @brief Is defined or not.
  */
  bool m_isDefined;

  /*!
    @brief Turn on/off chattiness
  */
  bool m_verbose;

  /*!
    @brief Lower-left corner of computational domain
  */
  RealVect m_probLo;

  /*!
    @brief Finest grid level
  */
  int m_finestLevel;

  /*!
    @brief Deposition radius
  */
  int m_radius;

  /*!
    @brief Grids on each level
  */
  Vector<RefCountedPtr<EBLevelGrid>> m_eblgs;

  /*!
    @brief Coarsened grids. 
    @note m_eblgsCoFi[i] contains the coarsening of m_eblgs[i] and m_eblgsCoFi[0] contains the nullpointer
  */
  Vector<RefCountedPtr<EBLevelGrid>> m_eblgsCoFi;

  /*!
    @brief Refined grids. 
    @note m_eblgsCoFi[i] contains the refinement of m_eblgs[i] and m_eblgsCoFi[m_finestLevel-1] contains the nullpointer
  */
  Vector<RefCountedPtr<EBLevelGrid>> m_eblgsFiCo;

  /*!
    @brief Refinement ratios between levels
  */
  Vector<int> m_refRat;

  /*!
    @brief Grid resolutions.
  */
  Vector<Real> m_dx;

  /*!
    @brief Deposition weights for each irregular cell. This contains 1.0/sum(alpha * dx^(D-1)) and the cells 
  */
  Vector<RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil>>>> m_depositionStencils;

  /*!
    @brief Level copier for adding data from valid+ghost -> valid. 
    @details Needed when adding ghost-deposited data back into valid cells
  */
  mutable Vector<Copier> m_levelCopiers;

  /*!
    @brief Define copiers
  */
  inline void
  defineLevelMotion() noexcept;

  /*!
   @brief Define deposition weights
  */
  inline void
  defineDepositionStencils() noexcept;

  /*!
    @brief Add the invalid coarse data (i.e., underneath the fine grid) to the fine level. 
    @details This does conservative interpolation to the fine grid
    @param[inout] a_meshData EB surface data
  */
  inline void
  addInvalidCoarseDataToFineData(EBAMRIVData& a_meshData) const noexcept;

  /*!
    @brief Add the ghosted fine-level data to the coarse data. 
    @details This does conservative coarsening of the fine-grid ghosted data to the coarse grid.
    @param[inout] a_meshData EB surface data
  */
  inline void
  addFineGhostDataToValidCoarData(EBAMRIVData& a_meshData) const noexcept;
};

#include <CD_NamespaceFooter.H>

#include <CD_EBAMRSurfaceDepositionImplem.H>

#endif
