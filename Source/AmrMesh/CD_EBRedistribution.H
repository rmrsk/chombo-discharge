/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBRedistribution.H
  @brief  Declaration of a class which can perform redistribution in an EB-AMR context.
  @author Robert Marskar
*/

#ifndef CD_EBRedistribution_H
#define CD_EBRedistribution_H

// Chombo includes
#include <EBLevelGrid.H>
#include <LevelData.H>
#include <EBCellFAB.H>
#include <BaseIVFAB.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for redistribution ala Chombo's flux redistribution. 
  @details This uses the same algorithm as the core Chombo code, but includes much faster define functions and smaller memory footprint. Because
  this class is used in a non-subcycling context, it has way fewer bells and whistles.
*/
class EBRedistribution
{
public:
  /*!
    @brief Weak constructor. Need to call define afterwards
  */
  EBRedistribution() noexcept;

  /*!
    @brief Copy constructor is not allowed.
  */
  EBRedistribution(const EBRedistribution& a_other) = delete;

  /*!
    @brief Full constructor. Calls the define function and puts object in usable state. 
    @details If there is a finer level, a_eblgFine and a_eblgCoarsenedFine must be defined. If there is a coarser grid level then
    a_eblgCoar and a_eblgRefinedCoar must be defined. 
    @param[in] a_eblgCoar             Coarser grid level
    @param[in] a_eblgRefinedCoar      Refined coarse grid
    @param[in] a_eblg                 Grids on this level
    @param[in] a_eblgCoarsenedFined   Coarsened fine grid
    @param[in] a_eblgCoarsenedFined   Coarsened fine grid
    @param[in] a_refToCoar            Refinement factor between this level and the coarse level
    @param[in] a_refToFine            Refinement factor between this level and the fine level
    @param[in] a_redistributionRadius Redistribution radius. Must be >= 1
  */
  EBRedistribution(const EBLevelGrid& a_eblgCoar,
                   const EBLevelGrid& a_eblgRefinedCoar,
                   const EBLevelGrid& a_eblg,
                   const EBLevelGrid& a_eblgCoarsenedFine,
                   const EBLevelGrid& a_eblgFine,
                   const int          a_refToCoar,
                   const int          a_refToFine,
                   const int          a_redistributionRadius) noexcept;

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~EBRedistribution() noexcept;

  /*!
    @brief Define fucntion. Puts object in usable state. 
    @details If there is a finer level, a_eblgFine and a_eblgCoarsenedFine must be defined. If there is a coarser grid level then
    a_eblgCoar and a_eblgRefinedCoar must be defined. 
    @param[in] a_eblgCoar             Coarser grid level
    @param[in] a_eblgRefinedCoar      Refined coarse grid
    @param[in] a_eblg                 Grids on this level
    @param[in] a_eblgCoarsenedFined   Coarsened fine grid
    @param[in] a_eblgCoarsenedFined   Coarsened fine grid
    @param[in] a_refToCoar            Refinement factor between this level and the coarse level
    @param[in] a_refToFine            Refinement factor between this level and the fine level
    @param[in] a_redistributionRadius Redistribution radius. Must be >= 1
  */
  virtual void
  define(const EBLevelGrid& a_eblgCoar,
         const EBLevelGrid& a_eblgRefinedCoar,
         const EBLevelGrid& a_eblg,
         const EBLevelGrid& a_eblgCoarsenedFine,
         const EBLevelGrid& a_eblgFine,
         const int          a_refToCoar,
         const int          a_refToFine,
         const int          a_redistributionRadius) noexcept;

protected:
  /*!
    @brief Is defined or not
  */
  bool m_isDefined;

  /*!
    @brief Has coarse level or not
  */
  bool m_hasCoar;

  /*!
    @brief Has fine level or not
  */
  bool m_hasFine;

  /*!
    @brief Refinement factor between this level and the coarse level
  */
  int m_refToCoar;

  /*!
    @brief Refinement factor between this level and the fine level
  */
  int m_refToFine;

  /*!
    @brief Grids on coarser level
  */
  EBLevelGrid m_eblgCoar;

  /*!
    @brief Grids on this level that are a refinement of the coarser level
  */
  EBLevelGrid m_eblgRefinedCoar;

  /*!
    @brief Grids on this level
  */
  EBLevelGrid m_eblg;

  /*!
    @brief Grids on this level that are a coarsening of the fine level
  */
  EBLevelGrid m_eblgCoarsenedFine;

  /*!
    @brief Fine grids
  */
  EBLevelGrid m_eblgFine;

  /*!
    @brief Stencils for redistribution into valid grid cells on the coarse level. 
    @note Defined on this level (m_eblg), but the stencils contain VoFs on the coarse level
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_redistStencilsCoar;

  /*!
    @brief Stencils for redistribution into valid grid cells on this level. 
    @note Defined on this level (m_eblg), and the stnecils contain VoFs on this level also.
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_redistStencils;

  /*!
    @brief Stencils for redistribution into valid grid cells on the fine level. 
    @note Defined on this level (m_eblg), but the stencils contain VoFs on the fine level
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_redistStencilsFine;
};

#include <CD_NamespaceFooter.H>
