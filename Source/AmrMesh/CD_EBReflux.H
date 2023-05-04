/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBReflux.H
  @brief  Declaration of a class which can reflux over the coarse-fine interface. 
  @author Robert Marskar
*/

#ifndef CD_EBReflux_H
#define CD_EBReflux_H

// Std includes
#include <map>

// Chombo includes
#include <Copier.H>
#include <EBLevelGrid.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class which can do refluxing across a coarse-fine interface. 
  @details This is like Chombos EBFluxRegister, but trimmed for reduced memory overhead and performance. 
*/
class EBReflux
{
public:
  /*!
    @brief Refluxing type. 
  */
  enum class Type
  {
    Conservative,
    Arithmetic
  };
  /*!
    @brief Disallowed constructor. 
  */
  EBReflux() noexcept;

  /*!
    @brief Disallowed copy constructor.
  */
  EBReflux(const EBReflux& a_other) = delete;

  /*!
    @brief Full constructor. Calls the define function. 
    @param[in] a_eblg         Grids
    @param[in] a_eblgFine     Fine grids
    @param[in] a_eblgCoFi     Coarsened fine grids
    @param[in] a_refRat       Refinement ratio between grids. 
  */
  EBReflux(const EBLevelGrid& a_eblg,
           const EBLevelGrid& a_eblgFine,
           const EBLevelGrid& a_eblgCoFi,
           const int          a_refRat) noexcept;

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~EBReflux() noexcept;

  /*!
    @brief Define method. Puts object in usable state. 
    @param[in] a_eblg     Grids on this level
    @param[in] a_eblgFine Fine grids
    @param[in] a_eblgCoFi Coarsened fine grids
    @param[in] a_refRat    Refinement ratio between coarse and fine grids.
  */
  virtual void
  define(const EBLevelGrid& a_eblg,
         const EBLevelGrid& a_eblgFine,
         const EBLevelGrid& a_eblgCoFi,
         const int          a_refRat) noexcept;

  /*!
    @brief Reflux into the coarse data. 
    @details This subtracts the coarse-fluxes from the EB and adds in the sum offine-grid fluxes. Multiplication by area 
    fractions and grid resolutions are done in this routine. 
    @param[inout] a_data      Data
    @param[in]    a_flux      Fluxes
    @param[in]    a_fineFlux  Fine-grid flux. 
    @param[in]    a_variables Variables to reflux.
  */
  virtual void
  reflux(LevelData<EBCellFAB>&       a_Lphi,
         const LevelData<EBFluxFAB>& a_flux,
         const LevelData<EBFluxFAB>& a_fineFlux,
         const Interval              a_variables,
         const Real                  a_scaleCoarFlux,
         const Real                  a_scaleFineFlux) const noexcept;

protected:
  /*!
    @brief Is defined or not
  */
  bool m_isDefined;

  /*!
    @brief Fine grids
  */
  EBLevelGrid m_eblgFine;

  /*!
    @brief Grid on this level. 
  */
  EBLevelGrid m_eblg;

  /*!
    @brief Coarsened fine grids. 
  */
  EBLevelGrid m_eblgCoFi;

  /*!
    @brief Regular coarse-fine regions for coarse-grid patches
  */
  LayoutData<std::map<std::pair<int, Side::LoHiSide>, DenseIntVectSet>> m_regularCoarseFineRegions;

  /*!
    @brief Cut-cells in each patch that lie on the coarse-fine interface
    @note Defined over the coarse grids
  */
  mutable LayoutData<std::map<std::pair<int, Side::LoHiSide>, VoFIterator>> m_irregularCoarseFineRegions;

  /*!
    @brief Refinement factor between fine and coarse level. 
  */
  int m_refRat;

  /*!
    @brief Define CF regions
  */
  virtual void
  defineRegionsCF() noexcept;

  /*!
    @brief 
  */
  virtual void
  coarsenFluxesCF(LevelData<EBFluxFAB>&       a_coarFluxes,
                  const LevelData<EBFluxFAB>& a_fineFluxes,
                  const int                   a_coarVar,
                  const int                   a_fineVar) const noexcept;

  virtual void
  refluxIntoCoarse(LevelData<EBCellFAB>&       a_Lphi,
                   const LevelData<EBFluxFAB>& a_oldFluxes,
                   const LevelData<EBFluxFAB>& a_newFluxes,
                   const int                   a_phiVar,
                   const int                   a_oldFluxVar,
                   const int                   a_newFluxVar,
                   const Real                  a_scaleCoarFlux,
                   const Real                  a_scaleFineFlux) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
