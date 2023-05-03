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
    @brief Disallowed constructor. 
  */
  EBReflux() noexcept;

  /*!
    @brief Disallowed copy constructor.
  */
  EBReflux(const EBReflux& a_other) = delete;

  /*!
    @brief Full constructor. Calls the define function. 
    @param[in] a_eblgFine     Fine grids
    @param[in] a_eblgCoFi     Coarsened fine grids
    @param[in] a_eblgCoar     Coarse grids
    @param[in] a_refRat       Refinement ratio between grids. 
  */
  EBReflux(const EBLevelGrid& a_eblgFine,
           const EBLevelGrid& a_eblgCoFi,
           const EBLevelGrid& a_eblgCoar,
           const int          a_refRat) noexcept;

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~EBReflux() noexcept;

  /*!
    @brief Define method. Puts object in usable state. 
    @param[in] a_eblgFine     Fine grids
    @param[in] a_eblgCoFi     Coarsened fine grids
    @param[in] a_eblgCoar     Coarse grids
    @param[in] a_refRat       Refinement ratio between coarse and fine grids.
  */
  virtual void
  define(const EBLevelGrid& a_eblgFine,
         const EBLevelGrid& a_eblgCoFi,
         const EBLevelGrid& a_eblgCoar,
         const int          a_refRat) noexcept;

  /*!
    @brief Reflux into the coarse data. 
    @param[in]    a_phiCoar   Coarse-grid data
    @param[in]    a_variables Variables to reflux
  */
  virtual void
  reflux(LevelData<EBCellFAB>& a_phiCoar, const Interval a_variables) const noexcept;

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
    @brief Coarse grids
  */
  EBLevelGrid m_eblgCoar;

  /*!
    @brief Coarsened fine grids. 
  */
  EBLevelGrid m_eblgCoFi;

  /*!
    @brief Refinement factor between fine and coarse level. 
  */
  int m_refRat;

  /*!
    @brief Number of ghost cells to fill across coarse-fine interface. 
    @note Only applies the the cut-cells. 
  */
  int m_ghostCF;
};

#include <CD_NamespaceFooter.H>

#endif
