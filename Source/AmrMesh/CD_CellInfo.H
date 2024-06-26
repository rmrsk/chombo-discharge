/* chombo-discharge
 * Copyright © 2024 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CellInfo.H
  @brief  Simple class for holding some quantities relevant in a grid cell. 
  @author Robert Marskar
*/

#ifndef CD_CellInfo_H
#define CD_CellInfo_H

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for the cell-information that is often queried when merging particles inside a cell. 
*/
class CellInfo
{
public:
  /*!
    @brief Default constructor - creates a regular cell
  */
  CellInfo() = delete;

  /*!
    @brief Constructor which creates a regular cell 
    @param[in] a_gridIndex     Grid index
    @param[in] a_dx            Grid resolution
  */
  CellInfo(const IntVect a_gridIndex, const Real a_dx) noexcept;

  /*!
    @brief Constructor which automatically constructs an approximation to the minimum box enclosing the valid region of the cell
    @param[in] a_gridIndex     Grid index
    @param[in] a_dx            Grid resolution
    @param[in] a_volFrac       Volume fraction. 
    @param[in] a_bndryCentroid Boundary centroid
    @param[in] a_bndryNormal   Boundary normal
  */
  CellInfo(const IntVect  a_gridIndex,
           const Real     a_dx,
           const Real     a_volFrac,
           const RealVect a_bndryCentroid,
           const RealVect a_bndryNormal) noexcept;

  /*!
    @brief Destructor
  */
  virtual ~CellInfo() noexcept;

  /*!
    @brief Get the grid index
  */
  IntVect&
  getGridIndex() noexcept;

  /*!
    @brief Get the grid index
  */
  const IntVect&
  getGridIndex() const noexcept;

  /*!
    @brief Get the grid spacing
  */
  Real&
  getDx() noexcept;

  /*!
    @brief Get the grid spacing
  */
  const Real&
  getDx() const noexcept;

  /*!
    @brief Get the volume fraction
  */
  Real&
  getVolFrac() noexcept;

  /*!
    @brief Get the volume fraction
  */
  const Real&
  getVolFrac() const noexcept;

  /*!
    @brief Get the boundary centroid
  */
  RealVect&
  getBndryCentroid() noexcept;

  /*!
    @brief Get the boundary centroid
  */
  const RealVect&
  getBndryCentroid() const noexcept;

  /*!
    @brief Get the boundary normal
  */
  RealVect&
  getBndryNormal() noexcept;

  /*!
    @brief Get the boundary normal
  */
  const RealVect&
  getBndryNormal() const noexcept;

  /*!
    @brief Get the lower valid region
  */
  RealVect&
  getValidLo() noexcept;

  /*!
    @brief Get the lower valid region
  */
  const RealVect&
  getValidLo() const noexcept;

  /*!
    @brief Get the upper valid region
  */
  RealVect&
  getValidHi() noexcept;

  /*!
    @brief Get the upper valid region
  */
  const RealVect&
  getValidHi() const noexcept;

protected:
  /*!
    @brief Grid index
  */
  IntVect m_gridIndex;

  /*!
    @brief Grid spacing
  */
  Real m_dx;

  /*!
    @brief Volume fraction
  */
  Real m_volFrac;

  /*!
    @brief Boundary centroid
  */
  RealVect m_bndryCentroid;

  /*!
    @brief EB boundary normal
  */
  RealVect m_bndryNormal;

  /*!
    @brief Valid region in each direction (lower limit)
  */
  RealVect m_validLo;

  /*!
    @brief Valid region in each direction (upper limit)
  */
  RealVect m_validHi;
};

#include <CD_NamespaceFooter.H>

#endif
