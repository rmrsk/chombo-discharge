/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBGhostCellInterpolator.H
  @brief  Declaration of a class which interpolates ghost cells across the coarse-fine interface. 
  @author Robert Marskar
*/

#ifndef CD_EBGhostCellInterpolator_H
#define CD_EBGhostCellInterpolator_H

// Std includes
#include <map>

// Chombo includes
#include <EBLevelGrid.H>
#include <QuadCFInterp.H>
#include <AggStencil.H>

// Our includes
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class which can interpolate ghost cells across the coarse-fine interface.

  To select an interpolation method, use the Type argument. 
*/
class EBGhostCellInterpolator
{
public:
  /*!
    @brief Type of interpolation method supported
  */
  /*!
    @brief Type of interpolation methods supported.
    PWC = Piecewise constant, ignoring the embedded boundary. This is not conservative. 
    ConservativePWC = Piecewise constant, including the EB volume/area fractions.
    ConservativeMinMod = Min-mod slope for interpolation, including the EB volume/area fractions.
    ConservativeMonotonizedCentral = Monotonized central slope for interpolation, including the EB volume/area. 
  */
  enum Type
  {
    PWC,
    MinMod,
    MonotonizedCentral,
    Superbee,
    ConservativePWC,
    ConservativeMinMod,
    ConservativeMonotonizedCentral,
    ConservativeSuperbee,
  };

  /*!
    @brief Disallowed constructor. 
  */
  EBGhostCellInterpolator() noexcept;

  /*!
    @brief Disallowed copy constructor.
  */
  EBGhostCellInterpolator(const EBGhostCellInterpolator& a_other) = delete;

  /*!
    @brief Full constructor. Calls the define function. 
    @param[in] a_eblgFine     Fine grids
    @param[in] a_eblgCoFi     Coarsened fine grids
    @param[in] a_eblgCoar     Coarse grids
    @param[in] a_dataLocation Interpration of data centering. Either on the cell center or the cell centroid. 
    @param[in] a_ghostVector  (Minimum) number of ghost cells in input/output data. 
    @param[in] a_refRat       Refinement ratio between coarse and fine grids.
    @param[in] a_ghostCF      Width of the ghost region to be filled. Only relevant near the EBCF. 
  */
  EBGhostCellInterpolator(const EBLevelGrid& a_eblgFine,
                          const EBLevelGrid& a_eblgCoFi,
                          const EBLevelGrid& a_eblgCoar,
                          const IntVect&     a_ghostVector,
                          const int          a_refRat,
                          const int          a_ghostCF) noexcept;

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~EBGhostCellInterpolator() noexcept;

  /*!
    @brief Define method. Puts object in usable state. 
    @param[in] a_eblgFine     Fine grids
    @param[in] a_eblgCoFi     Coarsened fine grids
    @param[in] a_eblgCoar     Coarse grids
    @param[in] a_dataLocation Interpration of data centering. Either on the cell center or the cell centroid. 
    @param[in] a_ghostVector  (Minimum) number of ghost cells in input/output data. 
    @param[in] a_refRat       Refinement ratio between coarse and fine grids.
    @param[in] a_ghostCF      Width of the ghost region to be filled. Only relevant near the EBCF. 
  */
  virtual void
  define(const EBLevelGrid& a_eblgFine,
         const EBLevelGrid& a_eblgCoFi,
         const EBLevelGrid& a_eblgCoar,
         const IntVect&     a_ghostVector,
         const int          a_refRat,
         const int          a_ghostCF) noexcept;

  /*!
    @brief Do inhomogeneous interpolation
    @param[inout] a_phiFine   Fine phi
    @param[in]    a_phiCoar   Coarse phi
    @param[in]    a_variables Variables to interpolate
    @param[in]    a_type      Interpolation type
  */
  virtual void
  interpolate(LevelData<EBCellFAB>&       a_phiFine,
              const LevelData<EBCellFAB>& a_phiCoar,
              const Interval              a_variables,
              const Type                  a_interpType) const noexcept;

protected:
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
    @brief Is defined or not
  */
  bool m_isDefined;

  /*!
    @brief Minimum number of ghost cells in input data. 
  */
  IntVect m_ghostVector;

  /*!
    @brief Refinement factor between fine and coarse level. 
  */
  int m_refRat;

  /*!
    @brief Number of ghost cells to fill across coarse-fine interface. 
    @note Only applies the the cut-cells. 
  */
  int m_ghostCF;

  /*!
    @brief Define buffer storage
  */
  virtual void
  defineBuffer() const noexcept;

  /*!
    @brief Delete buffer storage
  */
  virtual void
  undefineBuffer() const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
