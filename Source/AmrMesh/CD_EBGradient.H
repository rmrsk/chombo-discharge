/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBGradient.H
  @brief  Declaration of a class which computes gradients in an EBAMR context. 
  @author Robert Marskar
*/

#ifndef CD_EBGradient_H
#define CD_EBGradient_H

// Std includes
#include <map>

// Chombo includes
#include <EBLevelGrid.H>

// Our includes
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for computing gradients in an EBAMR context.
  @details Although computing the gradient might seem like an easy thing to do, the presence of a coarse-fine boundary together with an EB complicates
  both the stencil availability and the consistent handling of the gradient on the coarse level. This class introduces a two-level view of such stencil
  computations. This class uses regular second order finite differencing in regular regions, including those covered by a finer level. Near the EB, however,
  the stencil might reach into cut-cells that are covered by the finer level. The invalid data might be filled with bogus data, or at the very best coarsened
  from the fine level. Since there is no consistent way of doing that coarsening (both conserative and direct both make sense), this class switches to a least
  squares based approach for such cases. In those particular corner cases we only reach into valid regions in the domain. Because of that, the directional
  derivatives in those cells are in the form

       d(phi)/dx = sum(vofCoar) w(vofCoar)*phi(vofCoar) + sum(vofFine)[w(vofFine)*phi(vofFine)]

  i.e., the gradient is a two-level stencil.

  Near the domain boundaries this routine switches to first order stencils. 
  @note If there is no finer level, this class will compute regular second-order stencils away from the EB and switch to compact width-one stencils if
  there are not enough cells available near the EB. 
  @note This class requires at least one ghost cell. 
  @note Coarse data underneath the fine grid is considered bogus data!
*/
class EBGradient {
public:

  /*!
    @brief Specification of cell-center or cell-centroid type of interpolation
  */
  using CellLocation = Location::Cell;

  /*!
    @brief Disallowed constructor. 
  */
  EBGradient() = delete;

  /*!
    @brief Disallowed copy constructor.
  */
  EBGradient(const EBGradient& a_other) = delete;

  /*!
    @brief Full constructor, assumes that there is a finer level. 
    @param[in] a_eblg         Grids
    @param[in] a_eblgFine     Fine grids
    @param[in] a_dataLocation Interpration of data centering. Either on the cell center or the cell centroid. 
    @param[in] a_dx           Resolution
    @param[in] a_refRat       Refinement ratio
    @param[in] a_order        Interpolation order for modified gradients. 
    @param[in] a_weighting    Weighting factor least squares. 
  */
  EBGradient(const EBLevelGrid& a_eblg,
	     const EBLevelGrid& a_eblgFine,
	     const CellLocation a_dataLocation,
	     const Real         a_dx,
	     const int          a_refRat,
	     const int          a_order,
	     const int          a_weighting = 0);

  /*!
    @brief Compute gradient. 
    @param[out] a_gradient Gradient of input scalar
    @param[in]  a_phi      Input scalar
  */
  void computeLevelGradient(LevelData<EBCellFAB>&       a_gradient,
			    const LevelData<EBCellFAB>& a_phi) const;  

  /*!
    @brief Compute gradient. 
    @param[out] a_gradient Gradient of input scalar
    @param[in]  a_phi      Input scalar
  */
  void computeGradient(LevelData<EBCellFAB>&       a_gradient,
		       const LevelData<EBCellFAB>& a_phi) const;

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~EBGradient();
  
protected:

  /*!
    @brief Scalar component to reach into
  */  
  static constexpr int m_comp = 0;

  /*!
    @brief Number of components that are supported
  */
  static constexpr int m_nComp = 1;

  /*!
    @brief How to interpret data locations. 
  */
  CellLocation m_dataLocation;

  /*!
    @brief Grids
  */
  EBLevelGrid m_eblg;

  /*!
    @brief Fine grids
  */
  EBLevelGrid m_eblgFine;

  /*!
    @brief True if we have a finer level (m_eblgFine is defined)
  */
  bool m_hasFine;

  /*!
    @brief Resolution
  */
  Real m_dx;

  /*!
    @brief Fine resolution
  */
  Real m_dxFine;  

  /*!
    @brief Order
  */
  int m_order;

  /*!
    @brief Weight (for least squares)
  */
  int m_weighting;

  /*!
    @brief Refinement factor between fine and coarse level. 
  */
  int m_refRat;

  /*!
    @brief Stencils for boundaries and cut-cells
  */
  LayoutData<BaseIVFAB<VoFStencil> > m_bndryStencils;

  /*!
    @brief Iterator for boundary cells (either domain or EB)
  */
  mutable LayoutData<VoFIterator> m_bndryIterator;

  /*!
    @brief EBCF stencils
  */
  LayoutData<BaseIVFAB<VoFStencil> > m_ebcfStencils;

  /*!
    @brief EBCF iterators. 
  */
  mutable LayoutData<VoFIterator> m_ebcfIterator;

  /*!
    @brief Map of cells that lie on the coarse-fine interface
    @details This is built in defineMasks(). If true, the cell abuts the fine grid. 
  */
  LayoutData<DenseIntVectSet> m_coarseFineRegion;

  /*!
    @brief Map of cells that are covered by the finer level.
    @details This is built in defineMasks(). If true, the cell is covered by a finer level. 
  */
  LayoutData<DenseIntVectSet> m_invalidRegion;

  /*!
    @brief Define level stencils. 
  */
  void defineLevelStencils();

  /*!
    @brief Define stencils for interpolation across the EBCF refinement boundary. 
  */
  void defineStencilsEBCF();

  /*!
    @brief Make masks
  */
  void defineMasks();
};

#include <CD_NamespaceFooter.H>

#endif
