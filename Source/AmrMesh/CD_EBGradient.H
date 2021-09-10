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
    @brief Full constructor. 
    @param[in] a_eblg         Grids
    @param[in] a_eblgFine     Fine grids
    @param[in] a_dataLocation Interpration of data centering. Either on the cell center or the cell centroid. 
    @param[in] a_refRat       Refinement ratio
    @param[in] a_order        Interpolation order for modified gradients. 
    @param[in] a_weighting    Weighting factor least squares. 
  */
  EBGradient(const EBLevelGrid& a_eblg,
	     const EBLevelGrid& a_eblgFine,
	     const CellLocation a_dataLocation,
	     const int          a_refRat,
	     const int          a_order,
	     const int          a_weighting = 0);
			  

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~EBGradient();

  
protected:

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
    @brief Refined coarse grids
  */
  EBLevelGrid m_eblgCoFi;

  /*!
    @brief Refinement factor between fine and coarse level. 
  */
  int m_refRat;

  /*!
    @brief Interpolation order. 
  */
  int m_order;

  /*!
    @brief Least squares weighting factor
  */
  int m_weighting;

  /*!
    @brief Ghost cells that will be interpolated near the cut-cells. 
    @note This does not include the ghost cells on the refinement boundary away from the cut-cells. 
  */
  LayoutData<IntVectSet> m_ghostCells;

  /*!
    @brief Fine-grid interpolation stencils.
  */
  LayoutData<BaseIVFAB<VoFStencil> > m_fineStencils;

  /*!
    @brief Coarse-grid interpolation stencils.
  */
  LayoutData<BaseIVFAB<VoFStencil> > m_coarStencils;

  /*!
    @brief Regular coarse-fine interface boxes. 
    @note This is a map to the "regular" ghost cells that will be interpolated using a standard
    multigrid interpolation stencil (away from the EB). 
  */
  LayoutData<std::map<std::pair<int, Side::LoHiSide>, Box> > m_cfivs;

  /*!
    @brief Defines the CFIVS. 
  */
  void defineCFIVS();

  /*!
    @brief Define ghost cells to be filled in range m_ghostCF from each cut-cell. 
  */
  void defineGhostRegions();

  /*!
    @brief Define grids
    @param[in] a_eblgFine Fine grid
    @param[in] a_eblgCoar Coarse grid
  */
  void defineGrids(const EBLevelGrid& a_eblgFine, const EBLevelGrid& a_eblgCoar);

  /*!
    @brief Define data holders. 
  */
  void defineBuffers();

  /*!
    @brief Define stencils for interpolation across the EBCF refinement boundary. 
  */
  void defineStencilsEBCF();

  /*!
    @brief Build the AggStencil objects for faster stencil aggregation
  */
  void makeAggStencils();

  /*!
    @brief Build an interpolation stencil
    @param[out] a_stencilFine    Fine-grid stencil
    @param[out] a_stencilCoar    Coarse-grid stencil
    @param[in]  a_dataLocation   Data location
    @param[in]  a_ghostVofFine   Ghost cell which is to be filled on the coarse grid. 
    @param[in]  a_ghostVofCoar   Vof corresponding to a_ghostVofFine on the coarse grid
    @param[in]  a_ebisBoxFine    Fine grid EBISBox
    @param[in]  a_ebisBxoCoar    Coarse grid EBISBox
    @param[in]  a_validFineCells Fine-grid region which we are allowed to use in the interpolation. 
    @param[in]  a_validCoarCells Coarse-grid region which we are allowed to use in the interpolation. 
    @param[in]  a_dxFine         Grid resolution on the fine level
    @param[in]  a_dxCoar         Grid resolution on the coarse level
    @param[in]  a_order          Interpolation order
    @param[in]  a_weight         Equation weighting for least squares reconstruction. 
    @return Returns true if a stencil was found, in which case it is defined in a_stencilFine and a_stencilCoar. 
  */
  bool getStencil(VoFStencil&            a_stencilFine,
		  VoFStencil&            a_stencilCoar,
		  const CellLocation&    a_dataLocation,
		  const VolIndex&        a_ghostVofFine,
		  const VolIndex&        a_ghostVofCoar,
		  const EBISBox&         a_ebisboxFine,
		  const EBISBox&         a_ebisboxCoar,
		  const DenseIntVectSet& a_validFineCells,
		  const DenseIntVectSet& a_validCoarCells,
		  const Real&            a_dxFine,
		  const Real&            a_dxCoar,
		  const int&             a_order,
		  const int&             a_weight);
};

#include <CD_NamespaceFooter.H>

#endif
