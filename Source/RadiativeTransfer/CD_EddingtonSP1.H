/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EddingtonSP1.H
  @brief  Declaration of a first order Eddington solver for radiative transfer
  @author Robert Marskar
*/

#ifndef CD_EddingtonSP1_H
#define CD_EddingtonSP1_H

// Std includes
#include <random>

// Chombo includes
#include <AMRMultiGrid.H>
#include <AMRTGA.H>
#include <BiCGStabSolver.H>
#include <EBBackwardEuler.H>
#include <EBSimpleSolver.H>

// Our includes
#include <CD_RtSolver.H>
#include <CD_EbHelmholtzOpFactory.H>
#include <CD_EddingtonSP1DomainBc.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Radiative tranfer equation solver in the SP1 (diffusion) approximation. Also includes stochastically drawn Photons. 
*/
class EddingtonSP1 : public RtSolver {
public:

  /*!
    @brief Default function for space-time dependence of domain boundary conditions. 
    @details This is a utility function 
    @param[in] a_position Position in space
    @param[in] a_time     Time
    @return This returns one (1.0)
  */
  static Real s_defaultDomainBcFunction(const RealVect a_position, const Real a_time);

  /*!
    @brief Constructor
  */
  EddingtonSP1();

  /*!
    @brief Destructor
  */
  virtual ~EddingtonSP1();

  /*!
    @brief Advance RTE onto state a_phi
  */
  virtual bool advance(const Real a_dt, EBAMRCellData& a_phi, const EBAMRCellData& a_source, const bool a_zerophi = false);

  /*!
    @brief Parse class options
  */
  virtual void parseOptions();

  /*!
    @brief Parse class options
  */
  virtual void parseRuntimeOptions();

  /*!
    @brief Set the boundary condition function on a wall. 
    @param[in] a_dir Coordinate direction.
    @param[in] a_side Side in the coordinate direction.
    @param[in] a_function Boundary condition function.
    @details This sets a boundary condition for a particular side on a wall. The user must also specify how to use this BC in the input script. 
  */
  virtual void setDomainBcWallFunction(const int a_dir, const Side::LoHiSide a_side, const EddingtonSP1DomainBc::BcFunction& a_function);
  
  /*!
    @brief Allocate internal storage
  */
  virtual void allocateInternals();
  
  /*!
    @brief Cache state
  */
  virtual void preRegrid(const int a_base, const int a_oldFinestLevel);

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocateInternals();

  /*!
    @brief Regrid function for this class
  */
  virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Register operators
  */
  virtual void registerOperators();

  /*!
    @brief Compute the boundary flux. For Eddington, the boundary flux is = c*phi/2
  */
  virtual void computeBoundaryFlux(EBAMRIVData& a_ebFlux, const EBAMRCellData& a_phi);

  /*!
    @brief Compute the domain flux. For Eddington, the domain flux is = c*phi/2
  */
  virtual void computeDomainFlux(EBAMRIFData& a_domainflux, const EBAMRCellData& a_phi);

  /*!
    @brief Compute the flux
    @details For diffusive models, the flux will be something like grad(a_phi). For higher-order models, the flux 
    will be contained in a_phi (somehow);
  */
  virtual void computeFlux(EBAMRCellData& a_flux, const EBAMRCellData& a_phi);

  /*!
    @brief Get isotropic part. 
  */
  virtual void computeDensity(EBAMRCellData& a_isotropic, const EBAMRCellData& a_phi);

  /*!
    @brief Write plot file
  */
  virtual void writePlotFile();

  /*!
    @brief Write checkpoint data into handle
  */
  virtual void writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const;

  /*!
    @brief Read checkpoint data from handle
  */
  virtual void readCheckpointLevel(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Return required number of ghost cells
  */
  virtual int queryGhost() const;

protected:

  /*!
    @brief Enum class for supported bottom solvers in multigrid.
  */
  enum class BottomSolver {
    Simple,
    BiCGStab,
    GMRES
  };

  /*!
    @brief Enum for supported relaxation methods in multigrid. 
  */
  enum class RelaxationMethod {
    Jacobi,
    GaussSeidel,
    GSRBFast
  };

  /*!
    @brief Enum for multigrid cycle types. 
  */
  enum class MultigridType{
    FAS,
    VCycle,
    FCycle,
  };

  /*!
    @brief Enum for boundary condition types on EBs
  */
  enum class EbBcType{
    Dirichlet,
    Neumann,
    Robin
  };
  
  /*!
    @brief Relaxation type for gmg
  */
  RelaxationMethod m_multigridRelaxMethod;

  /*!
    @brief GMG multigrid type
  */
  MultigridType m_multigridType;

  /*!
    @brief Needs setup
  */
  bool m_needsMultigridSetup;

  /*!
    @brief Has deeper MG levels or not
  */
  bool m_hasDeeperMultigridLevels;

  /*!
    @brief Use tga for transient solevs
  */
  bool m_useTGA;

  /*!
    @brief Use the source term to draw stochastic Photons
  */
  bool m_stochastic_Photons;

  /*!
    @brief Verbosity for geometric multigrid
  */
  int m_multigridVerbosity;

  /*!
    @brief Precoarsening for multigrid
  */
  int m_numCoarseningsBeforeAggregation;

  /*!
    @brief Number of smoothings before averaging
  */
  int m_multigridPreSmooth;

  /*!
    @brief Number of smoothings before averaging
  */
  int m_multigridPostSmooth;

  /*!
    @brief Number of smoothing before bottom solver
  */
  int m_multigridBottomSmooth;

  /*!
    @brief Maximum number of iterations
  */
  int m_multigridMaxIterations;

  /*!
    @brief Minimum number of iterations
  */
  int m_multigridMinIterations;

  /*!
    @brief Bottom solver type
  */
  BottomSolver m_bottomSolver;

  /*!
    @brief Number of smoothing for bottom solver
  */
  int m_numSmoothingsForSimpleSolver;

  /*!
    @brief Set bottom drop depth
  */
  int m_numCellsBottomDrop;

  /*!
    @brief RNG seed
  */
  int m_seed;

  /*!
    @brief 
  */
  Real m_multigridTolerance;
  
  /*!
    @brief 
  */
  Real m_multigridHang;

  /*!
    @brief Reflection coefficient
  */
  Real m_reflectionCoefficientOne;

  /*!
    @brief Reflection coefficient
  */
  Real m_reflectionCoefficientTwo;

  /*!
    @brief RNG
  */
  std::mt19937_64* m_rng;

  /*!
    @brief Geometric multigrid solver
  */
  RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > m_multigridSolver;

  /*!
    @brief Conjugate gradient solver bottom MG level
  */
  BiCGStabSolver<LevelData<EBCellFAB> > m_bicgstab;

  /*!
    @brief TGA solver
  */
  RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > > m_tgaSolver;

  /*!
    @brief Euler solver
  */
  RefCountedPtr<EBBackwardEuler> m_eulerSolver;

  /*!
    @brief Operator factory
  */
  RefCountedPtr<EbHelmholtzOpFactory> m_operatorFactory;

  /*!
    @brief MG levelsgrids
  */
  Vector<EBLevelGrid> m_mg_levelgrids;

  /*!
    @brief multi-fluid simple solver
  */
  EBSimpleSolver m_simpleSolver;
  
  /*!
    @brief Residue for stationary solver
  */
  EBAMRCellData m_resid;

  /*!
    @brief a-coefficient
  */
  EBAMRCellData m_aco;

  /*!
    @brief b-coefficient
  */
  EBAMRFluxData m_bco;

  /*!
    @brief b-coefficient
  */
  EBAMRIVData m_bco_irreg;

  /*!
    @brief Wrapper calss. 
  */
  EddingtonSP1DomainBc m_domainBc;

  /*!
    @brief Actual functions on domain edges (faces). 
  */
  std::map<EddingtonSP1DomainBc::Wall, EddingtonSP1DomainBc::BcFunction> m_domainBcFunctions;

  /*!
    @brief Boundary condition on the embedded boundary
  */
  std::pair<EbBcType, Real> m_ebbc;

  /*!
    @brief Set up geometric multigrid
  */
  virtual void setupMultigrid();

  /*!
    @brief Set kappa
  */
  virtual void setMultigridCoefficients();

  /*!
    @brief Set the a-coefficient and b-coefficient
  */
  virtual void setACoefAndBCoef();

  /*!
    @brief Set both aco and bco (irregular version)
  */
  virtual void setACoefAndBCoefBox(EBCellFAB&       a_aco,
				   BaseIVFAB<Real>& a_bco_irreg,
				   const Box        a_cellbox,
				   const RealVect   a_origin,
				   const Real       a_dx,
				   const int        a_lvl,
				   const DataIndex& a_dit);

  /*!
    @brief Define the GMG levels
  */
  virtual void defineDeeperMultigridLevels();
  
  /*!
    @brief Set the operator factory
  */
  virtual void setupOperatorFactory();

  /*!
    @brief Set the GMG solver
  */
  virtual void setupMultigridSolver();

  /*!
    @brief Make domain bc string
  */
  virtual std::string makeBcString(const int a_dir, const Side::LoHiSide a_side) const;

  /*!
    @brief Returns BC type based on string.
    @param[in] a_str Boundary condition string. Must be "dirichlet number", "neumann number", "robin_number", "dirichlet_custom", "neumann_custom", or "robin_custom".
  */
  virtual EddingtonSP1DomainBc::BcType parseBcString(const std::string a_str) const;

  /*!
    @brief Set up the TGA solver
  */
  virtual void setupTGA();
  
  /*!
    @brief Set up the Euler solver
  */
  virtual void setupEuler();

  /*!
    @brief Parse domain BC settings
  */
  virtual void parseDomainBc();

    /*!
    @brief Parse domain BC settings
  */
  virtual void parseEbBc();

  /*!
    @brief Parse whether or not this is a stationary solver
  */
  virtual void parseStationary();   

  /*!
    @brief Parse plot variables
  */
  virtual void parsePlotVariables();

  /*!
    @brief Parse multigrid settings
  */
  virtual void parseMultigridSettings();

  /*!
    @brief Parse reflection coefficients for Robin bcs
  */
  virtual void parseReflection();

  /*!
    @brief Set default domain BC functions.
  */
  virtual void setDefaultDomainBcFunctions();
};

#include <CD_NamespaceFooter.H>

#endif