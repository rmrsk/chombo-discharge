/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBGradient.cpp
  @brief  Implementation of CD_EBGradient.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>
#include <EBArith.H>
#include <CH_Timer.H>

// Our includes
#include <CD_Timer.H>
#include <CD_EBGradient.H>
#include <CD_gradientF_F.H>
#include <CD_LoadBalancing.H>
#include <CD_NamespaceHeader.H>

constexpr int EBGradient::m_comp;
constexpr int EBGradient::m_nComp;

EBGradient::EBGradient(const EBLevelGrid& a_eblg,
		       const EBLevelGrid& a_eblgFine,
		       const CellLocation a_dataLocation,
		       const Real         a_dx,
		       const int          a_refRat,
		       const int          a_order,
		       const int          a_weighting){
  CH_TIME("EBGradient::EBGradient");

  CH_assert(a_order    >  0);
  CH_assert(a_weight   >= 0);  
  CH_assert(a_refRat%2 == 0);

  m_eblg         = a_eblg;
  m_eblgFine     = a_eblgFine;
  m_dataLocation = a_dataLocation;
  m_dx           = a_dx;
  m_refRat       = a_refRat;
  m_order        = a_order;
  m_weighting    = a_weighting;

  if(a_eblgFine.isDefined()){
    m_hasFine  = true;
    m_dxFine   = m_dx/m_refRat;
  }
  else{
    m_hasFine  = false;
    m_dxFine   = 1;
  }

  Timer timer("EBGradient::EBGradient");

  timer.startEvent("Define level stencils");
  this->defineLevelStencils();
  timer.stopEvent("Define level stencils");

  if(m_hasFine){
    timer.startEvent("Define masks");
    this->defineMasks();
    timer.stopEvent("Define masks");        
    timer.startEvent("Define EBCF stencils");    
    this->defineStencilsEBCF();
    timer.stopEvent("Define EBCF stencils");    
  }

  ParmParse pp("EBGradient");
  bool profile = false;
  pp.query("profile", profile);

  if(profile){
    timer.eventReport(pout(), false);
  }
}

EBGradient::~EBGradient(){
  CH_TIME("EBGradient::~EBGradient");  
}

void EBGradient::computeLevelGradient(LevelData<EBCellFAB>&       a_gradient,
				      const LevelData<EBCellFAB>& a_phi) const {
  CH_TIME("EBGradient::computeLevelGradient");

  CH_assert(a_gradient.nComp() == SpaceDim);
  CH_assert(a_phi.     nComp() == 1       );

  const DisjointBoxLayout& dbl   = m_eblg.getDBL();
  const EBISLayout&        ebisl = m_eblg.getEBISL();

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    EBCellFAB&       grad    = a_gradient[dit()];
    const EBCellFAB& phi     = a_phi     [dit()];
    const EBISBox&   ebisBox = ebisl     [dit()];        
    const Box        cellBox = dbl       [dit()];

    const bool isAllCovered = ebisBox.isAllCovered();

    if(!isAllCovered){
      // Compute regular cells
      BaseFab<Real>&       gradFAB = grad.getSingleValuedFAB();
      const BaseFab<Real>& phiFAB  = phi. getSingleValuedFAB();
      FORT_GRADIENT(CHF_FRA(gradFAB),
		    CHF_CONST_FRA1(phiFAB, m_comp),
		    CHF_CONST_REAL(m_dx),
		    CHF_BOX(cellBox));

      // Now do the boundary stencils.
      const BaseIVFAB<VoFStencil>& bndryStencils = m_bndryStencils[dit()];

      VoFIterator& vofit = m_bndryIterator[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit){
	const VolIndex&   vof  = vofit();
	const VoFStencil& sten = bndryStencils(vof, m_comp);

	for (int dir = 0; dir < SpaceDim; dir++){
	  grad(vof, dir) = 0.0;
	}

	for (int i = 0; i < sten.size(); i++){
	  const VolIndex& ivof    = sten.vof(i);
	  const Real&     iweight = sten.weight(i);
	  const int&      ivar    = sten.variable(i); // Note: For the gradient stencil the component is the direction. 

	  grad(vof, ivar) += phi(ivof, m_comp)*iweight;
	}
      }
    }
    else{
      for (int dir = 0; dir < SpaceDim; dir++){
	grad.setCoveredCellVal(0.0, dir);
      }
    }
  }
}


void EBGradient::computeGradient(LevelData<EBCellFAB>&       a_gradient,
				 const LevelData<EBCellFAB>& a_phi) const {
  CH_TIME("EBGradient::computeGradient(no finer)");

  this->computeLevelGradient(a_gradient, a_phi);

  // Do corrections near the EBCF. 
  if(m_hasFine){

    const DisjointBoxLayout& dbl = m_eblg.getDBL();

    for (DataIterator dit(dbl); dit.ok(); ++dit){
      EBCellFAB&       grad = a_gradient[dit()];
      const EBCellFAB& phi  = a_phi[dit()];

      // Do corrections near the EBCF. 
      VoFIterator& vofit = m_ebcfIterator[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit){
	const VolIndex&   vof  = vofit();
	const VoFStencil& sten = m_ebcfStencils[dit()](vof, m_comp);

	for (int dir = 0; dir < SpaceDim; dir++){
	  grad(vof, dir) = 0.0;
	}

	for (int i = 0; i < sten.size(); i++){
	  const VolIndex& ivof    = sten.vof(i);
	  const Real&     iweight = sten.weight(i);
	  const int&      ivar    = sten.variable(i); // Note: For the gradient stencil the component is the direction. 

	  grad(vof, ivar) += phi(ivof, m_comp)*iweight;
	}	
      }
    }
  }
}

void EBGradient::defineLevelStencils(){
  CH_TIME("EBGradient::defineLevelStencils");

  const DisjointBoxLayout& dbl    = m_eblg.getDBL();
  const EBISLayout&        ebisl  = m_eblg.getEBISL();
  const ProblemDomain&     domain = m_eblg.getDomain();

  m_bndryStencils.define(dbl);
  m_bndryIterator.define(dbl);  

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box      cellBox  = dbl       [dit()];
    const EBISBox& ebisbox  = ebisl     [dit()];
    const EBGraph& ebgraph  = ebisbox.getEBGraph();

    // Determine which cells that need explicit stencils. This is certainly all the cut-cells
    // and the ones at the domain boundaries.
    IntVectSet bndryIVS = ebisbox.getIrregIVS(cellBox);
    for (int dir = 0; dir < SpaceDim; dir++){
      Box loBox;
      Box hiBox;
      int hasLo;
      int hasHi;
	
      EBArith::loHi(loBox, hasLo, hiBox, hasHi, domain, cellBox, dir);
	
      if(hasLo){
	bndryIVS |= IntVectSet(loBox);
      }
      if(hasHi){
	bndryIVS |= IntVectSet(hiBox);
      }
    }

    VoFIterator&           bndryIterator = m_bndryIterator[dit()];        
    BaseIVFAB<VoFStencil>& bndryStencils = m_bndryStencils[dit()];

    bndryIterator.define(bndryIVS, ebgraph);
    bndryStencils.define(bndryIVS, ebgraph, m_nComp);

    for (bndryIterator.reset(); bndryIterator.ok(); ++bndryIterator){
      const VolIndex& vof = bndryIterator();
      VoFStencil& stencil = bndryStencils(vof, m_comp);
      
      stencil.clear();

      for (int dir = 0; dir < SpaceDim; dir++){
	VoFStencil derivDirStencil;
	EBArith::getFirstDerivStencilWidthOne(derivDirStencil, vof, ebisbox, dir, m_dx, nullptr, dir);
	stencil += derivDirStencil;	
      }
    }
  }
}

void EBGradient::defineStencilsEBCF(){
  CH_TIME("EBGradient::defineLevelStencils");

  const DisjointBoxLayout& dbl        = m_eblg.    getDBL();
  const DisjointBoxLayout& dblFine    = m_eblgFine.getDBL();
  
  const EBISLayout&        ebisl      = m_eblg.    getEBISL();
  const EBISLayout&        ebislFine  = m_eblgFine.getEBISL();
  
  const ProblemDomain&     domain     = m_eblg.    getDomain();
  const ProblemDomain&     domainFine = m_eblgFine.getDomain();

  // Define iterators and stencils for EBCF
  m_ebcfStencils.define(dbl);
  m_ebcfIterator.define(dbl);  

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box      cellBox = dbl  [dit()];
    const EBISBox& ebisBox = ebisl[dit()];
    const EBGraph& ebgraph = ebisBox.getEBGraph();

    const bool isAllRegular = ebisBox.isAllRegular();
    const bool isAllCovered = ebisBox.isAllCovered();
    const bool isIrregular  = !isAllRegular && !isAllCovered;

    const DenseIntVectSet& coarseFineRegion = m_coarseFineRegion[dit()];
    const DenseIntVectSet& invalidRegion    = m_invalidRegion   [dit()];    

    // Determine cells where we need to drop order. 
    IntVectSet ebcfIVS;
    if(isIrregular){

      // Iterate through the coarse-fine region and check if the finite difference stencil reaches into a cut-cell.
      for (DenseIntVectSetIterator divsIt(coarseFineRegion); divsIt.ok(); ++divsIt){
	const IntVect iv = divsIt();

	// Get all the vofs in this cell and check if we can use centered differencing in them. We can do this
	bool needSpecialStencil = false;
	
	const Vector<VolIndex> vofs = ebisBox.getVoFs(iv);
	for (int i = 0; i < vofs.size(); i++){
	  const VolIndex& vof = vofs[i];

	  // Get a regular finite-difference stencil for cell-centered data. 
	  VoFStencil gradSten;
	  for (int dir = 0; dir < SpaceDim; dir++){
	    
	    VoFStencil derivSten;
	    EBArith::getFirstDerivStencilWidthOne(derivSten, vof, ebisBox, dir, m_dx, nullptr, dir);
	    gradSten += derivSten;
	  }

	  // Check if we reach into stuff we shouldn't.
	  for (int j = 0; j < gradSten.size(); j++){
	    const VolIndex& ivof = gradSten.vof(j);
	    const IntVect   gid  = ivof.gridIndex();
	      
	    const bool isCoveredByFinerCell = invalidRegion[gid];
	    const bool isIrregularCell      = ebisBox.isIrregular(gid);

	    if(isCoveredByFinerCell && isIrregularCell){
	      needSpecialStencil = true;
	    }
	  }
	}

	if(needSpecialStencil){
	  ebcfIVS |= iv;
	}
      }
    }

    // Define the iterator and stencils for places where we need to drop order. 
    VoFIterator&           vofitEBCF    = m_ebcfIterator[dit()];
    BaseIVFAB<VoFStencil>& stencilsEBCF = m_ebcfStencils[dit()];

    vofitEBCF.   define(ebcfIVS, ebgraph);
    stencilsEBCF.define(ebcfIVS, ebgraph, m_nComp);

    for (vofitEBCF.reset(); vofitEBCF.ok(); ++vofitEBCF){
      const VolIndex& vof = vofitEBCF();
      const IntVect   iv  = vof.gridIndex();
      
      VoFStencil& stencil = stencilsEBCF(vof, m_comp);
      stencil.clear();

      for (int dir = 0; dir < SpaceDim; dir++){
	VoFStencil derivStencil;
	
	// Check for cells to the low/high side of this one.
	const IntVect ivLo = iv - BASISV(dir);
	const IntVect ivHi = iv + BASISV(dir);

	const bool isLoCellCovered    = invalidRegion[ivLo];
	const bool isHiCellCovered    = invalidRegion[ivHi];

	const bool isLoCellIrregular  = ebisBox.isIrregular(ivLo);
	const bool isHiCellIrregular  = ebisBox.isIrregular(ivHi);

	const Vector<VolIndex> loVoFs = ebisBox.getVoFs(vof, dir, Side::Lo, 1);
	const Vector<VolIndex> hiVoFs = ebisBox.getVoFs(vof, dir, Side::Hi, 1);

	const int numLoVoFs = loVoFs.size();
	const int numHiVoFs = hiVoFs.size();

	const bool useLoCell = (numLoVoFs > 0) && !(isLoCellIrregular && isLoCellCovered);
	const bool useHiCell = (numHiVoFs > 0) && !(isHiCellIrregular && isHiCellCovered);			

	// Use centered differencing in dir-direction if we can. Otherwise, drop order to forward/backward differences. 
	if(useLoCell && useHiCell){
	  for (int i = 0; i < numHiVoFs; i++){
	    derivStencil.add(hiVoFs[i], 1./numHiVoFs, dir);
	  }

	  for (int i = 0; i < numLoVoFs; i++){
	    derivStencil.add(loVoFs[i], -1./numLoVoFs, dir);
	  }

	  derivStencil *= 1./(2.0*m_dx);
	}
	else if(useLoCell){ 
	  derivStencil.add(vof, 1.0, dir);

	  for (int i = 0; i < numLoVoFs; i++){
	    derivStencil.add(loVoFs[i], -1./numLoVoFs, dir);
	  }

	  derivStencil *= 1./m_dx;
	}
	else if(useHiCell){ 
	  derivStencil.add(vof, -1.0, dir);

	  for (int i = 0; i < numHiVoFs; i++){
	    derivStencil.add(hiVoFs[i], 1./numHiVoFs, dir);
	  }

	  derivStencil *= 1./m_dx;
	}
	
	stencil += derivStencil;	
      }
    }
  }
}

void EBGradient::defineMasks(){
  CH_TIME("EBGradient::defineMasks");

  CH_assert(m_hasFine);

  constexpr Real     zero   = 0.0;
  constexpr Real     one    = 1.0;    
  const     Interval interv = Interval(m_comp, m_comp);

  // Handle to computational grids. 
  const DisjointBoxLayout& dbl        = m_eblg.    getDBL();
  const DisjointBoxLayout& dblFine    = m_eblgFine.getDBL();
  
  const EBISLayout&        ebisl      = m_eblg.    getEBISL();
  const EBISLayout&        ebislFine  = m_eblgFine.getEBISL();
  
  const ProblemDomain&     domain     = m_eblg.    getDomain();
  const ProblemDomain&     domainFine = m_eblgFine.getDomain();

  // Create some temporary storage. We coarsen the fine grid and create a BoxLayoutData<FArrayBox> object where each
  // Box is grown by one. We do this because we want to put some data in the ghost cells outside the valid region
  // of the fine level. 
  DisjointBoxLayout dblCoFi;
  coarsen(dblCoFi, dblFine, m_refRat);

  Vector<int> coFiRanks      = dblCoFi.procIDs();  
  Vector<Box> coFiBoxes      = dblCoFi.boxArray();
  Vector<Box> grownCoFiBoxes = dblCoFi.boxArray();
  for(int i = 0; i < coFiBoxes.size(); i++){
    grownCoFiBoxes[i].grow(1);
  }

  BoxLayout coFiLayout     (     coFiBoxes, coFiRanks);
  BoxLayout grownCoFiLayout(grownCoFiBoxes, coFiRanks);
  
  BoxLayoutData<FArrayBox> coFiMaskCF     (grownCoFiLayout, m_nComp); 
  BoxLayoutData<FArrayBox> coFiMaskInvalid(     coFiLayout, m_nComp); 

  // Build the mask regions. 
  for (DataIterator dit(coFiLayout); dit.ok(); ++dit){
    const Box cellBox  =      coFiLayout[dit()];
    const Box grownBox = grownCoFiLayout[dit()];

    // Set the "coarse-fine" mask. It is set to 1 in cells that are not valid cells. We do this by setting the FArrayBox data to one
    // in the entire patch, and then we iterate through the parts of the patch that overlap with valid grid boxes. That part of the
    // data is set to zero. We end up with a bunch of grid patches which hold a value of 1 in ghost cells that do not overlap with
    // the valid region of the grid. 
    coFiMaskCF[dit()].setVal(one);
    coFiMaskCF[dit()].setVal(zero, cellBox, m_comp, m_nComp);
    for (LayoutIterator lit = dblCoFi.layoutIterator(); lit.ok(); ++lit){
      const Box neighborBox = dblCoFi[lit()];
      const Box overlapBox  = neighborBox & grownBox;

      if(!overlapBox.isEmpty()){
	coFiMaskCF[dit()].setVal(zero, overlapBox, m_comp, m_nComp);
      }
    }

    // Set the invalid cell mask on the coarsened fine grids. We set it to 1 everywhere on the fine grid and 0 elsewhere. When we add this to the coarse
    // level we will add 1 into every cell that is covered by a finer level. 
    coFiMaskInvalid[dit()].setVal(one, cellBox, m_comp, m_nComp);
  }

  // Add to the coarMaskCF. After this, coarMaskCF will have a value of 1 in all cells that abut the fine level. Likewise, coarMaskInvalid
  // will have a value of 1 in all cells that are covered by a finer level.
  LevelData<FArrayBox> coarMaskCF     (dbl, m_nComp, IntVect::Zero); // CF region always lies on the valid region in the coarse grid. 
  LevelData<FArrayBox> coarMaskInvalid(dbl, m_nComp, IntVect::Unit); // Stencil's have width 1, so need one ghost layer to check for valid cells
  
  for (DataIterator dit(dbl); dit.ok(); ++dit){
    coarMaskCF     [dit()].setVal(zero);
    coarMaskInvalid[dit()].setVal(zero);
  }

  coFiMaskCF.     addTo(interv, coarMaskCF,      interv, ebisl.getDomain()); // Contains > 0 in coarse cells that abut the refinement bounary. 
  coFiMaskInvalid.addTo(interv, coarMaskInvalid, interv, ebisl.getDomain()); // Contains > 0 in coarse cells that are covered by a finer level.

  coarMaskInvalid.exchange();

  // Define the boolean (DenseIntVectSet) masks and iterate through the FArrayBoxData data in order to determine
  // if the mask should be set to true. 
  m_coarseFineRegion.define(dbl);
  m_invalidRegion.   define(dbl);
  
  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box cellBox  = dbl[dit()];
    const Box grownBox = grow(cellBox, 1);

    // Set to false by default. 
    m_coarseFineRegion[dit()] = DenseIntVectSet(cellBox,  false);
    m_invalidRegion   [dit()] = DenseIntVectSet(grownBox, false);
    
    DenseIntVectSet& coarseFineRegion = m_coarseFineRegion[dit()];
    DenseIntVectSet& invalidRegion    = m_invalidRegion   [dit()];

    const FArrayBox& fabMaskCF      = coarMaskCF     [dit()];
    const FArrayBox& fabMaskInvalid = coarMaskInvalid[dit()];

    for (BoxIterator bit(cellBox); bit.ok(); ++bit){
      const IntVect iv = bit();

      // Set the coarse-fine region. 
      if(fabMaskCF(iv, m_comp) > zero){
	coarseFineRegion |= iv;
      }
    }

    // Set the invalid mask. 
    for(BoxIterator bit(grownBox); bit.ok(); ++bit){
      const IntVect iv = bit();
      
      if(fabMaskInvalid(iv, m_comp) > zero){
	invalidRegion |= iv;
      }
    }

#if 1 // Bogus check
    for (BoxIterator bit(cellBox); bit.ok(); ++bit){
      const IntVect iv = bit();
      
      const DenseIntVectSet& coarseFineRegion = m_coarseFineRegion[dit()];
      const DenseIntVectSet& invalidRegion    = m_invalidRegion   [dit()];

      
      if(coarseFineRegion[iv] && invalidRegion[iv]){
	std::cout << "Logic bust in cell = " << iv << " on domain = " << domain << std::endl;
      }
    }
#endif
  }
}


#include <CD_NamespaceFooter.H>
