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
#include <CH_Timer.H>

// Our includes
#include <CD_Timer.H>
#include <CD_EBGradient.H>
#include <CD_gradientF_F.H>

#include <CD_NamespaceHeader.H>

EBGradient::EBGradient(const EBLevelGrid& a_eblg,
		       const EBLevelGrid& a_eblgFine,
		       const CellLocation a_dataLocation,
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
  m_refRat       = a_refRat;
  m_order        = a_order;
  m_weighting    = a_weighting;

  Timer timer("EBGradient::EBGradient");

  
}

EBGradient::~EBGradient(){
  CH_TIME("EBGradient::~EBGradient");  
}

#include <CD_NamespaceFooter.H>
