/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbIntersection.H
  @brief  Enum for putting some logic into EB intersection tests. 
  @author Robert Marskar
*/

#ifndef CD_EbIntersection_H
#define CD_EbIntersection_H

#include <CD_NamespaceHeader.H>

/*!
  @brief Enum for putting some logic into how we think about intersection between particles and EBs. 
*/
enum class EbIntersection
{
  Bisection,
  Raycast
};

#include <CD_NamespaceFooter.H>

#endif
