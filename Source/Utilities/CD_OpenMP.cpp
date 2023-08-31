/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_OpenMP.cpp
  @brief  Implementation of CD_OpenMP.H
  @author Robert Marskar
*/

// Our includes
#include <CD_OpenMP.H>
#include <CD_NamespaceHeader.H>

// Thread-safe reduction operator for taking the union of IntVectSet
void
ThreadSafeIVSUnion(IntVectSet& ivsInOut, const IntVectSet& ivsIn)
{
  ivsInOut |= ivsIn;
}

#include <CD_NamespaceFooter.H>
