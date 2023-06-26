#ifndef _CD_EBHELMHOLTZOPF_F_H_
#define _CD_EBHELMHOLTZOPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C" {

#ifndef GUARDAPPLYOPREGULAR
#define GUARDAPPLYOPREGULAR
// Prototype for Fortran procedure applyopregular ...
//
void FORTRAN_NAME(APPLYOPREGULAR, applyopregular)(CHFp_FRA1(lphi),
                                                  CHFp_CONST_FRA1(phi),
                                                  CHFp_CONST_FRA1(aco),
                                                  CHFp_CONST_FRA1(b0),
                                                  CHFp_CONST_FRA1(b1),
                                                  CHFp_CONST_FRA1(b2),
                                                  CHFp_CONST_REAL(alpha),
                                                  CHFp_CONST_REAL(beta),
                                                  CHFp_CONST_REAL(dx),
                                                  CHFp_BOX(region));

#define FORT_APPLYOPREGULAR FORTRAN_NAME(inlineAPPLYOPREGULAR, inlineAPPLYOPREGULAR)
#define FORTNT_APPLYOPREGULAR FORTRAN_NAME(APPLYOPREGULAR, applyopregular)

inline void FORTRAN_NAME(inlineAPPLYOPREGULAR, inlineAPPLYOPREGULAR)(CHFp_FRA1(lphi),
                                                                     CHFp_CONST_FRA1(phi),
                                                                     CHFp_CONST_FRA1(aco),
                                                                     CHFp_CONST_FRA1(b0),
                                                                     CHFp_CONST_FRA1(b1),
                                                                     CHFp_CONST_FRA1(b2),
                                                                     CHFp_CONST_REAL(alpha),
                                                                     CHFp_CONST_REAL(beta),
                                                                     CHFp_CONST_REAL(dx),
                                                                     CHFp_BOX(region))
{
  CH_TIMELEAF("FORT_APPLYOPREGULAR");
  FORTRAN_NAME(APPLYOPREGULAR, applyopregular)
  (CHFt_FRA1(lphi),
   CHFt_CONST_FRA1(phi),
   CHFt_CONST_FRA1(aco),
   CHFt_CONST_FRA1(b0),
   CHFt_CONST_FRA1(b1),
   CHFt_CONST_FRA1(b2),
   CHFt_CONST_REAL(alpha),
   CHFt_CONST_REAL(beta),
   CHFt_CONST_REAL(dx),
   CHFt_BOX(region));
}
#endif // GUARDAPPLYOPREGULAR
}

#endif
