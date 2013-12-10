/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : Coupling Interface                                *
 *                                                                     *
 *            Written by Shin'ichi Ezure (RIST)                        *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using Hight End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/




#ifndef INC_HECMW_COUPLE_COPY_C2F
#define INC_HECMW_COUPLE_COPY_C2F

#include "hecmw_couple_startup.h"

extern int
HECMW_couple_copy_c2f_init(struct hecmw_couple_value *_couple_value);

extern int
HECMW_couple_copy_c2f_finalize(void);

#endif	/* INC_HECMW_COUPLE_COPY_C2f */
