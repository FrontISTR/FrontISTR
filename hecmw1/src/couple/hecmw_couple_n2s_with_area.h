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




#ifndef INC_HECMW_COUPLE_N2S_WITH_AVERAGE
#define INC_HECMW_COUPLE_N2S_WITH_AVERAGE

#include "hecmw_struct.h"
#include "hecmw_couple_weight.h"

extern struct hecmw_couple_weight_list *
HECMW_couple_n2s_with_area(const struct hecmwST_local_mesh *mesh,
		const struct hecmw_couple_boundary *boundary);

#endif	/* INC_HECMW_COUPLE_N2S_WITH_AVERAGE */
