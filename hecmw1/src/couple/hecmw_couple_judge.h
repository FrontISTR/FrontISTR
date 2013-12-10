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




#ifndef INC_HECMW_COUPLE_JUDGE
#define INC_HECMW_COUPLE_JUDGE

#include "hecmw_struct.h"

extern int
HECMW_couple_judge_tet1(const struct hecmwST_local_mesh *local_mesh, int elem, int surf_id,
		double coord_px, double coord_py, double coord_pz, double *dot_product, double *distance);
extern int
HECMW_couple_judge_hex1(const struct hecmwST_local_mesh *local_mesh, int elem, int surf_id,
		double coord_px, double coord_py, double coord_pz, double *dot_product, double *distance);

#endif	/* INC_HECMW_COUPLE_JUDGE */
