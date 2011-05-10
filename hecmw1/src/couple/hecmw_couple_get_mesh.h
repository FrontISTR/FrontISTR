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




#ifndef INC_HECMW_COUPLE_GET_MESH
#define INC_HECMW_COUPLE_GET_MESH

#include "hecmw_struct.h"

extern struct hecmwST_local_mesh *
HECMW_couple_get_mesh(char *name_ID, char *unit_ID);

#endif	/* INC_HECMW_COUPLE_GET_MESH */
