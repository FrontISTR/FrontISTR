/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2009/02/16                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by GOTO, Kazuya (AdvanceSoft)                    *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/



#ifndef HECMW_DIST_REFINE_INCLUDED
#define HECMW_DIST_REFINE_INCLUDED

#include "hecmw_struct.h"

extern int HECMW_dist_refine(struct hecmwST_local_mesh **mesh, int refine, const char *cad_filename, const char *part_filename);

#endif
