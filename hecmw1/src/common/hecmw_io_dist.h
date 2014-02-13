/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuaki Sakane (RIST)                         *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/

#ifndef INC_HECMW_IO_DIST
#define INC_HECMW_IO_DIST

#include "hecmw_struct.h"

extern struct hecmwST_local_mesh *HECMW_get_dist_mesh( char *fname );

extern int HECMW_put_dist_mesh( const struct hecmwST_local_mesh *mesh, char *fname );

#endif
