/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuaki Sakane (RIST)                         *
 *                       Shin'ichi Ezure (RIST)                        *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/



#ifndef INC_HECMW_ETYPE
#define INC_HECMW_ETYPE

#include "hecmw_util.h"

extern int
HECMW_get_etype_UTIL2HECMW( int etype );

extern int
HECMW_get_etype_HECMW2UTIL( int etype );

extern int
HECMW_get_etype_GeoFEM2HECMW( int etype );

extern int
HECMW_get_max_node( int etype );

extern int
HECMW_get_max_edge( int etype );

extern int
HECMW_get_max_surf( int etype );

extern int
HECMW_get_max_tsuf( int etype );

extern int
HECMW_get_max_qsuf( int etype );

extern char
*HECMW_get_ucd_label( int etype );


extern int
HECMW_is_etype_rod(int etype);


extern int
HECMW_is_etype_surface(int etype);


extern int
HECMW_is_etype_solid(int etype);


extern int
HECMW_is_etype_interface(int etype);


extern int
HECMW_is_etype_beam(int etype);


extern int
HECMW_is_etype_shell(int etype);


extern int
HECMW_is_etype_link(int etype);


extern int
HECMW_is_etype_33struct(int etype);

#endif

