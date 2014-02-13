/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : HEC-MW Utility                                    *
 *                                                                     *
 *            Written by Shin'ichi Ezure (RIST)                        *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/


#ifndef INC_HECMW_MESH_EDGE_INFO
#define INC_HECMW_MESH_EDGE_INFO

#include "hecmw_struct.h"
#include "hecmw_part_struct.h"

extern int
HECMW_mesh_edge_info( struct hecmwST_local_mesh *mesh, struct hecmw_part_edge_data *edge_data );

#endif  /* INC_HECMW_MESH_EDGE_INFO */
