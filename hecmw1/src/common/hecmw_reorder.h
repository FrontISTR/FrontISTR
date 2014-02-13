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



#ifndef INC_HECMW_REORDER
#define INC_HECMW_REORDER

#include "hecmw_struct.h"

extern int
HECMW_reorder_node_mpc( struct hecmwST_local_mesh *local_mesh );
extern int
HECMW_reorder_node_dof( struct hecmwST_local_mesh *local_mesh );
extern int
HECMW_reorder_elem_type( struct hecmwST_local_mesh *local_mesh );
extern int
HECMW_reorder( struct hecmwST_local_mesh *local_mesh );

#endif
