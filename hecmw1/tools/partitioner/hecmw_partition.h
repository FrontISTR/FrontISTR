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


#ifndef INC_HECMW_PARTITION
#define INC_HECMW_PARTITION

#include "hecmw_struct.h"
#include "hecmw_part_struct.h"

extern struct hecmwST_local_mesh
*HECMW_partition_inner( struct hecmwST_local_mesh *global_mesh,
                        struct hecmw_part_cont_data *cont_data );

extern struct hecmwST_local_mesh
*HECMW_partition( struct hecmwST_local_mesh *local_mesh );

#endif  /* INC_HPWMC_PARTITION */
