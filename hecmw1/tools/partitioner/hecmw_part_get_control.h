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


#ifndef INC_HECMW_PART_GET_CONTROL
#define INC_HECMW_PART_GET_CONTROL

#include "hecmw_part_struct.h"

extern int
HECMW_part_set_ctrl_file_name( char *fname );

extern struct hecmw_part_cont_data
*HECMW_part_get_control( );

extern void
HECMW_part_free_control( struct hecmw_part_cont_data *cont_data );

#endif  /* INC_HECMW_PART_GET_CONTROL */
