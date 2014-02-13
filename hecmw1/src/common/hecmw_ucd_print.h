/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Shin'ichi Ezure (RIST)                        *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/



#ifndef INC_HECMW_UCD_PRINT
#define INC_HECMW_UCD_PRINT

#include "hecmw_struct.h"
#include "hecmw_result.h"

extern int
HECMW_ucd_print( const struct hecmwST_local_mesh *mesh,
                 const struct hecmwST_result_data *result,
                 const char *ofname );
extern int
HECMW_ucd_legacy_print( const struct hecmwST_local_mesh *mesh,
                        const struct hecmwST_result_data *result,
                        const char *ofname );

#endif

