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



#include <stdio.h>
#include "hecmw_util.h"
#include "hecmw_init.h"
/* #include "hecmw_couple_info.h"  2007/12/27 S.Ito   */


int
HECMW_init_ex(int *argc, char ***argv, const char *ctrlfile)
{
	if(HECMW_comm_init(argc, argv)) return -1;
	HECMW_log(HECMW_LOG_DEBUG, "Initilalizing..."); 
	if(ctrlfile == NULL) ctrlfile = HECMW_CTRL_FILE;
	if(HECMW_ctrl_init_ex(ctrlfile)) return -1;
/*     if(HECMW_couple_comm_init() != HECMW_SUCCESS) return -1;  2007/12/27 S.Ito */
	return 0;
}



int
HECMW_init(int *argc, char ***argv)
{
	return HECMW_init_ex(argc, argv, HECMW_CTRL_FILE);
}

