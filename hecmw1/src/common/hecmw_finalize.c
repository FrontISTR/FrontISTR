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



#include "hecmw_finalize.h"
#include "hecmw_util.h"

int 
HECMW_finalize(void)
{
	HECMW_log(HECMW_LOG_DEBUG, "Finalizing..."); 

	HECMW_ctrl_finalize();

#ifndef HECMW_SERIAL
	MPI_Finalize();

	HECMW_log(HECMW_LOG_DEBUG, "MPI finalized");
#endif

	return 0;
}
