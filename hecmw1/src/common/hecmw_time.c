/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Noboru Imai (Univ. of Tokyo)                  *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/


#include "hecmw_time.h"


double HECMW_Wtime(void)
{
#ifndef HECMW_SERIAL
	double t;
	t = MPI_Wtime();
	return t;
#else
	struct timeb t;
	double sec;
	ftime(&t);
	sec = t.time + (double)t.millitm * 1e-3;
	return sec;
#endif
}

double HECMW_Wtick(void)
{
#ifndef HECMW_SERIAL
	return MPI_Wtick();
#else
	return 1e-3;
#endif
}	

/* interface for fortran */


double hecmw_wtime_fi(void) { return HECMW_Wtime(); }
double hecmw_wtime_fi_(void) { return HECMW_Wtime(); }
double hecmw_wtime_fi__(void) { return HECMW_Wtime(); }
double HECMW_WTIME_FI(void) { return HECMW_Wtime(); }
double HECMW_WTIME_FI_(void) { return HECMW_Wtime(); }
double HECMW_WTIME_FI__(void) { return HECMW_Wtime(); }


double hecmw_wtick_fi(void) { return HECMW_Wtick(); }
double hecmw_wtick_fi_(void) { return HECMW_Wtick(); }
double hecmw_wtick_fi__(void) { return HECMW_Wtick(); }
double HECMW_WTICK_FI(void) { return HECMW_Wtick(); }
double HECMW_WTICK_FI_(void) { return HECMW_Wtick(); }
double HECMW_WTICK_FI__(void) { return HECMW_Wtick(); }




