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



#ifndef HECMW_ERROR_INCLUDED
#define HECMW_ERROR_INCLUDED

#include <stdarg.h>


extern int HECMW_set_verror(int errorno, const char *fmt, va_list ap);


extern int HECMW_set_error(int errorno, const char *fmt, ...);


extern int HECMW_get_error(char **errmsg);


extern int HECMW_get_errno(void);


extern char *HECMW_get_errmsg(void);

#endif
