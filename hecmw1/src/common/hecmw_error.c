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



#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "hecmw_util.h"
#include "hecmw_error.h"


static int hecmw_errno;


static char hecmw_errmsg[HECMW_MSG_LEN+1];


int
HECMW_set_verror(int errorno, const char *fmt, va_list ap)
{
	char errmsg[HECMW_MSG_LEN+1];

	hecmw_errno = errorno;

	HECMW_snprintf(hecmw_errmsg, sizeof(hecmw_errmsg), "%s", HECMW_strmsg(errorno));
	HECMW_vsnprintf(errmsg, sizeof(errmsg), fmt, ap);

	if(strlen(errmsg) > 0) {
		HECMW_snprintf(hecmw_errmsg+strlen(hecmw_errmsg), sizeof(hecmw_errmsg)-strlen(hecmw_errmsg), " (%s)", errmsg);
	}

	HECMW_print_error();

	return 0;
}


int
HECMW_set_error(int errorno, const char *fmt, ...)
{
	int rc;
	va_list ap;

	rc = 0;
	va_start(ap, fmt);
	rc = HECMW_set_verror(errorno, fmt, ap);
	va_end(ap);

	return rc;
}


int
HECMW_get_error(char **errmsg)
{
	if(errmsg) {
		*errmsg = hecmw_errmsg;
	}
	return hecmw_errno;
}


int
HECMW_get_errno(void)
{
	return hecmw_errno;
}


char *
HECMW_get_errmsg(void)
{
	return hecmw_errmsg;
}

