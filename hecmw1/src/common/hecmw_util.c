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
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <ctype.h>
#include <errno.h>
#include "hecmw_config.h"
#include "hecmw_util.h"

void
HECMW_fprintf(FILE *fp, char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	vfprintf(fp, fmt, ap);
	va_end(ap);
}


void
HECMW_printerr(char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);
}


char *
HECMW_get_date(void)
{
	int rc;
	time_t now;
	static char static_buf[100];

	if(time(&now) == (time_t)-1) return NULL;
	rc = strftime(static_buf, sizeof(static_buf),
					"%b %d %H:%M:%S", localtime(&now));
	return rc ? static_buf : NULL;
}


void
HECMW_assert_(int cond, char *cond_str, char *file, int line)
{
	if(!cond) {
		fprintf(stderr, "%s:%d: Assersion `%s' failed.\n", file, line, cond_str);
#ifdef HECMW_SERIAL
		abort();
#else
		MPI_Abort(MPI_COMM_WORLD, HECMW_EXIT_ERROR);
#endif
	}
}


int
HECMW_check_condition_( int cond, char *cond_str, int isabort, char *file, int line )
{

  if( cond ) return 0;

	if( isabort ) {
    fprintf( stderr, "%s:%d: Assertion `%s' falied.\n", file, line, cond_str );
#ifdef HECMW_SERIAL
		abort( );
#else
		MPI_Abort( MPI_COMM_WORLD, HECMW_EXIT_ERROR );
#endif
	}
	return 1;
}


void
HECMW_abort(HECMW_Comm comm)
{
#ifdef HECMW_SERIAL
	abort();
#else
/*	HECMW_comm_is_initialized() ? MPI_Abort(comm, HECMW_EXIT_ERROR) : abort(); */
	if( HECMW_comm_is_initialized() ) {
		MPI_Abort(comm, HECMW_EXIT_ERROR);
	} else {
		abort();
	}
#endif
}


char *
HECMW_toupper(char *s)
{
	char *p;

	if(s == NULL) return NULL;

	for(p=s; *p; p++) {
		*p = (char)toupper(*p);
	}
	return s;
}


char *
HECMW_tolower(char *s)
{
	char *p;

	if(s == NULL) return NULL;

	for(p=s; *p; p++) {
		*p = (char)tolower(*p);
	}
	return s;
}


void
HECMW_print_error(void)
{
	HECMW_log(HECMW_LOG_ERROR, HECMW_get_errmsg());
}


void
HECMW_print_vmsg(int loglv, int msgno, const char *fmt, va_list ap)
{
	char msg[HECMW_MSG_LEN+1];
	char vmsg[HECMW_MSG_LEN+1];

	HECMW_snprintf(msg, sizeof(msg), "%s", HECMW_strmsg(msgno));
	HECMW_vsnprintf(vmsg, sizeof(vmsg), fmt, ap);
	if(strlen(vmsg) > 0) {
		HECMW_snprintf(msg+strlen(msg), sizeof(msg)-strlen(msg), " (%s)", vmsg);
	}
	HECMW_log(loglv, msg);
}


void
HECMW_print_msg(int loglv, int msgno, const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	HECMW_print_vmsg(loglv, msgno, fmt, ap);
	va_end(ap);
}


int
HECMW_vsnprintf(char *str, size_t size, const char *format, va_list ap)
{
#ifdef WIN32_MSVC
	return _vsnprintf(str, size, format, ap);
#else
	return vsnprintf(str, size, format, ap);
#endif
}


int
HECMW_snprintf(char *str, size_t size, const char *format, ...)
{
	va_list ap;
	int rtc;

	va_start(ap, format);
#ifdef WIN32_MSVC
	rtc = _vsnprintf(str, size, format, ap);
#else
	rtc = vsnprintf(str, size, format, ap);
#endif
	va_end(ap);

	return rtc;
}
