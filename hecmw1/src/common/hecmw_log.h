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

 

#ifndef HECMW_LOG_INCLUDED
#define HECMW_LOG_INCLUDED

#include <stdarg.h>


#define HECMW_LOG_MAX 10

/* ログレベル定義 */

#define HECMW_LOG_NONE	0

#define HECMW_LOG_ERROR	1

#define HECMW_LOG_WARN	2

#define HECMW_LOG_INFO	4

#define HECMW_LOG_DEBUG	8

#define HECMW_LOG_ALL	(HECMW_LOG_ERROR|HECMW_LOG_WARN|HECMW_LOG_INFO|HECMW_LOG_DEBUG)

/* オプション */

#define HECMW_LOG_PERROR	1

#define HECMW_LOG_OPTALL	(HECMW_LOG_PERROR)


extern int HECMW_openlog(const char *logfile, int loglv, int options);


extern int HECMW_closelog(int id);


extern int HECMW_vlog(int loglv, const char *fmt, va_list ap);


extern int HECMW_log(int loglv, const char *fmt , ...);


extern void HECMW_setloglv(int loglv);


extern void HECMW_log_set_enable(int from, int to, int true_or_false);

#endif
