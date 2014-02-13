/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2007/06/29                                        *
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

 

#ifndef HECMW_UTIL_INCLUDED
#define HECMW_UTIL_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>
#include "hecmw_config.h"
#include "hecmw_init.h"
#include "hecmw_finalize.h"
#include "hecmw_malloc.h"
#include "hecmw_log.h"
#include "hecmw_msg.h"
#include "hecmw_lib_fc.h"
#include "hecmw_comm.h"
#include "hecmw_control.h"
#include "hecmw_error.h"
#include "hecmw_time.h"

#ifndef HECMW_SERIAL
#include "mpi.h"
#endif

#ifdef DEBUG
#define HECMW_DEBUG(args) \
	(HECMW_printerr("DEBUG: "), HECMW_printerr args, HECMW_printerr("\n"))
#ifdef HECMW_SERIAL
#define HECMW_assert(cond) assert(cond)
#else
#define HECMW_assert(cond) HECMW_assert_((cond) ? 1:0, #cond,__FILE__,__LINE__) 
#endif
#else
#define HECMW_DEBUG(args) ((void)0)
#define HECMW_assert(cond) ((void)0)
#endif

#define HECMW_check_condition( cond, isabort ) HECMW_check_condition_( (cond)?1:0, #cond, isabort, __FILE__, __LINE__ );


extern void HECMW_fprintf(FILE *fp, char *fmt, ...);


extern void HECMW_printerr(char *fmt, ...);


extern char *HECMW_get_date(void);


extern void HECMW_assert_(int cond, char *cond_str, char *file, int line);


extern int HECMW_check_condition_( int cond, char *cond_str, int isabort, char *file, int line );


extern void HECMW_abort(HECMW_Comm comm);



extern char *HECMW_toupper(char *s);

extern char *HECMW_tolower(char *s);

 
extern void HECMW_print_error(void);

 
extern void HECMW_print_msg(int loglv, int msgno, const char *fmt, ...);

 
extern void HECMW_print_vmsg(int loglv, int msgno, const char *fmt, va_list ap);


extern int HECMW_snprintf(char *msg, size_t size, const char *format, ...);

extern int HECMW_vsnprintf(char *msg, size_t size, const char *format, va_list ap);
#endif

