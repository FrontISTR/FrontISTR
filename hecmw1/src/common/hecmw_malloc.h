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



#ifndef HECMW_MALLOC_INCLUDED
#define HECMW_MALLOC_INCLUDED

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef HECMW_MALLOC
#define HECMW_malloc(size) HECMW_malloc_(size, __FILE__, __LINE__)
#define HECMW_calloc(nmemb, size) HECMW_calloc_(nmemb, size, __FILE__, __LINE__)
#define HECMW_realloc(ptr, size) HECMW_realloc_(ptr, size, __FILE__, __LINE__)
#define HECMW_strdup(s) HECMW_strdup_(s, __FILE__, __LINE__)
#define HECMW_free(ptr) HECMW_free_(ptr, __FILE__, __LINE__)
#else
#define HECMW_malloc(size) malloc(size)
#define HECMW_calloc(nmemb, size) calloc(nmemb, size)
#define HECMW_realloc(ptr, size) realloc(ptr, size)
#define HECMW_strdup(s) strdup(s)
#define HECMW_free(ptr) free(ptr)
#endif



extern void *HECMW_malloc_(size_t size, char *file, int line);

extern void *HECMW_calloc_(size_t nmemb, size_t size, char *file, int line);


extern void *HECMW_realloc_(void *ptr, size_t size, char *file, int line);


extern char *HECMW_strdup_(const char *s, char *file, int line);


extern void HECMW_free_(void *ptr, char *file, int line);


extern int HECMW_check_memleak(void);


extern long HECMW_get_memsize(void);


extern void HECMW_set_autocheck_memleak(int flag);


extern int HECMW_list_meminfo(FILE *fp);

#endif
