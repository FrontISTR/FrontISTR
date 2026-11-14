/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_VIS_MEM_UTIL_H_INCLUDED
#define HECMW_VIS_MEM_UTIL_H_INCLUDED

#include "hecmw_vis_SF_geom.h"

#if defined(__cplusplus) && __cplusplus >= 201103L
#define HECMW_VIS_NORETURN [[noreturn]]
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
#define HECMW_VIS_NORETURN _Noreturn
#elif defined(__GNUC__)
#define HECMW_VIS_NORETURN __attribute__((__noreturn__))
#else
#define HECMW_VIS_NORETURN
#endif

HECMW_VIS_NORETURN void HECMW_vis_memory_exit(char *var);
HECMW_VIS_NORETURN void HECMW_vis_print_exit(char *var);
void mfree(void *pointer);
Point *alloc_verts(int num);
Polygon *alloc_polygons(int num);

#endif /* HECMW_VIS_MEM_UTIL_H_INCLUDED */
