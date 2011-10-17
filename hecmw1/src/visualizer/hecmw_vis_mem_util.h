#ifndef HECMW_VIS_MEM_UTIL_H_INCLUDED
#define HECMW_VIS_MEM_UTIL_H_INCLUDED

#include "hecmw_vis_SF_geom.h"
void HECMW_vis_memory_exit(char *var);
void HECMW_vis_print_exit(char *var);
void mfree(void *pointer);
Point *alloc_verts(int num);
Polygon *alloc_polygons(int num);

#endif /* HECMW_VIS_MEM_UTIL_H_INCLUDED */









