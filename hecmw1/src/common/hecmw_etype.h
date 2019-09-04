/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_ETYPE
#define INC_HECMW_ETYPE

#include "hecmw_util.h"

extern int HECMW_get_etype_UTIL2HECMW(int etype);

extern int HECMW_get_etype_HECMW2UTIL(int etype);

extern int HECMW_get_etype_GeoFEM2HECMW(int etype);

extern int HECMW_get_max_node(int etype);

extern int HECMW_get_max_edge(int etype);

extern int HECMW_get_max_surf(int etype);

extern int HECMW_get_max_tsuf(int etype);

extern int HECMW_get_max_qsuf(int etype);

extern char *HECMW_get_ucd_label(int etype);

extern int HECMW_get_etype_class(int etype);

extern int HECMW_get_etype_shape(int etype);

extern int HECMW_get_etype_vtk_shape(int etype);

extern int HECMW_is_etype_rod(int etype);

extern int HECMW_is_etype_surface(int etype);

extern int HECMW_is_etype_solid(int etype);

extern int HECMW_is_etype_interface(int etype);

extern int HECMW_is_etype_beam(int etype);

extern int HECMW_is_etype_shell(int etype);

extern int HECMW_is_etype_link(int etype);

extern int HECMW_is_etype_33struct(int etype);

extern int HECMW_is_etype_truss(int etype);

extern int HECMW_is_etype_patch(int etype);

extern const int *HECMW_get_surf_nodes(int etype, int sid, int *surf_etype);

#endif
