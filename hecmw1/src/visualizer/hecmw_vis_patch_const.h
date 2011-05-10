#ifndef HECMW_VIS_PATCH_CONST_H_INCLUDED
#define HECMW_VIS_PATCH_CONST_H_INCLUDED

#include "hecmw_vis_SF_geom.h"

/*----------------------------------------------------------------------
#     Subroutines in this file on isosurface generation by Marching Cubes is based
	  on the revision of Dr. Yuriko Takeshima's codes when she was working part time in RIST
#---------------------------------------------------------------------- */

int	merge_vol_iso(int, Cell *, double, Cube_polygons *,
		double, Cube_polygons *, int, int *,
		Point **, Point **, Point **, Polygon **);
int get_edge_index(int bound_index, int vert_index);
int get_vert_ident(int pindex, Cell *cell, double fvalue, int *sum_verts,
		Point **CS_verts_tail, Point **CS_verts_refer,
		Point **CS_verts_head, int bdflag);
int search_verts_table(Fgeom *pgeom, Point **CS_verts_refer,
		int bdflag, int table_no);
int add_verts_table(Point **CS_verts_tail, int table_no, int pident,
		double pfield, double cdata, Fgeom *pgeom,
		int bdflag);
int check_vertex(Cell *cell, double fvalue, Cube_polygons *cube_polys);

#endif /* HECMW_VIS_PATCH_CONST_H_INCLUDED */
