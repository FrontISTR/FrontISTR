#ifndef HECMW_VIS_SURFACE_COMPUTE_H_INCLUDED
#define HECMW_VIS_SURFACE_COMPUTE_H_INCLUDED

#include "hecmw_struct.h"
#include "hecmw_result.h"
#include "hecmw_util.h"
#include "hecmw_vis_SF_geom.h"

/*----------------------------------------------------------------------
#     Subroutines in this file on isosurface generation for hexahedra by Marching Cubes is based
	  on the revision of Dr. Yuriko Takeshima's codes when she was working part time in RIST
#---------------------------------------------------------------------- */


int HECMW_vis_surface_compute(Surface *sff, struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data,
		int *bdflag, int *sum_v, int *sum_t,int *tvertex, int *tpatch,
		double *minc, double *maxc, Result *result, int sf_i, int mynode,  HECMW_Comm VIS_COMM );
int HECMW_vis_chk_bounds(struct hecmwST_local_mesh *mesh, int *bdflag);
void find_isoline(int isomum, int sum_polys, double mincolor, double maxcolor, double *vcoord,
		int *plist, double *vcolor, Isohead *isohead);
int find_line_segment(double f[3][3], double c[3], double isocolor, double iso_p[6]);
void line_find(double isocolor, double c[3], Fgeom g[3], int k, Isohead *isohead);
void isoline_free(int isonum,Isohead *isohead);

#endif /* HECMW_VIS_SURFACE_COMPUTE_H_INCLUDED */
