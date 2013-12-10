#ifndef HECMW_VIS_VOXEL_GEN_H_INCLUDED
#define HECMW_VIS_VOXEL_GEN_H_INCLUDED

#include "hecmw_util.h"

void voxel_gen(double range[6], double c_range[2], int nv[3], double *voxel_dxyz, double *voxel_orig_xyz, int *level,
		int *voxel_n_neighbor_pe, int **voxel_neighbor_pe, HECMW_Comm VIS_COMM, int vox_on, int display_range_on, double display_range[6]);

#endif /* HECMW_VIS_VOXEL_GEN_H_INCLUDED */











