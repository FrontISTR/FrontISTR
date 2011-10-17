#ifndef HECMW_VIS_MESH_EXTENT_H_INCLUDED
#define HECMW_VIS_MESH_EXTENT_H_INCLUDED

#include "hecmw_struct.h"
#include "hecmw_util.h"

int calc_extent(struct hecmwST_local_mesh *mesh, double *extent);
int calc_voxel_level(int n_voxel, struct hecmwST_local_mesh *mesh, double *voxel_dxyz, double *voxel_orig_xyz,
		double *extent, int *level, HECMW_Comm VIS_COMM);

#endif /* HECMW_VIS_MESH_EXTENT_H_INCLUDED */









