#ifndef HECMW_VIS_NEW_REFINE_H_INCLUDED
#define HECMW_VIS_NEW_REFINE_H_INCLUDED

#include "hecmw_struct.h"
#include "hecmw_util.h"

void transform_face_node(int face, int node[4]);
void refinement(struct hecmwST_local_mesh *mesh, double *node1,
		int n_voxel, double *voxel_dxyz, double *voxel_orig_xyz, int *level,
		int *voxel_n_neighbor_pe, int **voxel_neighbor_pe, double *extent,
		int my_rank, HECMW_Comm VIS_COMM,  int *empty_flag, double *var);

#endif /* HECMW_VIS_NEW_REFINE_H_INCLUDED */







