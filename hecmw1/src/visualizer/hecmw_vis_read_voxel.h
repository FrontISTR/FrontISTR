#ifndef HECMW_VIS_READ_VOXEL_H_INCLUDED
#define HECMW_VIS_READ_VOXEL_H_INCLUDED

int read_voxel_file(char *filename,  int n_voxel, double *voxel_dxyz, double *voxel_orig_xyz, int *level,
		int *voxel_n_neighbor_pe, int **voxel_neighbor_pe);

#endif /* HECMW_VIS_READ_VOXEL_H_INCLUDED */













