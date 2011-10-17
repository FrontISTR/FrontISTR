#ifndef HECMW_VIS_UCD_TRANS_H_INCLUDED
#define HECMW_VIS_UCD_TRANS_H_INCLUDED

#include "hecmw_struct.h"
#include "hecmw_result.h"
#include "hecmw_util.h"
#include "hecmw_vis_ray_trace.h"

void transform_ucd_pvr(struct hecmwST_result_data *data, double *node1,  struct hecmwST_local_mesh *mesh,
		Parameter_vr *vr,int my_rank, int pe_size, HECMW_Comm VIS_COMM, double *voxel_dxyz, double *voxel_orig_xyz,
		int *level, int *voxel_n_neighbor_pe, int **voxel_neighbor_pe, int voxel_on, int display_range_on,
		int init_flag, int num_of_pvr);

#endif /* HECMW_VIS_UCD_TRANS_H_INCLUDED */






