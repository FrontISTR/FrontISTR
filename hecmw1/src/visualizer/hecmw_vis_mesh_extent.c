/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : Visualization                                     *
 *                                                                     *
 *            Written by Li Chen (Univ. of Tokyo)                      *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/

#include "hecmw_vis_mesh_extent.h"

#include <stdlib.h>
#include "hecmw_malloc.h"


int calc_extent(struct hecmwST_local_mesh *mesh, double *extent)
{
	int		nodeID;
	int		i, j;

	double	x_min, x_max, y_min, y_max, z_min, z_max;
	double	x, y, z;

	for (i = 0; i < mesh->n_elem; i++) {
		nodeID = mesh->elem_node_item[mesh->elem_node_index[i]] - 1;
		x_min = x_max = mesh->node[nodeID*3];
		y_min = y_max = mesh->node[nodeID*3+1];
		z_min = z_max = mesh->node[nodeID*3+2];
		for (j = mesh->elem_node_index[i]+1; j < mesh->elem_node_index[i+1]; j++) {
			nodeID = mesh->elem_node_item[j] - 1;
			x = mesh->node[nodeID*3];
			y = mesh->node[nodeID*3+1];
			z = mesh->node[nodeID*3+2];
			if (x_min > x) x_min = x;
			if (x_max < x) x_max = x;
			if (y_min > y) y_min = y;
			if (y_max < y) y_max = y;
			if (z_min > z) z_min = z;
			if (z_max < z) z_max = z;
		}
		extent[i*6] = x_min;
		extent[i*6+1] = x_max;
		extent[i*6+2] = y_min;
		extent[i*6+3] = y_max;
		extent[i*6+4] = z_min;
		extent[i*6+5] = z_max;
	}

	return 1;
}


int calc_voxel_level(int n_voxel, struct hecmwST_local_mesh *mesh, double *voxel_dxyz, double *voxel_orig_xyz,
		double *extent, int *level, HECMW_Comm VIS_COMM)
{
	int	i, j;
	int	flag;
	double	*extent_xyz;
	double 	*result_xyz;
	int           pesize;
	int	rank;


	extent_xyz = (double *)HECMW_calloc(n_voxel*3, sizeof(double));
	result_xyz = (double *)HECMW_calloc(n_voxel*3, sizeof(double));

	for (i = 0; i < n_voxel*3; i++) {
		extent_xyz[i] = voxel_dxyz[i];
	}

	HECMW_Comm_rank(VIS_COMM, &rank);
	HECMW_Comm_size(VIS_COMM, &pesize);
	for (i = 0; i < mesh->n_elem; i++) {
		if(mesh->elem_type[i]<400) {
			for (j = 0; j < n_voxel; j++) {
				flag = 1;
				if ((extent[i*6+1] < voxel_orig_xyz[0])
						|| (voxel_orig_xyz[0]+voxel_dxyz[0]< extent[i*6])) {
					flag = 0;
				} else if ((extent[i*6+3] < voxel_orig_xyz[1])
						|| (voxel_orig_xyz[1]+voxel_dxyz[1] < extent[i*6+2])) {
					flag = 0;
				} else if ((extent[i*6+5] < voxel_orig_xyz[2])
						|| (voxel_orig_xyz[2]+voxel_dxyz[2] < extent[i*6+4])) {
					flag = 0;
				}
				if (flag) {
					if (extent_xyz[j*3] > (extent[i*6+1] - extent[i*6]))
						extent_xyz[j*3] = extent[i*6+1] - extent[i*6];
					if (extent_xyz[j*3+1] > (extent[i*6+3] - extent[i*6+2]))
						extent_xyz[j*3+1] = extent[i*6+3] - extent[i*6+2];
					if (extent_xyz[j*3+2] > (extent[i*6+5] - extent[i*6+4]))
						extent_xyz[j*3+2] = extent[i*6+5] - extent[i*6+4];
				}
			}
		}
	}
	/*  for (i = 0; i < n_voxel; i++) {
    fprintf(stderr, "Extent %d: %lf %lf %lf\n",
	    i, extent_x[i], extent_y[i], extent_z[i]);
  }
	 */
	HECMW_Barrier(VIS_COMM);
	if(pesize>1) {
		HECMW_Allreduce(extent_xyz, result_xyz, n_voxel*3, HECMW_DOUBLE, HECMW_MIN,
				VIS_COMM);

	}
	else {
		for(i=0;i<n_voxel;i++){
			result_xyz[i]=extent_xyz[i];
		}
	}
	HECMW_Barrier(VIS_COMM);
	/*
  for (i = 0; i < vox->n_voxel; i++) {
    if (result_x[i] != vox->info[i].dx) {
      vox->info[i].level[0]
	= (int)(log(vox->info[i].dx / result_x[i])/log(2.0));
    } else vox->info[i].level[0] = 0;
    if (result_y[i] != vox->info[i].dy) {
      vox->info[i].level[1]
	= (int)(log(vox->info[i].dy / result_y[i])/log(2.0));
    } else vox->info[i].level[1] = 0;
    if (result_z[i] != vox->info[i].dz) {
      vox->info[i].level[2]
	= (int)(log(vox->info[i].dz / result_z[i])/log(2.0));
    } else vox->info[i].level[2] = 0;
  }
	 */
	for (i = 0; i < n_voxel; i++) {
		if (result_xyz[i*3] !=voxel_dxyz[i*3]) {
			level[i*3]
			      = (int)(voxel_dxyz[i*3] / result_xyz[i*3]);
		} else level[i*3] = 0;
		if (result_xyz[i*3+1] != voxel_dxyz[i*3+1]) {
			level[i*3+1]
			      = (int)(voxel_dxyz[i*3+1] / result_xyz[i*3+1]);
		} else level[i*3+1] = 0;
		if (result_xyz[i*3+2] != voxel_dxyz[i*3+2]) {
			level[i*3+2]
			      = (int)(voxel_dxyz[i*3+2] / result_xyz[i*3+2]);
		} else level[i*3+2] = 0;
	}


	HECMW_free(extent_xyz);
	HECMW_free(result_xyz);

	return 1;
}
