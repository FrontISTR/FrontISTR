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

#include "hecmw_vis_read_voxel.h"

#include <stdio.h>
#include "hecmw_vis_mem_util.h"


int read_voxel_file(char *filename,  int n_voxel, double *voxel_dxyz, double *voxel_orig_xyz, int *level,
		int *voxel_n_neighbor_pe, int **voxel_neighbor_pe)
{
	FILE		*fp;
	int		i, j;
	int		ret;
	if ((fp = fopen(filename, "r")) == NULL)
		HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0010: Cannot open voxel file");
	for (j = 0; j < n_voxel; j++) {
		/* read parallel information */
		ret = fscanf(fp, "%lf %lf %lf", &voxel_orig_xyz[j*3], &voxel_orig_xyz[j*3+1], &voxel_orig_xyz[j*3+2]);
		if (ret != 3) HECMW_vis_print_exit("ERROR: voxel file format error\n");
		ret = fscanf(fp, "%lf %lf %lf", &voxel_dxyz[j*3], &voxel_dxyz[j*3+1], &voxel_dxyz[j*3+2]);
		if (ret != 3) HECMW_vis_print_exit("ERROR: voxel file format error\n");

		ret = fscanf(fp, "%d", &voxel_n_neighbor_pe[j]);
		if (ret != 1) HECMW_vis_print_exit("ERROR: voxel file format error\n");

		for (i = 0; i < voxel_n_neighbor_pe[j]; i++) {
			ret = fscanf(fp, "%d", &voxel_neighbor_pe[j][i]);
			if (ret != 1) HECMW_vis_print_exit("ERROR: voxel file format error\n");
		}

		level[j*3+0] = 0;
		level[j*3+1] = 0;
		level[j*3+2] = 0;
	}
	fclose(fp);

	return 1;
}
#ifdef later
int write_voxel_file(Voxel_data *vox)
{
	int	i, j;

	for (i = 0; i < vox->n_voxel; i++) {
		fprintf(stderr, "%d: %lf %lf %lf(%lf %lf %lf)\n", i, vox->info[i].dx,
				vox->info[i].dy ,vox->info[i].dz,
				vox->info[i].orig_x, vox->info[i].orig_y,
				vox->info[i].orig_z);
		for (j = 0; j < vox->info[i].n_neighbor_pe; j++) {
			fprintf(stderr, "%d ", vox->info[i].neighbor_pe[j]);
		}
		fprintf(stderr, "\n");
	}

	return 1;
}

#endif
