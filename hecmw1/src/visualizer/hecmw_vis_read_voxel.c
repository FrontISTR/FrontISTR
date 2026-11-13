/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include "hecmw_vis_read_voxel.h"

#include <stdio.h>
#include "hecmw_vis_mem_util.h"

int read_voxel_file(char *filename, int n_voxel, double *voxel_dxyz,
                    double *voxel_orig_xyz, int *level,
                    int *voxel_n_neighbor_pe, int **voxel_neighbor_pe) {
  FILE *fp;
  int i, j;
  int ret;
  if ((fp = fopen(filename, "r")) == NULL)
    HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0010: Cannot open voxel file");
  for (j = 0; j < n_voxel; j++) {
    /* read parallel information */
    ret = fscanf(fp, "%lf %lf %lf", &voxel_orig_xyz[j * 3],
                 &voxel_orig_xyz[j * 3 + 1], &voxel_orig_xyz[j * 3 + 2]);
    if (ret != 3) HECMW_vis_print_exit("ERROR: voxel file format error\n");
    ret = fscanf(fp, "%lf %lf %lf", &voxel_dxyz[j * 3], &voxel_dxyz[j * 3 + 1],
                 &voxel_dxyz[j * 3 + 2]);
    if (ret != 3) HECMW_vis_print_exit("ERROR: voxel file format error\n");

    ret = fscanf(fp, "%d", &voxel_n_neighbor_pe[j]);
    if (ret != 1) HECMW_vis_print_exit("ERROR: voxel file format error\n");

    for (i = 0; i < voxel_n_neighbor_pe[j]; i++) {
      ret = fscanf(fp, "%d", &voxel_neighbor_pe[j][i]);
      if (ret != 1) HECMW_vis_print_exit("ERROR: voxel file format error\n");
    }

    level[j * 3 + 0] = 0;
    level[j * 3 + 1] = 0;
    level[j * 3 + 2] = 0;
  }
  fclose(fp);

  return 1;
}
