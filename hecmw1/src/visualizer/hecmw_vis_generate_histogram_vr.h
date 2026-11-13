/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_VIS_GENERATE_HISTOGRAM_VR_H_INCLUDED
#define HECMW_VIS_GENERATE_HISTOGRAM_VR_H_INCLUDED

#include "hecmw_util.h"

void find_color_minmax_vr(double *var, int *empty_flag, int nx, int ny, int nz,
                          double *mincolor, double *maxcolor);
void generate_histogram_graph_vr(double tmincolor, double tmaxcolor,
                                 double *var, int *empty_flag, int nx, int ny,
                                 int nz, int mynode, int pesize,
                                 HECMW_Comm VIS_COMM, int color_system_type);
void generate_interval_point_vr(double tmincolor, double tmaxcolor, double *var,
                                int *empty_flag, int nx, int ny, int nz,
                                int pesize, HECMW_Comm VIS_COMM,
                                double *interval_point);
void output_histogram_vr(double tmincolor, double tmaxcolor, double *var,
                         int *empty_flag, int nx, int ny, int nz, int mynode,
                         int pesize, HECMW_Comm VIS_COMM);
void find_minmax_vr(double *voxel_dxyz, double *voxel_orig_xyz, int mynode,
                    double range[6]);

#endif /* HECMW_VIS_GENERATE_HISTOGRAM_VR_H_INCLUDED */
