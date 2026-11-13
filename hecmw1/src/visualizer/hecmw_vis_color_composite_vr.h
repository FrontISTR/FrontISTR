/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_VIS_COLOR_COMPOSITE_VR_H_INCLUDED
#define HECMW_VIS_COLOR_COMPOSITE_VR_H_INCLUDED

extern void compute_color_vr(
    int current_ijk[3], int color_mapping_style, double *interval_point,
    int transfer_function_style, double opa_value, int num_of_features,
    double *fea_point, int interval_mapping_num, int color_system_type,
    double k_ads[3], int r_level[3], double orig_xyz[3], double r_dxyz[3],
    double *var, double accum_rgba[4], double mincolor, double maxcolor,
    double in_point[3], double out_point[3], double tav_length);

#endif /* HECMW_VIS_COLOR_COMPOSITE_VR_H_INCLUDED */
