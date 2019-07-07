/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_VIS_PSF_RENDERING_H_INCLUDED
#define HECMW_VIS_PSF_RENDERING_H_INCLUDED

/*#define  ResSize 1
#define  PixelSize 1
#define  DIGN_PE  2
 */
#define BAR_WIDTH 10
#define NUM_CONTROL_PARAS2 49
#define SQR(x) (x) * (x)
#define PI 3.1415926

typedef struct _rendering_parameter_struct {
  int xr;
  int yr;
  int projection_style;
  int num_of_lights;

  double view_point_d[3];
  double screen_point[3];
  double up[3];
  double k_ads[3];
  int color_mapping_style;
  double *light_point;
  int interval_mapping_num;
  double *interval_point; /* 2:mincolor, maxcolor 3: interval_mapping_num*2
                             (value, mark_value) */
  int transfer_function_style;
  /*  1: constant    input: value
2: first-order derivatives  input: none
  3: feature points  input:  num_of_featurepoints, point[num]
  4: feature intervals  input: num_of_intervals  point[num*2]
  5: distance inverse
  6: distance proportional
  7: look-up table      input: name of the look-up table file
   */
  double opa_value;
  int num_of_features;
  double *fea_point;
  char name_lookup[128];

  int rotate_style;
  int color_mapping_bar_on;
  int scale_marking_on;
  int rotate_num_of_frames;
  int deform_num_of_frames;
  double background_color[3];
  double font_color[3];
  int color_system_type;
  double font_size;
  int color_bar_style;
  int fixed_range_on;
  double range_value[2];
  int num_of_scale;
  int mark_0_on;
  int histogram_on;
  int boundary_line_on;
  double isoline_color[3];
  int time_mark_on;
  int fixed_scale_mark;
  int smooth_shading;
  double fixed_mesh_range[6];
  double start_time;
  double time_interval;
} Parameter_rendering;
/*
Parameter_rendering *sr;
 */

#endif /* HECMW_VIS_PSF_RENDERING_H_INCLUDED */
