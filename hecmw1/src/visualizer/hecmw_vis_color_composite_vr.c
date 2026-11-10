/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include "hecmw_vis_color_composite_vr.h"

#include <math.h>
#include "hecmw_vis_ray_trace.h"
#include "hecmw_malloc.h"

void compute_color_vr(int current_ijk[3], int color_mapping_style,
                      double *interval_point, int transfer_function_style,
                      double opa_value, int num_of_features, double *fea_point,
                      double view_point_d[3], int interval_mapping_num,
                      int color_system_type, int num_of_lights,
                      double *light_point, double k_ads[3], int r_level[3],
                      double orig_xyz[3], double r_dxyz[3], double *var,
                      double *grad_var, double accum_rgba[4], double mincolor,
                      double maxcolor, double grad_minmax[2],
                      double feap_minmax[2], double feai_minmax[2],
                      double dis_minmax[2], double *opa_table,
                      double in_point[3], double out_point[3],
                      double tav_length, int time_step, int print_flag) {
  int i, j, k;
  double a_current;
  double color[3];
  double value;
  double length, coff_i;
  int index;
  int m;
  double value2[8], vv[8 * 3], v;
  double dis[8], d;
  double t;
  double r, g, b;

  length =
      sqrt(SQR(out_point[0] - in_point[0]) + SQR(out_point[1] - in_point[1]) +
           SQR(out_point[2] - in_point[2]));
  coff_i = length / tav_length;
  index  = current_ijk[2] * r_level[0] * r_level[1] +
          current_ijk[1] * r_level[0] + current_ijk[0];

  /*---------------------start computing value of in_point
   * --------------------*/
  vv[0 * 3] = vv[4 * 3] = vv[7 * 3] = vv[3 * 3] =
      orig_xyz[0] + current_ijk[0] * r_dxyz[0];
  vv[1 * 3] = vv[5 * 3] = vv[6 * 3] = vv[2 * 3] =
      orig_xyz[0] + (current_ijk[0] + 1) * r_dxyz[0];
  vv[0 * 3 + 2] = vv[1 * 3 + 2] = vv[2 * 3 + 2] = vv[3 * 3 + 2] =
      orig_xyz[2] + current_ijk[2] * r_dxyz[2];
  vv[4 * 3 + 2] = vv[5 * 3 + 2] = vv[6 * 3 + 2] = vv[7 * 3 + 2] =
      orig_xyz[2] + (current_ijk[2] + 1) * r_dxyz[2];
  vv[0 * 3 + 1] = vv[1 * 3 + 1] = vv[5 * 3 + 1] = vv[4 * 3 + 1] =
      orig_xyz[1] + current_ijk[1] * r_dxyz[1];
  vv[3 * 3 + 1] = vv[2 * 3 + 1] = vv[6 * 3 + 1] = vv[7 * 3 + 1] =
      orig_xyz[1] + (current_ijk[1] + 1) * r_dxyz[1];

  for (m = 0; m < 8; m++) {
    i     = m % 2;
    j     = (m / 2) % 2;
    k     = (m / 2) / 2;
    index = (current_ijk[2] + k) * (r_level[0] + 1) * (r_level[1] + 1) +
            (current_ijk[1] + j) * (r_level[0] + 1) + current_ijk[0] + i;
    value2[k * 4 + j * 2 + i] = var[index];
  }

  for (i = 0; i < 8; i++) {
    dis[i] =
        sqrt(SQR(in_point[0] - vv[i * 3]) + SQR(in_point[1] - vv[i * 3 + 1]) +
             SQR(in_point[2] - vv[i * 3 + 2]));
  }
  d = 0.0;
  for (i = 0; i < 8; i++) d += 1.0 / (dis[i] + EPSILON);
  v      = 0.0;
  for (i = 0; i < 8; i++) {
    v += value2[i] / ((dis[i] + EPSILON) * d);
  }
  /*---------------------end computing value of in_point --------------------*/

  value = v;

  a_current = opa_value;
  if (transfer_function_style == 3) {
    a_current = opa_value;

    for (i = 0; i < num_of_features; i++) {
      t = fabs(value - fea_point[i * 3]);
      if (t < fea_point[i * 3 + 1])
        a_current += fea_point[i * 3 + 2] * (fea_point[i * 3 + 1] - t) /
                     fea_point[i * 3 + 1];
    }
  }

  if (transfer_function_style == 4) {
    a_current = opa_value;
    for (i = 0; i < num_of_features; i++) {
      if ((value >= fea_point[i * 3]) && (value <= fea_point[i * 3 + 1])) {
        a_current = fea_point[i * 3 + 2];
      }
    }
  }

  /*--------------map value to rgb -------------------   */
  if (color_mapping_style == 1) {
    if (fabs(maxcolor - mincolor) > EPSILON)
      value = (value - mincolor) / (maxcolor - mincolor);
  }

  if (color_mapping_style == 2) {
    mincolor = interval_point[0];
    maxcolor = interval_point[1];
    if (fabs(maxcolor - mincolor) > EPSILON)
      value = (value - mincolor) / (maxcolor - mincolor);
  }

  if (color_mapping_style == 3) {
    if (value < interval_point[0])
      value = 0.0;
    else if (value > interval_point[interval_mapping_num * 2])
      value = 1.0;
    else {
      for (i = 1; i < interval_mapping_num + 1; i++) {
        if ((value <= interval_point[i * 2]) &&
            (value > interval_point[(i - 1) * 2])) {
          value = (value - interval_point[(i - 1) * 2]) /
                      (interval_point[i * 2] - interval_point[(i - 1) * 2]) *
                      (interval_point[i * 2 + 1] -
                       interval_point[(i - 1) * 2 + 1]) +
                  interval_point[(i - 1) * 2 + 1];
          break;
        }
      }
    }
  }

  if (color_system_type == 1) {
    if (value < 0.0) value = 0.0;
    if (value > 1.0) value = 1.0;
    if (value <= 0.25) {
      r = 0.0;
      g = value * 4.0;
      b = 1.0;
    } else if ((value > 0.25) && (value <= 0.5)) {
      r = 0.0;
      g = 1.0;
      b = (0.5 - value) * 4.0;
    } else if ((value > 0.5) && (value <= 0.75)) {
      r = (value - 0.5) * 4.0;
      g = 1.0;
      b = 0.0;
    } else if (value > 0.75) {
      r = 1.0;
      g = (1.0 - value) * 4.0;
      b = 0.0;
    }

  }

  else if (color_system_type == 2) {
    if (value < 0.0) value = 0.0;
    if (value > 1.0) value = 1.0;
    if (value <= 0.2) {
      g = 0.0;
      b = 1.0;
      r = (0.2 - value) * 5.0;
    } else if ((value > 0.2) && (value <= 0.4)) {
      r = 0.0;
      b = 1.0;
      g = (value - 0.2) * 5.0;
    } else if ((value > 0.4) && (value <= 0.6)) {
      r = 0.0;
      g = 1.0;
      b = 1.0 - (value - 0.4) * 5.0;
    } else if ((value > 0.6) && (value <= 0.8)) {
      r = (value - 0.6) * 5.0;
      g = 1.0;
      b = 0.0;
    } else if (value > 0.0) {
      r = 1.0;
      g = 1.0 - (value - 0.8) * 5.0;
      b = 0.0;
    }
  } else if (color_system_type == 3) {
    r = g = b = value;
  }

  color[0] = r;
  color[1] = g;
  color[2] = b;
  /* --------------end of mapping value to rgb -----------------------  */

  r = 0.0;
  g = 0.0;
  b = 0.0;
/*	for(j=0;j<num_of_lights;j++) {
 */

  /*		r=color[0]*(ka+kd*cosalpha+ks*pow(costheta,6));
         g=color[1]*(ka+kd*cosalpha+ks*pow(costheta,6));
         b=color[2]*(ka+kd*cosalpha+ks*pow(costheta,6));
   */

  r += color[0] * k_ads[0] * coff_i;
  g += color[1] * k_ads[0] * coff_i;
  b += color[2] * k_ads[0] * coff_i;

  /*
 }
   */
  r *= a_current;
  g *= a_current;
  b *= a_current;
  if (accum_rgba[3] < 0.99) {
    accum_rgba[0] += r * (1.0 - accum_rgba[3]);
    accum_rgba[1] += g * (1.0 - accum_rgba[3]);
    accum_rgba[2] += b * (1.0 - accum_rgba[3]);
    accum_rgba[3] += a_current * (1.0 - accum_rgba[3]) * coff_i;
  }

  return;
}
