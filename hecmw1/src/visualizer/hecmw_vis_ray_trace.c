/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include "hecmw_vis_ray_trace.h"

#include <math.h>
#include "hecmw_vis_mem_util.h"
#include "hecmw_vis_color_composite_vr.h"

/*
void transform_face_node(int ijkn[3], int current_ijk[3], int face_sect[3],
double vv[3],
                                        int next_ijk[3]);
 */
static int find_out_point(double orig_xyz[3], double r_dxyz[3], int ijk[3],
                          double in_point[3], double ray_direction[3],
                          int face_sect[3], double out_point[3]);
static void find_next_cell(int ijkn[3], int current_ijk[3], int face_sect[3],
                           double vv[3], int next_ijk[3]);

static int find_first_inter_point(double point_o[3], double view_p[3],
                                  double ray_direction[3], double orig_xyz[3],
                                  double dxyz[3], int r_level[3],
                                  double first_p[3]) {
  int i, j, mincomp;
  double minmax[6], t[6], mint;

  for (i = 0; i < 3; i++) {
    minmax[i * 2]     = orig_xyz[i];
    minmax[i * 2 + 1] = orig_xyz[i] + dxyz[i];
  }
  for (i = 0; i < 3; i++) {
    if (fabs(ray_direction[i]) < EPSILON) {
      t[i * 2]     = 1.0E+17;
      t[i * 2 + 1] = 1.0E+17;
    } else {
      t[i * 2]     = (minmax[i * 2] - view_p[i]) / ray_direction[i];
      t[i * 2 + 1] = (minmax[i * 2 + 1] - view_p[i]) / ray_direction[i];
    }
  }

  mint    = 1.0E+17;
  mincomp = -1;
  for (i = 0; i < 6; i++) {
    if ((mint > t[i]) && (t[i] > EPSILON)) {
      for (j = 0; j < 3; j++) first_p[j] = view_p[j] + t[i] * ray_direction[j];
      if ((first_p[0] >= minmax[0] - EPSILON) &&
          (first_p[0] <= minmax[1] + EPSILON) &&
          (first_p[1] >= minmax[2] - EPSILON) &&
          (first_p[1] <= minmax[3] + EPSILON) &&
          (first_p[2] >= minmax[4] - EPSILON) &&
          (first_p[2] <= minmax[5] + EPSILON)) {
        mint    = t[i];
        mincomp = i;
      }
    }
  }
  if ((mincomp >= 0) && (mincomp <= 5)) {
    for (i = 0; i < 3; i++) first_p[i] = view_p[i] + mint * ray_direction[i];
  }
  /*	for(i=0;i<3;i++) {
                  if(fabs(out_point[i]-minmax[i*2])<EPSILON)
                          face_sect[i]=i*2;
                  else if(fabs(out_point[i]-minmax[i*2+1])<EPSILON)
                          face_sect[i]=i*2+1;
          }
   */

  return (mincomp);
}

int find_first_inter(double point_o[3], double view_point_d[3], int r_level[3],
                     double orig_xyz[3], double dxyz[3], double r_dxyz[3],
                     double ray_direction[3], double first_p[3], int ijk[3]) {
  int i;
  int ini_face, intersection, con1, next_ijk[3], face_sect[3];
  double out_point[3];

  ini_face = find_first_inter_point(point_o, view_point_d, ray_direction,
                                    orig_xyz, dxyz, r_level, first_p);
  if (ini_face == -1) { /*this ray has no intersection with dataset */
    intersection = 0;
  } else {
    intersection = 1;
    for (i   = 0; i < 3; i++)
      ijk[i] = (int)((first_p[i] - orig_xyz[i] - EPSILON) / r_dxyz[i]);
    if ((ijk[0] < 0) || (ijk[0] > r_level[0] - 1) || (ijk[1] < 0) ||
        (ijk[1] > r_level[1] - 1) || (ijk[2] < 0) || (ijk[2] > r_level[2] - 1))
      HECMW_vis_print_exit(
          "ERROR: HEC-MW-VIS-E2006:There is something wrong on finding the "
          "first point");
    for (i = 0; i < 3; i++) out_point[i] = 0.0;
    con1 = find_out_point(orig_xyz, r_dxyz, ijk, first_p, ray_direction,
                          face_sect, out_point);
    if (con1 == 0) {
      /* start voxel maybe at the neighbour voxel*/
      find_next_cell(r_level, ijk, face_sect, ray_direction, next_ijk);
      if ((next_ijk[0] < 0) || (next_ijk[0] > r_level[0] - 1) ||
          (next_ijk[1] < 0) || (next_ijk[1] > r_level[1] - 1) ||
          (next_ijk[2] < 0) || (next_ijk[2] > r_level[2] - 1))
        intersection = 0;

      for (i = 0; i < 3; i++) ijk[i] = next_ijk[i];
    }
  }
  return intersection;
}

static int find_out_point(double orig_xyz[3], double r_dxyz[3], int ijk[3],
                          double in_point[3], double ray_direction[3],
                          int face_sect[3], double out_point[3]) {
  int i, j, mincomp, intersection;
  double minmax[6], t[3], mint;
  int i1, j1;
  for (i = 0; i < 3; i++) face_sect[i] = -1;
  for (i = 0; i < 6; i++) {
    i1        = i / 2;
    j1        = i % 2;
    minmax[i] = orig_xyz[i1] + (ijk[i1] + j1) * r_dxyz[i1];
  }
  for (i = 0; i < 3; i++) {
    if (ray_direction[i] > 0)
      t[i] = (minmax[i * 2 + 1] - in_point[i]) / (ray_direction[i] + EPSILON);
    else if (ray_direction[i] < 0)
      t[i] = (in_point[i] - minmax[i * 2]) / (-ray_direction[i] + EPSILON);
  }

  mint    = 1.0E+17;
  mincomp = -1;
  for (i = 0; i < 3; i++) {
    if ((mint > t[i]) && (t[i] > EPSILON)) {
      for (j         = 0; j < 3; j++)
        out_point[j] = in_point[j] + t[i] * ray_direction[j];
      if ((out_point[0] >= minmax[0] - EPSILON) &&
          (out_point[0] <= minmax[1] + EPSILON) &&
          (out_point[1] >= minmax[2] - EPSILON) &&
          (out_point[1] <= minmax[3] + EPSILON) &&
          (out_point[2] >= minmax[4] - EPSILON) &&
          (out_point[2] <= minmax[5] + EPSILON)) {
        mint    = t[i];
        mincomp = i;
      }
    }
  }
  if ((mincomp >= 0) && (mincomp <= 2)) {
    intersection = 1;
    for (i         = 0; i < 3; i++)
      out_point[i] = in_point[i] + mint * ray_direction[i];
  } else {
    intersection = 0;
    for (i = 0; i < 3; i++) out_point[i] = in_point[i];
  }
  for (i = 0; i < 3; i++) {
    if (fabs(out_point[i] - minmax[i * 2]) < EPSILON)
      face_sect[i] = i * 2;
    else if (fabs(out_point[i] - minmax[i * 2 + 1]) < EPSILON)
      face_sect[i] = i * 2 + 1;
  }

  return (intersection);
}

static void find_next_cell(int ijkn[3], int current_ijk[3], int face_sect[3],
                           double vv[3], int next_ijk[3])
/*get the code of next cell*/
{
  int j, m;
  int face1_sect[3];

  for (j = 0; j < 3; j++) face1_sect[j] = face_sect[j];
  for (j = 0; j < 3; j++) next_ijk[j] = current_ijk[j];
  for (m = 0; m < 3; m++) {
    if (face_sect[m] != -1) {
      if (face_sect[m] == 0) {
        if (vv[0] < 0) {
          next_ijk[0]   = current_ijk[0] - 1;
          face1_sect[m] = 1;
          if (next_ijk[0] < 0) break;
        }
      } else if (face_sect[m] == 1) {
        if (vv[0] > 0) {
          next_ijk[0]   = current_ijk[0] + 1;
          face1_sect[m] = 0;
          if (next_ijk[0] > ijkn[0] - 1) break;
        }
      } else if (face_sect[m] == 2) {
        if (vv[1] < 0) {
          next_ijk[1]   = current_ijk[1] - 1;
          face1_sect[m] = 3;
          if (next_ijk[1] < 0) break;
        }
      } else if (face_sect[m] == 3) {
        if (vv[1] > 0) {
          next_ijk[1]   = current_ijk[1] + 1;
          face1_sect[m] = 2;
          if (next_ijk[1] > ijkn[1] - 1) break;
        }
      } else if (face_sect[m] == 4) {
        if (vv[2] < 0) {
          next_ijk[2]   = current_ijk[2] - 1;
          face1_sect[m] = 5;
          if (next_ijk[2] < 0) break;
        }
      } else if (face_sect[m] == 5) {
        if (vv[2] > 0) {
          next_ijk[2]   = current_ijk[2] + 1;
          face1_sect[m] = 4;
          if (next_ijk[2] > ijkn[2] - 1) break;
        }
      }
    }
  }

  for (m = 0; m < 3; m++) face_sect[m] = face1_sect[m];

  return;
}

static int find_next_point(double orig_xyz[3], double r_dxyz[3], int r_level[3],
                           int current_ijk[3], double in_point[3],
                           double out_point[3], double ray_direction[3],
                           int next_ijk[3]) {
  int flag, face_sect[3];

  flag = 1;
  flag = find_out_point(orig_xyz, r_dxyz, current_ijk, in_point, ray_direction,
                        face_sect, out_point);
  /*	con1=find_out_point(vd, current_ijk, in_point,ray_direction,face_sect,out_point);
	 */ if (flag == 0)
    HECMW_vis_print_exit(
        "ERROR: HEC-MW-VIS-E2007:There is some problem in finding intersection "
        "point");

  /*	    for(i=0;i<3;i++)
                 mid_point[i]=(out_point[i]+in_point[i])/2.0;
         for(i=0;i<3;i++)
         current_ijk[i]=(int)((mid_point[i]-vd->xyz0[i])/vd->dxyz[i]);
   */
  if (flag == 1) {
    find_next_cell(r_level, current_ijk, face_sect, ray_direction, next_ijk);
    if ((next_ijk[0] < 0) || (next_ijk[0] > r_level[0] - 1) ||
        (next_ijk[1] < 0) || (next_ijk[1] > r_level[1] - 1) ||
        (next_ijk[2] < 0) || (next_ijk[2] > r_level[2] - 1))
      flag = -1;
  }
  return (flag);
}

void ray_trace(int remove_0_display_on, int color_mapping_style,
               double *interval_point, int transfer_function_style,
               double opa_value, int num_of_features, double *fea_point,
               double view_point_d[3], int interval_mapping_num,
               int color_system_type, int num_of_lights, double *light_point,
               double k_ads[3], double orig_xyz[3], double dxyz[3],
               double r_dxyz[3], int r_level[3], int *empty_flag, double *var,
               double *grad_var, double first_p[3], int first_ijk[3],
               double ray_direction[3], double mincolor, double maxcolor,
               double accum_rgba[4], double grad_minmax[2],
               double feap_minmax[2], double feai_minmax[2],
               double dis_minmax[2], double *opa_table, double tav_length,
               int time_step, int test_i, int test_j)

/*void ray_trace(Parameter_vr *vr, VR_data *vd, Tree_pointer *root_tree,
   Tree_pointer *voxel_p,
                           double point_o[3], double ray_direction[3], double
   mincolor, double maxcolor,
                           int connect[8][6], int local_face[8][6], double
   accum_rgba[4], double grad_minmax[2],
                           double feap_minmax[2], double feai_minmax[2],
                           double dis_minmax[2], double *opa_table, double
   tav_length, int time_step)
 */
{
  int i, j;
  double in_point[3], out_point[3];
  int flag, single_empty_flag, value_flag, index, current_ijk[3], next_ijk[3];
  int print_flag;

  for (i = 0; i < 4; i++) accum_rgba[i] = 0.0;
  flag                                  = 1;
  for (i = 0; i < 3; i++) in_point[i] = first_p[i];
  for (i = 0; i < 3; i++) current_ijk[i] = first_ijk[i];
  while ((accum_rgba[3] < 0.99) && (flag == 1)) {
    /*
while((accum_rgba[3]<0.99) && (flag==1) && (accum_rgba[0]<0.99) &&
(accum_rgba[1]<0.99)
&& (accum_rgba[2]<0.99)) {
     */
    flag = find_next_point(orig_xyz, r_dxyz, r_level, current_ijk, in_point,
                           out_point, ray_direction, next_ijk);

    single_empty_flag = 1;
    value_flag        = 1;
    /*     for(k=0;k<2;k++)
     for(j=0;j<2;j++)
             for(i=0;i<2;i++) {

                     index=(current_ijk[2]+k)*(r_level[0]+1)*(r_level[1]+1)+(current_ijk[1]+j)*(r_level[0]+1)+
                             current_ijk[0]+i;
     */
    for (i = 0; i < 8; i++) {
      j     = (int)i % 4;
      index = (current_ijk[2] + (int)(i / 4)) * (r_level[0] + 1) *
                  (r_level[1] + 1) +
              (current_ijk[1] + (int)(j / 2)) * (r_level[0] + 1) +
              current_ijk[0] + j % 2;
      if (empty_flag[index] == 0) single_empty_flag = 0;
      /*                        if(fabs(var[index])>EPSILON)
         value_flag=1;
       */
    }
    /*     value_flag=1;
if(remove_0_display_on==1) {
     value_flag=0;
for(i=0;i<8;i++) {
j
             if(fabs(var[index])>EPSILON) {
                    value_flag=1;
                    break;
                     }
             }
}

overlap_flag=1;
for(k=0;k<3;k++) {
     if(fabs(in_point[k]-(vd->xyz0[k]+vd->dxyz[k]))<EPSILON) {
             overlap_flag=0;
             break;
     }
}
     */

    /*
if((empty_flag==1) && (value_flag==1))
     */
    if ((single_empty_flag == 1) && (flag != 0) && (value_flag == 1))
      compute_color_vr(current_ijk, color_mapping_style, interval_point,
                       transfer_function_style, opa_value, num_of_features,
                       fea_point, view_point_d, interval_mapping_num,
                       color_system_type, num_of_lights, light_point, k_ads,
                       r_level, orig_xyz, r_dxyz, var, grad_var, accum_rgba,
                       mincolor, maxcolor, grad_minmax, feap_minmax,
                       feai_minmax, dis_minmax, opa_table, in_point, out_point,
                       tav_length, time_step, print_flag);
    if (flag == 1) {
      for (i = 0; i < 3; i++) in_point[i] = out_point[i];
      for (i = 0; i < 3; i++) current_ijk[i] = next_ijk[i];
    }
    /*	 if((vd->empty_flag[current_voxel->cell_id-vd->nxyz[0]*vd->nxyz[1]*vd->nxyz[2]]==1)
&& ((fabs(vd->var[current_voxel->cell_id])>EPSILON) ||
     (fabs(vd->grad_var[current_voxel->cell_id*3])>EPSILON) ||
     (fabs(vd->grad_var[current_voxel->cell_id*3+1])>EPSILON) ||
     (fabs(vd->grad_var[current_voxel->cell_id*3+2])>EPSILON)))
     */
    /*	  if((vd->empty_flag[current_voxel->cell_id-vd->nxyz[0]*vd->nxyz[1]*vd->nxyz[2]]==1)
       && (fabs(vd->var[current_voxel->cell_id])>EPSILON))
     */
  }
  return;
}
