/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include "hecmw_vis_rendering.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hecmw_vis_bmp.h"
#include "hecmw_vis_mem_util.h"
#include "hecmw_vis_comm_util.h"
#include "hecmw_vis_define_parameters.h"
#include "hecmw_vis_subimage_composite_sf.h"
#include "hecmw_vis_color_composite_sf.h"
#include "hecmw_vis_font_texture.h"
#include "hecmw_vis_generate_histogram_vr.h"
#include "hecmw_vis_generate_histogram_sf.h"
#include "hecmw_vis_surface_compute.h"
#include "hecmw_vis_color_mapping.h"
#include "hecmw_malloc.h"

double *image;
int ms, ns;

void HECMW_vis_rendering_surface(struct surface_module *sf,
                                 Parameter_rendering *sr,
                                 struct hecmwST_local_mesh *mesh,
                                 struct hecmwST_result_data *data, int tvertex,
                                 int tpatch, int *color_list, double *minvalue,
                                 double *maxvalue, Result *result,
                                 char *outfile, int stat_para[NUM_CONTROL_PSF],
                                 HECMW_Comm VIS_COMM, int timestep)

{
  char fname[128];
  int i, j, k, m, ii, iii, m1, m2, mmm, mmm1, jjj;
  int pesize, mynode;
  FILE *outfp;
  double view_n[3], nview_n, vertex[3 * 8];
  double point_s[3], point_o[3];
  double range[6], xd, yd;
  double coff_matrix[3][3], n_vertex[3 * 8], inv_matrix[3][3];
  double view_point[3], scr_area[4], x, y;
  double tmincolor, tmaxcolor, trange[6], svertex[3 * 8], sn_vertex[3 * 8],
      sscr_area[4];
  int intersection;
  double accum_rgba[4];
  /* double       *accum_opa; */
  int pix_num, pix0_num, pixn;
  double *n_subimage, *n_subopa, *subimage, *subopa;
  double *subimage1, *subopa1;
  HECMW_Status stat;
  int starti, startj, endi, endj, start_x, start_y;
  double minx, miny, minz, maxx, maxy, maxz, tminx, tminy, tminz, tmaxx, tmaxy,
      tmaxz;
  int num_img, base_x, base_y;
  double o_v_point[3], m_scr[4];
  double mincolor, maxcolor, *mivalue, *mavalue, color_iso[3], org_mincolor,
      org_maxcolor;

  BITMAPFILEHEADER header; /* File header */
  unsigned char r, g, b;
  int ri, gi, bi;
  BITMAPINFOHEADER info;
  int scale, scale_type;
  double *z_buffer, scr_f[4], f[3][3], n_f[3][3], inter_point[3], normal[9];
  double value, c_value[3], depth;
  double deltac, iso_p[6], isocolor, n_iso[6], pixel_d[2][2], slip_y;
  double *iso_value, *iso_rgb;
  int pixel_i[2][2], flag, tmp_i, last_j, next_j;
  int normal_flag;
  double *v_normal;
  double *p_normal;
  int *v_num;

  double f_f[4], norm, v_f[9];
  int end_timestep;
  double start_time, time_interval;
  int len_step;
  double f_disp[3][3];
  int deform_flag, count_dash_num;

  HECMW_Comm_rank(VIS_COMM, &mynode);
  HECMW_Comm_size(VIS_COMM, &pesize);
  if (mynode == 0) {
    /*      t1=HECMW_Wtime();
fprintf(stderr, "Start rendering patches\n");
     */
  }
  /*    start_time=0.0;
time_interval=0.01;
start_timestep=1;
end_timestep=1000;
   */
  mivalue = (double *)HECMW_calloc(data->nn_component, sizeof(double));
  mavalue = (double *)HECMW_calloc(data->nn_component, sizeof(double));
  if ((mivalue == NULL) || (mavalue == NULL))
    HECMW_vis_memory_exit("mivalue, mavalue");
  if (pesize > 1) {
    HECMW_Allreduce(minvalue, mivalue, data->nn_component, HECMW_DOUBLE,
                    HECMW_MIN, VIS_COMM);
    HECMW_Allreduce(maxvalue, mavalue, data->nn_component, HECMW_DOUBLE,
                    HECMW_MAX, VIS_COMM);
  } else {
    for (i = 0; i < data->nn_component; i++) {
      mivalue[i] = minvalue[i];
      mavalue[i] = maxvalue[i];
    }
  }

  for (i = 0; i < data->nn_component; i++) {
    if (color_list[i] == 1) {
      org_mincolor = tmincolor = mivalue[i];
      org_maxcolor = tmaxcolor = mavalue[i];
      if (mynode == 0)
        fprintf(stderr,
                "For data attribute %s, the minimum value is %lf, the maximum "
                "value is %lf\n",
                data->node_label[i], tmincolor, tmaxcolor);
      if (sr->fixed_range_on == 1) {
        mivalue[i] = sr->range_value[0];
        mavalue[i] = sr->range_value[1];
        tmincolor  = sr->range_value[0];
        tmaxcolor  = sr->range_value[1];
      }
    }
  }

  if (sr->histogram_on == 1)
    output_histogram_sf(sf, color_list, data, mivalue, mavalue, result, mynode,
                        pesize, VIS_COMM);
  if (sr->histogram_on == 2)
    generate_histogram_graph_sf(sf, color_list, data, mivalue, mavalue, result,
                                mynode, pesize, VIS_COMM,
                                sr->color_system_type);
  if (sr->color_mapping_style == 4) {
    sr->interval_mapping_num = 10;
    sr->interval_point       = (double *)HECMW_calloc(22, sizeof(double));
    generate_interval_point_sf(sf, color_list, data, mivalue, mavalue, result,
                               mynode, pesize, VIS_COMM, sr->interval_point);
  }
  if (stat_para[70] == 0) {
    if (sr->boundary_line_on == 1)
      find_minmax_sf(mesh, mynode, range);
    else
      find_patch_minmax_sf(result, sf, range);
    for (i = 0; i < 6; i++) trange[i] = 0.0;
    minx                              = range[0];
    maxx                              = range[1];
    miny                              = range[2];
    maxy                              = range[3];
    minz                              = range[4];
    maxz                              = range[5];
    HECMW_Barrier(VIS_COMM);
    if (pesize > 1) {
      HECMW_Allreduce(&minx, &tminx, 1, HECMW_DOUBLE, HECMW_MIN, VIS_COMM);

      HECMW_Allreduce(&maxx, &tmaxx, 1, HECMW_DOUBLE, HECMW_MAX, VIS_COMM);
      HECMW_Allreduce(&miny, &tminy, 1, HECMW_DOUBLE, HECMW_MIN, VIS_COMM);
      HECMW_Allreduce(&maxy, &tmaxy, 1, HECMW_DOUBLE, HECMW_MAX, VIS_COMM);
      HECMW_Allreduce(&minz, &tminz, 1, HECMW_DOUBLE, HECMW_MIN, VIS_COMM);
      HECMW_Allreduce(&maxz, &tmaxz, 1, HECMW_DOUBLE, HECMW_MAX, VIS_COMM);
    } else {
      tminx = minx;
      tmaxx = maxx;
      tminy = miny;
      tmaxy = maxy;
      tminz = minz;
      tmaxz = maxz;
    }
    fprintf(stderr,
            "trange: minx=%lf maxx=%lf miny=%lf maxy=%lf minz=%lf maxz=%lf\n",
            tminx, tmaxx, tminy, tmaxy, tminz, tmaxz);
  } else {
    find_patch_minmax_sf(result, sf, range);

    for (i = 0; i < 6; i++) trange[i] = sr->fixed_mesh_range[i];
    tminx                             = trange[0];
    tmaxx                             = trange[1];
    tminy                             = trange[2];
    tmaxy                             = trange[3];
    tminz                             = trange[4];
    tmaxz                             = trange[5];

    fprintf(stderr,
            "Fixed range is used for dynamic dataset: minx=%lf maxx=%lf "
            "miny=%lf maxy=%lf minz=%lf maxz=%lf\n",
            tminx, tmaxx, tminy, tmaxy, tminz, tmaxz);
  }
  deform_flag = 0;
  for (iii = 0; iii < sf[0].surface_style; iii++) {
    if (sf[iii + 1].deform_display_on == 1) {
      if (stat_para[69] == 0)
        get_deform_scale(sf, iii + 1, tmaxx - tminx, tmaxy - tminy,
                         tmaxz - tminz, mesh, data, pesize, VIS_COMM);
      else
        sf[iii + 1].disp_scale = sf[iii + 1].real_disp_scale;
      deform_flag              = 1;
    }
  }
  if (deform_flag == 1) {
    find_new_patch_minmax_sf(result, sf, range);
    if (stat_para[70] == 0) {
      for (i = 0; i < 6; i++) trange[i] = 0.0;
      minx                              = range[0];
      maxx                              = range[1];
      miny                              = range[2];
      maxy                              = range[3];
      minz                              = range[4];
      maxz                              = range[5];
      HECMW_Barrier(VIS_COMM);
      if (pesize > 1) {
        HECMW_Allreduce(&minx, &tminx, 1, HECMW_DOUBLE, HECMW_MIN, VIS_COMM);

        HECMW_Allreduce(&maxx, &tmaxx, 1, HECMW_DOUBLE, HECMW_MAX, VIS_COMM);
        HECMW_Allreduce(&miny, &tminy, 1, HECMW_DOUBLE, HECMW_MIN, VIS_COMM);
        HECMW_Allreduce(&maxy, &tmaxy, 1, HECMW_DOUBLE, HECMW_MAX, VIS_COMM);
        HECMW_Allreduce(&minz, &tminz, 1, HECMW_DOUBLE, HECMW_MIN, VIS_COMM);
        HECMW_Allreduce(&maxz, &tmaxz, 1, HECMW_DOUBLE, HECMW_MAX, VIS_COMM);
      } else {
        tminx = minx;
        tmaxx = maxx;
        tminy = miny;
        tmaxy = maxy;
        tminz = minz;
        tmaxz = maxz;
      }
    }
  }

  if (stat_para[26] == 0) {
    sr->light_point = (double *)HECMW_calloc(3, sizeof(double));
    if (sr->light_point == NULL) HECMW_vis_memory_exit("sr: light_point");

    sr->light_point[0] = (tminx + tmaxx) / 2.0;
    sr->light_point[1] = tmaxy + 0.1 * (tmaxy - tminy);
    sr->light_point[2] = tmaxz + (tmaxz - tminz) * 2.0;
  }
  if (stat_para[27] == 0) {
    sr->view_point_d[0] = (tminx + tmaxx) / 2.0;
    sr->view_point_d[1] = tmaxy + 1.5 * (tmaxy - tminy);
    sr->view_point_d[2] = tmaxz + 1.5 * (tmaxz - tminz);
  }
  if (stat_para[28] == 0) {
    sr->screen_point[0] = (tminx + tmaxx) / 2.0;
    sr->screen_point[1] = (tminy + tmaxy) / 2.0;
    sr->screen_point[2] = (tminz + tmaxz) / 2.0;
  }

  trange[0] = tminx;
  trange[1] = tmaxx;
  trange[2] = tminy;
  trange[3] = tmaxy;
  trange[4] = tminz;
  trange[5] = tmaxz;

  HECMW_Barrier(VIS_COMM);
  /*----------------------------------------
Define the sampling and projection parameters
   */
  /* First, define the equation of screen */
  if (mesh->elem_type[0] < 300) {
    sr->up[0]           = 0.0;
    sr->up[1]           = 1.0;
    sr->up[2]           = 0.0;
    sr->view_point_d[0] = (tminx + tmaxx) / 2.0;
    sr->view_point_d[1] = (tminy + tmaxy) / 2.0;
    sr->view_point_d[2] = -(tmaxx - tminx) * 1000000.0;
    sr->screen_point[0] = sr->view_point_d[0];
    sr->screen_point[1] = sr->view_point_d[1];
    sr->screen_point[2] = 0.0;
  }
  /*  if(mynode==0) {
fprintf(stderr, "up_direction=%lf %lf %lf\n", sr->up[0], sr->up[1], sr->up[2]);
  fprintf(stderr, "rotate_style=%d num_of_frames=%d\n", sr->rotate_style,
sr->rotate_num_of_frames);
}
   */
  if (mynode == 0) {
    fprintf(stderr,
            "viewpoint=%lf %lf %lf screen_point=%lf %lf %lfup=%lf %lf %lf\n",
            sr->view_point_d[0], sr->view_point_d[1], sr->view_point_d[2],
            sr->screen_point[0], sr->screen_point[1], sr->screen_point[2],
            sr->up[0], sr->up[1], sr->up[2]);
    fprintf(stderr, "rotate_style=%d num_of_frames=%d\n", sr->rotate_style,
            sr->rotate_num_of_frames);
  }
  if (sr->rotate_style == 0)
    num_img = 1;
  else if ((sr->rotate_style >= 1) && (sr->rotate_style <= 4)) {
    num_img = sr->rotate_num_of_frames;
    for (i                = 0; i < 3; i++)
      sr->screen_point[i] = (trange[i * 2] + trange[i * 2 + 1]) / 2.0;
    /*find the range of projection */
    for (i = 0; i < 3; i++) o_v_point[i] = sr->view_point_d[i];
    m_scr[0] = m_scr[2] = 1.0E+17;
    m_scr[1] = m_scr[3] = -1.0E+17;

    for (ii = 0; ii < num_img; ii++) {
      /*	   if(ii!=0)
       */
      view1_parameter_define(
          ii, sr->rotate_num_of_frames, sr->rotate_style, sr->view_point_d,
          sr->screen_point, sr->num_of_lights, sr->light_point, sr->up, trange);
      for (i      = 0; i < 3; i++)
        view_n[i] = sr->screen_point[i] - sr->view_point_d[i];
      nview_n     = sqrt(SQR(view_n[0]) + SQR(view_n[1]) + SQR(view_n[2]));
      if (fabs(nview_n) > EPSILON)
        for (i = 0; i < 3; i++) view_n[i] /= nview_n;
      transform_range_vertex(trange, vertex);

      get_frame_transform_matrix(sr->view_point_d, sr->screen_point, sr->up,
                                 coff_matrix);
      transform_frame(sr->screen_point, vertex, coff_matrix, n_vertex);
      for (i          = 0; i < 3; i++)
        view_point[i] = sr->view_point_d[i] - sr->screen_point[i];
      transform2_frame(coff_matrix, view_point);
      find_projection_range(view_point, n_vertex, scr_area);
      if (scr_area[0] < m_scr[0]) m_scr[0] = scr_area[0];
      if (scr_area[1] > m_scr[1]) m_scr[1] = scr_area[1];
      if (scr_area[2] < m_scr[2]) m_scr[2] = scr_area[2];
      if (scr_area[3] > m_scr[3]) m_scr[3] = scr_area[3];
    }

    for (i = 0; i < 3; i++) sr->view_point_d[i] = o_v_point[i];
  }

  for (ii = 0; ii < num_img; ii++) {
    /*	   if(ii!=0)
     */
    if (num_img > 1)
      if (mynode == 0) {
        fprintf(stderr, "Start rendering rotation frame %d\n", ii);
        fprintf(
            stderr,
            "viewpoint=%lf %lf %lf screen_point=%lf %lf %lf up=%lf %lf %lf\n",
            sr->view_point_d[0], sr->view_point_d[1], sr->view_point_d[2],
            sr->screen_point[0], sr->screen_point[1], sr->screen_point[2],
            sr->up[0], sr->up[1], sr->up[2]);
        fprintf(stderr, "rotate_style=%d num_of_frames=%d\n", sr->rotate_style,
                sr->rotate_num_of_frames);
      }
    if (sr->rotate_style > 0)
      view_parameter_define(ii, sr->rotate_num_of_frames, sr->rotate_style,
                            sr->view_point_d, sr->screen_point, sr->up,
                            sr->num_of_lights, sr->light_point, trange);
    for (i      = 0; i < 3; i++)
      view_n[i] = sr->screen_point[i] - sr->view_point_d[i];
    nview_n     = sqrt(SQR(view_n[0]) + SQR(view_n[1]) + SQR(view_n[2]));
    if (fabs(nview_n) > EPSILON)
      for (i = 0; i < 3; i++) view_n[i] /= nview_n;
    /*   for(i=0;i<3;i++)
p_screen[i]=view_n[i];
p_screen[3]=-view_n[0]*vr->screen_point[0]-view_n[1]*vr->screen_point[1]-
view_n[2]*vr->screen_point[2];
     */

    /* Second, find the projection of the dataset in each PE */
    transform_range_vertex(trange, vertex);
    if (mynode == 0)
      fprintf(stderr, "viewpoint=%lf %lf %lf screen_point=%lf %lf %lf\n",
              sr->view_point_d[0], sr->view_point_d[1], sr->view_point_d[2],
              sr->screen_point[0], sr->screen_point[1], sr->screen_point[2]);
    get_frame_transform_matrix(sr->view_point_d, sr->screen_point, sr->up,
                               coff_matrix);
    transform_frame(sr->screen_point, vertex, coff_matrix, n_vertex);
    /*  if(mynode==0){
for(i=0;i<3;i++)
for(j=0;j<3;j++)
}
     */
    for (i          = 0; i < 3; i++)
      view_point[i] = sr->view_point_d[i] - sr->screen_point[i];
    transform2_frame(coff_matrix, view_point);
    find_projection_range(view_point, n_vertex, scr_area);
    /* for multi-PEs, should find minx, maxx, miny, maxy in scr_area for all PE
     */
    if (ii == 0) {
      if ((sr->color_mapping_bar_on == 0) && (sr->scale_marking_on == 0) &&
          (sr->time_mark_on == 0)) {
        start_x = 10;
        start_y = 10;

        if (num_img == 1) {
          xd = (scr_area[1] - scr_area[0]) / (sr->xr - 20);
          yd = (scr_area[3] - scr_area[2]) / (sr->yr - 20);
        }
        if (num_img > 1) {
          xd = (m_scr[1] - m_scr[0]) / (sr->xr - 20);
          yd = (m_scr[3] - m_scr[2]) / (sr->yr - 20);
        }

        if (xd >= yd) {
          yd = xd;
          if (num_img == 1)
            sr->yr =
                (int)((int)(20 + (scr_area[3] - scr_area[2]) / yd) / 10) * 10 +
                10;
          else if (num_img > 1)
            sr->yr =
                (int)((int)(20 + (m_scr[3] - m_scr[2]) / yd) / 10) * 10 + 10;
        }
        if (xd < yd) {
          xd = yd;

          if (num_img == 1)
            sr->xr =
                (int)((int)(20 + (scr_area[1] - scr_area[0]) / xd) / 8) * 8 + 8;
          else if (num_img > 1)
            sr->xr = (int)((int)(20 + (m_scr[1] - m_scr[0]) / xd) / 8) * 8 + 8;
        }
      } else if ((sr->color_mapping_bar_on == 0) &&
                 (sr->scale_marking_on == 0) && (sr->time_mark_on == 1)) {
        char buf3[128];
        sprintf(buf3, "%d",
                (int)(start_time + time_interval * (double)end_timestep));
        len_step   = strlen(buf3);
        scale_type = 2;
        scale      = (int)sr->font_size;
        start_x    = 10;
        start_y    = 10 + 7 * scale + 5;
        if (sr->yr < start_y + 10)
          HECMW_vis_print_exit(
              "The resolution for yr is not enough, please specify a larger "
              "one");
        if (sr->xr < start_x + (len_step + 3 + 3) * scale)
          HECMW_vis_print_exit(
              "The resolution for xr is not enough, please specify a larger "
              "one ");
        if (num_img == 1) {
          xd = (scr_area[1] - scr_area[0]) / (sr->xr - (start_x + 10));
          yd = (scr_area[3] - scr_area[2]) / (sr->yr - (start_y + 10));
        } else if (num_img > 1) {
          xd = (m_scr[1] - m_scr[0]) / (sr->xr - (start_x + 10));
          yd = (m_scr[3] - m_scr[2]) / (sr->yr - (start_y + 10));
        }
        if (xd >= yd) {
          yd = xd;
          if (num_img == 1)
            sr->yr =
                (int)((int)(start_y + 10 + (scr_area[3] - scr_area[2]) / yd) /
                      10) *
                    10 +
                10;
          else if (num_img > 1)
            sr->yr =
                (int)((int)(start_y + 10 + (m_scr[3] - m_scr[2]) / yd) / 10) *
                    10 +
                10;
        }
        if (xd < yd) {
          xd = yd;

          if (num_img == 1)
            sr->xr =
                (int)((int)(20 + (scr_area[1] - scr_area[0]) / xd) / 8) * 8 + 8;
          else if (num_img > 1)
            sr->xr = (int)((int)(20 + (m_scr[1] - m_scr[0]) / xd) / 8) * 8 + 8;
        }
      }

      else if ((sr->color_mapping_bar_on == 1) && (sr->scale_marking_on == 0)) {
        start_x = 30;
        start_y = 10;
        if (num_img == 1) {
          xd = (scr_area[1] - scr_area[0]) / (sr->xr - 40);
          yd = (scr_area[3] - scr_area[2]) / (sr->yr - 20);
        } else if (num_img > 1) {
          xd = (m_scr[1] - m_scr[0]) / (sr->xr - 40);
          yd = (m_scr[3] - m_scr[2]) / (sr->yr - 20);
        }

        if (xd >= yd) {
          yd = xd;
          if (num_img == 1)
            sr->yr =
                (int)((int)(20 + (scr_area[3] - scr_area[2]) / yd) / 10) * 10 +
                10;
          else if (num_img > 1)
            sr->yr =
                (int)((int)(20 + (m_scr[3] - m_scr[2]) / yd) / 10) * 10 + 10;
        } else if (xd < yd) {
          xd = yd;
          if (num_img == 1)
            sr->xr =
                (int)((int)(40 + (scr_area[1] - scr_area[0]) / xd) / 8) * 8 + 8;
          else if (num_img > 1)
            sr->xr = (int)((int)(40 + (m_scr[1] - m_scr[0]) / xd) / 8) * 8 + 8;
        }
      } else if ((sr->color_mapping_bar_on == 1) &&
                 (sr->scale_marking_on == 1)) {
        scale = (int)sr->font_size;
        if ((sr->font_size - scale) < 0.5 - EPSILON)
          scale_type = 1;
        else
          scale_type = 2;
        if (scale_type == 1) {
          start_x                               = 20 + 45 * scale;
          if (sr->color_bar_style == 2) start_x = 20 + 45 * scale + 15;
        } else if (scale_type == 2) {
          start_x                               = 20 + 63 * scale;
          if (sr->color_bar_style == 2) start_x = 20 + 63 * scale + 15;
        }

        start_y = 10;
        if (sr->xr < start_x + 10) {
          fprintf(stderr,
                  "The image x_resolution cannot write such size charaters\n");
          HECMW_vis_print_exit(
              "Please reduce the font size or enlarge x_resolution and run "
              "again");
        }
        if (num_img == 1) {
          xd = (scr_area[1] - scr_area[0]) / (sr->xr - (start_x + 10));
          yd = (scr_area[3] - scr_area[2]) / (sr->yr - 20);
        } else if (num_img > 1) {
          xd = (m_scr[1] - m_scr[0]) / (sr->xr - (start_x + 10));
          yd = (m_scr[3] - m_scr[2]) / (sr->yr - 20);
        }

        if (xd >= yd) {
          yd = xd;
          if (num_img == 1)
            sr->yr =
                (int)((int)(20 + (scr_area[3] - scr_area[2]) / yd) / 10) * 10 +
                10;
          else if (num_img > 1)
            sr->yr =
                (int)((int)(20 + (m_scr[3] - m_scr[2]) / yd) / 10) * 10 + 10;
        } else if (xd < yd) {
          xd = yd;
          if (num_img == 1)
            sr->xr =
                (int)((int)(start_x + 10 + (scr_area[1] - scr_area[0]) / xd) /
                      8) *
                    8 +
                8;
          else if (num_img > 1)
            sr->xr =
                (int)((int)(start_x + 10 + (m_scr[1] - m_scr[0]) / xd) / 8) *
                    8 +
                8;
        }
      }
    }
    if (ii > 0) {
      start_x -= base_x;
      start_y -= base_y;
    }
    if (num_img > 1) {
      base_x = (int)((scr_area[0] - m_scr[0]) / xd);
      base_y = (int)((scr_area[2] - m_scr[2]) / yd);

      start_x += base_x;
      start_y += base_y;
    }
    /* find the subimage range for this PE */
    transform_range_vertex(range, svertex);
    /*   output_frame(vr, view_point, n_vertex, scr_area, xd, yd);
     */
    transform_frame(sr->screen_point, svertex, coff_matrix, sn_vertex);
    find_projection_range(view_point, sn_vertex, sscr_area);
    /*   starti=(int)((sscr_area[0]-scr_area[0]+EPSILON)/xd);
endi=(int)((sscr_area[1]-scr_area[0]-EPSILON)/xd)+1;
startj=(int)((sscr_area[2]-scr_area[2]+EPSILON)/yd);
endj=(int)((sscr_area[3]-scr_area[2]-EPSILON)/yd)+1;
if(starti<0) starti=0;
if(startj<0) startj=0;
if(endi>vr->xr) endi=vr->xr;
if(endj>vr->yr) endj=vr->yr;

starti=(int)((sscr_area[0]-scr_area[0])/xd)+start_x-1;
endi=(int)((sscr_area[1]-scr_area[0])/xd)+start_x+1;
startj=(int)((sscr_area[2]-scr_area[2])/yd)+start_y-1;
endj=(int)((sscr_area[3]-scr_area[2])/yd)+start_y+1;
     */

    HECMW_Barrier(VIS_COMM);
#ifdef later
    grad_minmax[0] = 1.0E17;
    grad_minmax[1] = -1.0E17;
    if (transfer_function_style == 2) {
    }
    feap_minmax[0] = 1.0E17;
    feap_minmax[1] = -1.0E17;
    if (transfer_function_style == 3) {
      find_feap_minmax(num_of_features, fea_point, tmincolor, tmaxcolor,
                       feap_minmax);
      HECMW_Barrier(VIS_COMM);
      HECMW_Allreduce(&feap_minmax[0], &tfeap_minmax[0], 1, HECMW_DOUBLE,
                      HECMW_MIN, VIS_COMM);
      HECMW_Allreduce(&feap_minmax[1], &tfeap_minmax[1], 1, HECMW_DOUBLE,
                      HECMW_MAX, VIS_COMM);
      feap_minmax[0] = tfeap_minmax[0];
      feap_minmax[1] = tfeap_minmax[1];
    }
    feai_minmax[0] = 1.0E17;
    feai_minmax[1] = -1.0E17;
    if (transfer_function_style == 4) {
      /*	   find_feai_minmax(num_of_features, fea_point, tmincolor,
tmaxcolor, feai_minmax);
HECMW_Barrier(VIS_COMM);
HECMW_Allreduce(&feai_minmax[0], &tfeai_minmax[0], 1, HECMW_DOUBLE, HECMW_MIN,
VIS_COMM);
HECMW_Allreduce(&feai_minmax[1], &tfeai_minmax[1], 1, HECMW_DOUBLE, HECMW_MAX,
VIS_COMM);
feai_minmax[0]=tfeai_minmax[0];  feai_minmax[1]=tfeai_minmax[1];
       */
    }
    dis_minmax[0] = 1.0E17;
    dis_minmax[1] = -1.0E17;
    if ((transfer_function_style == 5) || (transfer_function_style == 6))
      find_dis_minmax(view_point_d, vertex, dis_minmax);
#endif
    /* Build hierachical tree*/
    for (jjj = 0; jjj < sr->deform_num_of_frames; jjj++) {
      if (sr->deform_num_of_frames > 1)
        fprintf(stderr, "Start rendering deformation frame %d\n", jjj);
      image = (double *)HECMW_calloc(sr->xr * sr->yr * 3, sizeof(double));
      if (image == NULL) HECMW_vis_memory_exit("image");
      /*   accum_opa=(double *)HECMW_calloc(sr->xr*sr->yr, sizeof(double));
if(accum_opa==NULL) {
fprintf(stderr, "There is no enough memory: image\n");
exit(0);
}
       */
      /*}
       */
      /*---------new change for vec ---------*/
      for (j = 0; j < sr->yr * sr->xr * 3; j++) image[j] = 0.0;
      /*   for(j=0;j<sr->yr*sr->xr;j++)
accum_opa[j]=0.0;
       */
      /*   for(j=0;j<yr;j++)
for(i=0;i<xr;i++) {
image[(j*xr+i)*3] =image[(j*xr+i)*3+1] = image[(j*xr+i)*3+2]=0.0;
 accum_opa[j*xr+i]=0.0;
}
       */
      z_buffer = (double *)HECMW_calloc(sr->xr * sr->yr, sizeof(double));
      if (z_buffer == NULL) HECMW_vis_memory_exit("z_buffer");
      for (j = 0; j < sr->xr * sr->yr; j++) z_buffer[j] = 1.0E16;

      HECMW_Barrier(VIS_COMM);
      /*   t3=HECMW_Wtime();
if(mynode==0) {
fprintf(stderr, "The time for preprocessing of pvr is %lf, the total time up to
now is %lf\n", t3-t2, t3-t1);
t2=t3;
}
       */
      find_inverse_matrix(coff_matrix, inv_matrix);
      for (i = 0; i < data->nn_component; i++) {
        if (color_list[i] == 1) {
          mincolor = mivalue[i];
          maxcolor = mavalue[i];
        }
      }

      if (sr->color_mapping_style == 2) {
        mincolor = sr->interval_point[0];
        maxcolor = sr->interval_point[1];
      }
      if (sr->color_mapping_style == 3) {
        mincolor = sr->interval_point[0];
        maxcolor = sr->interval_point[2 * sr->interval_mapping_num];
      }
      for (iii = 0; iii < sf[0].surface_style; iii++) {
        if (sf[iii + 1].deform_display_on == 1) {
          if ((sf[iii + 1].initial_style == 2) ||
              (sf[iii + 1].initial_style == 3)) {
            if (sr->smooth_shading == 1) {
              v_normal = (double *)HECMW_calloc(result[iii].n_vertex * 3,
                                                sizeof(double));
              p_normal = (double *)HECMW_calloc(result[iii].n_patch * 4,
                                                sizeof(double));
              v_num = (int *)HECMW_calloc(result[iii].n_vertex, sizeof(int));
              if ((v_normal == NULL) || (v_num == NULL) || (p_normal == NULL))
                HECMW_vis_memory_exit("v_normal, p_normal");
              for (k = 0; k < result[iii].n_vertex; k++) v_num[k] = 0;
              for (k = 0; k < result[iii].n_vertex * 3; k++) v_normal[k] = 0.0;
              for (k = 0; k < result[iii].n_patch; k++) {
                for (i = 0; i < 3; i++) {
                  for (j = 0; j < 3; j++)
                    f[i][j] =
                        result[iii]
                            .vertex[(result[iii].patch[k * 3 + i] - 1) * 3 + j];
                }
                p_normal[k * 4 + 0] =
                    (f[1][1] - f[0][1]) * (f[2][2] - f[0][2]) -
                    (f[2][1] - f[0][1]) * (f[1][2] - f[0][2]);
                p_normal[k * 4 + 1] =
                    -(f[1][0] - f[0][0]) * (f[2][2] - f[0][2]) +
                    (f[2][0] - f[0][0]) * (f[1][2] - f[0][2]);
                p_normal[k * 4 + 2] =
                    (f[1][0] - f[0][0]) * (f[2][1] - f[0][1]) -
                    (f[2][0] - f[0][0]) * (f[1][1] - f[0][1]);
                norm = sqrt(p_normal[k * 4 + 0] * p_normal[k * 4 + 0] +
                            p_normal[k * 4 + 1] * p_normal[k * 4 + 1] +
                            p_normal[k * 4 + 2] * p_normal[k * 4 + 2]);
                if (fabs(norm) > EPSILON) {
                  for (j = 0; j < 3; j++) p_normal[k * 4 + j] /= norm;
                }
                p_normal[k * 4 + 3] = -p_normal[k * 4 + 0] * f[0][0] -
                                      p_normal[k * 4 + 1] * f[0][1] -
                                      p_normal[k * 4 + 2] * f[0][2];
                for (i = 0; i < 3; i++) {
                  v_num[result[iii].patch[k * 3 + i] - 1]++;
                  for (j = 0; j < 3; j++)
                    v_normal[(result[iii].patch[k * 3 + i] - 1) * 3 + j] +=
                        p_normal[k * 4 + j];
                }
              }
              for (k = 0; k < result[iii].n_vertex; k++) {
                if (v_num[k] != 0) {
                  for (j = 0; j < 3; j++) v_normal[k * 3 + j] /= v_num[k];
                }
              }
            }
            normal_flag = 1;
            for (k = 0; k < result[iii].n_patch; k++) {
              for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++)
                  f[i][j] =
                      result[iii]
                          .vertex[(result[iii].patch[k * 3 + i] - 1) * 3 + j];
                c_value[i] =
                    result[iii].color[result[iii].patch[k * 3 + i] - 1];
              }
              if (sr->smooth_shading == 1) {
                for (i = 0; i < 3; i++)
                  for (j = 0; j < 3; j++)
                    normal[i * 3 + j] =
                        v_normal[(result[iii].patch[k * 3 + i] - 1) * 3 + j];
              }

              transform_frame3(sr->screen_point, f, coff_matrix, n_f);
              find_projection_range2(view_point, n_f, scr_f);
              starti = (int)((scr_f[0] - scr_area[0]) / xd) + start_x - 1;
              endi   = (int)((scr_f[1] - scr_area[0]) / xd) + start_x + 1;
              startj = (int)((scr_f[2] - scr_area[2]) / yd) + start_y - 1;
              endj   = (int)((scr_f[3] - scr_area[2]) / yd) + start_y + 1;
              for (j = startj; j < endj; j++)
                for (i = starti; i < endi; i++) {
                  x = xd * (i - start_x + 0.5) + scr_area[0];
                  y = yd * (j - start_y + 0.5) + scr_area[2];
                  /*transform to the original frame */
                  point_s[0] = x;
                  point_s[1] = y;
                  point_s[2] = 0.0;
                  tranverse_transform(sr->screen_point, point_s, inv_matrix,
                                      point_o);

                  intersection = find_surface_point(
                      f, point_o, sr->view_point_d, inter_point, f_f, v_f,
                      c_value, &value, normal, normal_flag, sr->smooth_shading);
                  if (intersection == 1) {
                    depth = sqrt(SQR(inter_point[0] - sr->view_point_d[0]) +
                                 SQR(inter_point[1] - sr->view_point_d[1]) +
                                 SQR(inter_point[2] - sr->view_point_d[2]));
                    if (depth < z_buffer[j * sr->xr + i]) {
                      z_buffer[j * sr->xr + i] = depth;
                      if (sf[iii + 1].initial_style == 2) {
                        value = (mincolor + maxcolor) / 2.0;
                        compute_color_sf(
                            inter_point, value, normal, sr->color_mapping_style,
                            sr->interval_point, sr->view_point_d,
                            sr->interval_mapping_num, 3, sr->num_of_lights,
                            sr->light_point, sr->k_ads, accum_rgba, mincolor,
                            maxcolor, sf[iii + 1].display_method);
                      } else
                        compute_color_sf(
                            inter_point, value, normal, sr->color_mapping_style,
                            sr->interval_point, sr->view_point_d,
                            sr->interval_mapping_num, sr->color_system_type,
                            sr->num_of_lights, sr->light_point, sr->k_ads,
                            accum_rgba, mincolor, maxcolor,
                            sf[iii + 1].display_method);
                      image[(j * sr->xr + i) * 3]     = accum_rgba[0];
                      image[(j * sr->xr + i) * 3 + 1] = accum_rgba[1];
                      image[(j * sr->xr + i) * 3 + 2] = accum_rgba[2];
                    }
                  }
                }
            }
          } /* end if (initial_style=2 or 3) */
          if ((sf[iii + 1].deform_style == 2) ||
              (sf[iii + 1].deform_style == 3)) {
            normal_flag = 1;
            if (sr->smooth_shading == 1) {
              v_normal = (double *)HECMW_calloc(result[iii].n_vertex * 3,
                                                sizeof(double));
              p_normal = (double *)HECMW_calloc(result[iii].n_patch * 4,
                                                sizeof(double));
              v_num = (int *)HECMW_calloc(result[iii].n_vertex, sizeof(int));
              if ((v_normal == NULL) || (v_num == NULL) || (p_normal == NULL))
                HECMW_vis_memory_exit("v_normal_d, p_normal");
              for (k = 0; k < result[iii].n_vertex; k++) v_num[k] = 0;
              for (k = 0; k < result[iii].n_vertex * 3; k++) v_normal[k] = 0.0;
              for (k = 0; k < result[iii].n_patch; k++) {
                for (i = 0; i < 3; i++) {
                  if (stat_para[70] == 0) {
                    for (j = 0; j < 3; j++)
                      f[i][j] =
                          result[iii]
                              .vertex[(result[iii].patch[k * 3 + i] - 1) * 3 +
                                      j] +
                          result[iii]
                                  .disp[(result[iii].patch[k * 3 + i] - 1) * 3 +
                                        j] *
                              (sf[iii + 1].disp_scale /
                               (sr->deform_num_of_frames) * (jjj + 1));
                  } else {
                    for (j = 0; j < 3; j++)
                      f[i][j] =
                          result[iii]
                              .vertex[(result[iii].patch[k * 3 + i] - 1) * 3 +
                                      j] +
                          result[iii]
                                  .disp[(result[iii].patch[k * 3 + i] - 1) * 3 +
                                        j] *
                              sf[iii + 1].disp_scale;
                  }
                }
                p_normal[k * 4 + 0] =
                    (f[1][1] - f[0][1]) * (f[2][2] - f[0][2]) -
                    (f[2][1] - f[0][1]) * (f[1][2] - f[0][2]);
                p_normal[k * 4 + 1] =
                    -(f[1][0] - f[0][0]) * (f[2][2] - f[0][2]) +
                    (f[2][0] - f[0][0]) * (f[1][2] - f[0][2]);
                p_normal[k * 4 + 2] =
                    (f[1][0] - f[0][0]) * (f[2][1] - f[0][1]) -
                    (f[2][0] - f[0][0]) * (f[1][1] - f[0][1]);
                norm = sqrt(p_normal[k * 4 + 0] * p_normal[k * 4 + 0] +
                            p_normal[k * 4 + 1] * p_normal[k * 4 + 1] +
                            p_normal[k * 4 + 2] * p_normal[k * 4 + 2]);
                if (fabs(norm) > EPSILON) {
                  for (j = 0; j < 3; j++) p_normal[k * 4 + j] /= norm;
                }
                p_normal[k * 4 + 3] = -p_normal[k * 4 + 0] * f[0][0] -
                                      p_normal[k * 4 + 1] * f[0][1] -
                                      p_normal[k * 4 + 2] * f[0][2];
                for (i = 0; i < 3; i++) {
                  v_num[result[iii].patch[k * 3 + i] - 1]++;
                  for (j = 0; j < 3; j++)
                    v_normal[(result[iii].patch[k * 3 + i] - 1) * 3 + j] +=
                        p_normal[k * 4 + j];
                }
              }
              for (k = 0; k < result[iii].n_vertex; k++) {
                if (v_num[k] != 0) {
                  for (j = 0; j < 3; j++) v_normal[k * 3 + j] /= v_num[k];
                }
              }
            }

            for (k = 0; k < result[iii].n_patch; k++) {
              for (i = 0; i < 3; i++) {
                if (stat_para[70] == 0) {
                  for (j = 0; j < 3; j++)
                    f_disp[i][j] =
                        result[iii]
                            .vertex[(result[iii].patch[k * 3 + i] - 1) * 3 +
                                    j] +
                        result[iii]
                                .disp[(result[iii].patch[k * 3 + i] - 1) * 3 +
                                      j] *
                            (sf[iii + 1].disp_scale /
                             (sr->deform_num_of_frames) * (jjj + 1));
                } else {
                  for (j = 0; j < 3; j++)
                    f_disp[i][j] =
                        result[iii]
                            .vertex[(result[iii].patch[k * 3 + i] - 1) * 3 +
                                    j] +
                        result[iii]
                                .disp[(result[iii].patch[k * 3 + i] - 1) * 3 +
                                      j] *
                            sf[iii + 1].disp_scale;
                }

                c_value[i] =
                    result[iii].color[result[iii].patch[k * 3 + i] - 1];
              }
              if (sr->smooth_shading == 1) {
                for (i = 0; i < 3; i++)
                  for (j = 0; j < 3; j++)
                    normal[i * 3 + j] =
                        v_normal[(result[iii].patch[k * 3 + i] - 1) * 3 + j];
              }

              transform_frame3(sr->screen_point, f_disp, coff_matrix, n_f);
              find_projection_range2(view_point, n_f, scr_f);
              starti = (int)((scr_f[0] - scr_area[0]) / xd) + start_x - 1;
              endi   = (int)((scr_f[1] - scr_area[0]) / xd) + start_x + 1;
              startj = (int)((scr_f[2] - scr_area[2]) / yd) + start_y - 1;
              endj   = (int)((scr_f[3] - scr_area[2]) / yd) + start_y + 1;
              for (j = startj; j < endj; j++)
                for (i = starti; i < endi; i++) {
                  x = xd * (i - start_x + 0.5) + scr_area[0];
                  y = yd * (j - start_y + 0.5) + scr_area[2];
                  /*transform to the original frame */
                  point_s[0] = x;
                  point_s[1] = y;
                  point_s[2] = 0.0;
                  tranverse_transform(sr->screen_point, point_s, inv_matrix,
                                      point_o);

                  intersection = find_surface_point(
                      f_disp, point_o, sr->view_point_d, inter_point, f_f, v_f,
                      c_value, &value, normal, normal_flag, sr->smooth_shading);
                  if (intersection == 1) {
                    depth = sqrt(SQR(inter_point[0] - sr->view_point_d[0]) +
                                 SQR(inter_point[1] - sr->view_point_d[1]) +
                                 SQR(inter_point[2] - sr->view_point_d[2]));
                    if (depth < z_buffer[j * sr->xr + i]) {
                      z_buffer[j * sr->xr + i] = depth;

                      if (sf[iii + 1].deform_style == 2) {
                        value = (mincolor + maxcolor) / 2.0;
                        compute_color_sf(
                            inter_point, value, normal, sr->color_mapping_style,
                            sr->interval_point, sr->view_point_d,
                            sr->interval_mapping_num, 3, sr->num_of_lights,
                            sr->light_point, sr->k_ads, accum_rgba, mincolor,
                            maxcolor, sf[iii + 1].display_method);
                      } else
                        compute_color_sf(
                            inter_point, value, normal, sr->color_mapping_style,
                            sr->interval_point, sr->view_point_d,
                            sr->interval_mapping_num, sr->color_system_type,
                            sr->num_of_lights, sr->light_point, sr->k_ads,
                            accum_rgba, mincolor, maxcolor,
                            sf[iii + 1].display_method);
                      image[(j * sr->xr + i) * 3]     = accum_rgba[0];
                      image[(j * sr->xr + i) * 3 + 1] = accum_rgba[1];
                      image[(j * sr->xr + i) * 3 + 2] = accum_rgba[2];
                    }
                  }
                }
            }
          }
          if ((sf[iii + 1].initial_style == 1) ||
              (sf[iii + 1].initial_style == 4)) {
            normal_flag = 1;
            for (k = 0; k < result[iii].n_patch; k++) {
              for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++)
                  f[i][j] =
                      result[iii]
                          .vertex[(result[iii].patch[k * 3 + i] - 1) * 3 + j];
                c_value[i] =
                    result[iii].color[result[iii].patch[k * 3 + i] - 1];
              }
              for (mmm = 0; mmm < 3; mmm++) {
                mmm1                = mmm + 1;
                if (mmm1 == 3) mmm1 = 0;
                for (j = 0; j < 3; j++) {
                  iso_p[j]     = f[mmm][j];
                  iso_p[j + 3] = f[mmm1][j];
                }

                count_dash_num = -1;
                transform_frame4(sr->screen_point, iso_p, coff_matrix, n_iso);
                find_projection_range3(view_point, n_iso, pixel_d, iso_p);
                pixel_i[0][0] =
                    (int)((pixel_d[0][0] - scr_area[0]) / xd - 0.5) + start_x;
                pixel_i[0][1] =
                    (int)((pixel_d[0][1] - scr_area[2]) / yd - 0.5) + start_y;
                pixel_i[1][0] =
                    (int)((pixel_d[1][0] - scr_area[0]) / xd + 0.5) + start_x;
                pixel_i[1][1] =
                    (int)((pixel_d[1][1] - scr_area[2]) / yd + 0.5) + start_y;
                if (pixel_i[1][0] != pixel_i[0][0]) {
                  slip_y = (double)(pixel_i[1][1] - pixel_i[0][1]) /
                           (double)(pixel_i[1][0] - pixel_i[0][0]);
                  last_j = pixel_i[0][1];
                  for (i = pixel_i[0][0]; i <= pixel_i[1][0]; i++) {
                    next_j = (int)(pixel_i[0][1] +
                                   slip_y * (double)(i - pixel_i[0][0]));
                    if (last_j <= next_j) {
                      for (j = last_j; j <= next_j; j++) {
                        x          = xd * (i - start_x + 0.5) + scr_area[0];
                        y          = yd * (j - start_y + 0.5) + scr_area[2];
                        point_s[0] = x;
                        point_s[1] = y;
                        point_s[2] = 0.0;
                        tranverse_transform(sr->screen_point, point_s,
                                            inv_matrix, point_o);

                        intersection =
                            find_point_depth(f, point_o, sr->view_point_d, f_f,
                                             inter_point, normal_flag);
                        if (intersection == 1)
                          depth =
                              sqrt(SQR(inter_point[0] - sr->view_point_d[0]) +
                                   SQR(inter_point[1] - sr->view_point_d[1]) +
                                   SQR(inter_point[2] - sr->view_point_d[2]));
                        else
                          depth = 1.0E17;
                        if (depth <= z_buffer[j * sr->xr + i] + EPSILON) {
                          if (sf[iii + 1].initial_style == 4) {
                            count_dash_num++;
                            if ((count_dash_num % 4) < 1) {
                              image[(j * sr->xr + i) * 3] =
                                  sf[iii + 1].initial_line_color[0];
                              image[(j * sr->xr + i) * 3 + 1] =
                                  sf[iii + 1].initial_line_color[1];
                              image[(j * sr->xr + i) * 3 + 2] =
                                  sf[iii + 1].initial_line_color[2];
                              z_buffer[j * sr->xr + i] = depth;
                            }
                          } else if (sf[iii + 1].initial_style == 1) {
                            image[(j * sr->xr + i) * 3] =
                                sf[iii + 1].initial_line_color[0];
                            image[(j * sr->xr + i) * 3 + 1] =
                                sf[iii + 1].initial_line_color[1];
                            image[(j * sr->xr + i) * 3 + 2] =
                                sf[iii + 1].initial_line_color[2];
                            z_buffer[j * sr->xr + i] = depth;
                          }
                        }
                      }
                    } else {
                      for (j = last_j; j >= next_j; j--) {
                        x          = xd * (i - start_x + 0.5) + scr_area[0];
                        y          = yd * (j - start_y + 0.5) + scr_area[2];
                        point_s[0] = x;
                        point_s[1] = y;
                        point_s[2] = 0.0;
                        tranverse_transform(sr->screen_point, point_s,
                                            inv_matrix, point_o);

                        intersection =
                            find_point_depth(f, point_o, sr->view_point_d, f_f,
                                             inter_point, normal_flag);
                        if (intersection == 1)
                          depth =
                              sqrt(SQR(inter_point[0] - sr->view_point_d[0]) +
                                   SQR(inter_point[1] - sr->view_point_d[1]) +
                                   SQR(inter_point[2] - sr->view_point_d[2]));
                        else
                          depth = 1.0E17;
                        if (depth <= z_buffer[j * sr->xr + i] + EPSILON) {
                          if (sf[iii + 1].initial_style == 4) {
                            count_dash_num++;
                            if ((count_dash_num % 4) < 1) {
                              image[(j * sr->xr + i) * 3] =
                                  sf[iii + 1].initial_line_color[0];
                              image[(j * sr->xr + i) * 3 + 1] =
                                  sf[iii + 1].initial_line_color[1];
                              image[(j * sr->xr + i) * 3 + 2] =
                                  sf[iii + 1].initial_line_color[2];
                              z_buffer[j * sr->xr + i] = depth;
                            }
                          } else if (sf[iii + 1].initial_style == 1) {
                            image[(j * sr->xr + i) * 3] =
                                sf[iii + 1].initial_line_color[0];
                            image[(j * sr->xr + i) * 3 + 1] =
                                sf[iii + 1].initial_line_color[1];
                            image[(j * sr->xr + i) * 3 + 2] =
                                sf[iii + 1].initial_line_color[2];
                            z_buffer[j * sr->xr + i] = depth;
                          }
                        }
                      }
                    }

                    last_j = next_j;
                  }
                } else {
                  i = pixel_i[0][0];
                  if (pixel_i[0][1] <= pixel_i[1][1]) {
                    for (j = pixel_i[0][1]; j <= pixel_i[1][1]; j++) {
                      x = xd * (i - start_x + 0.5) + scr_area[0];
                      y = yd * (j - start_y + 0.5) + scr_area[2];
                      /*transform to the original frame */
                      point_s[0] = x;
                      point_s[1] = y;
                      point_s[2] = 0.0;
                      tranverse_transform(sr->screen_point, point_s, inv_matrix,
                                          point_o);

                      intersection =
                          find_point_depth(f, point_o, sr->view_point_d, f_f,
                                           inter_point, normal_flag);
                      if (intersection == 1)
                        depth = sqrt(SQR(inter_point[0] - sr->view_point_d[0]) +
                                     SQR(inter_point[1] - sr->view_point_d[1]) +
                                     SQR(inter_point[2] - sr->view_point_d[2]));
                      else
                        depth = 1.0E17;
                      if (depth <= z_buffer[j * sr->xr + i] + EPSILON) {
                        if (sf[iii + 1].initial_style == 4) {
                          count_dash_num++;
                          if ((count_dash_num % 4) < 1) {
                            image[(j * sr->xr + i) * 3] =
                                sf[iii + 1].initial_line_color[0];
                            image[(j * sr->xr + i) * 3 + 1] =
                                sf[iii + 1].initial_line_color[1];
                            image[(j * sr->xr + i) * 3 + 2] =
                                sf[iii + 1].initial_line_color[2];
                            z_buffer[j * sr->xr + i] = depth;
                          }
                        } else if (sf[iii + 1].initial_style == 1) {
                          image[(j * sr->xr + i) * 3] =
                              sf[iii + 1].initial_line_color[0];
                          image[(j * sr->xr + i) * 3 + 1] =
                              sf[iii + 1].initial_line_color[1];
                          image[(j * sr->xr + i) * 3 + 2] =
                              sf[iii + 1].initial_line_color[2];
                          z_buffer[j * sr->xr + i] = depth;
                        }
                      }
                    }
                  } else if (pixel_i[0][1] > pixel_i[1][1]) {
                    for (j = pixel_i[1][1]; j <= pixel_i[0][1]; j++) {
                      x = xd * (i - start_x + 0.5) + scr_area[0];
                      y = yd * (j - start_y + 0.5) + scr_area[2];
                      /*transform to the original frame */
                      point_s[0] = x;
                      point_s[1] = y;
                      point_s[2] = 0.0;
                      tranverse_transform(sr->screen_point, point_s, inv_matrix,
                                          point_o);

                      intersection =
                          find_point_depth(f, point_o, sr->view_point_d, f_f,
                                           inter_point, normal_flag);
                      if (intersection == 1)
                        depth = sqrt(SQR(inter_point[0] - sr->view_point_d[0]) +
                                     SQR(inter_point[1] - sr->view_point_d[1]) +
                                     SQR(inter_point[2] - sr->view_point_d[2]));
                      else
                        depth = 1.0E17;
                      if (depth <= z_buffer[j * sr->xr + i] + EPSILON) {
                        if (sf[iii + 1].initial_style == 4) {
                          count_dash_num++;
                          if ((count_dash_num % 4) < 1) {
                            image[(j * sr->xr + i) * 3] =
                                sf[iii + 1].initial_line_color[0];
                            image[(j * sr->xr + i) * 3 + 1] =
                                sf[iii + 1].initial_line_color[1];
                            image[(j * sr->xr + i) * 3 + 2] =
                                sf[iii + 1].initial_line_color[2];
                            z_buffer[j * sr->xr + i] = depth;
                          }
                        } else if (sf[iii + 1].initial_style == 1) {
                          image[(j * sr->xr + i) * 3] =
                              sf[iii + 1].initial_line_color[0];
                          image[(j * sr->xr + i) * 3 + 1] =
                              sf[iii + 1].initial_line_color[1];
                          image[(j * sr->xr + i) * 3 + 2] =
                              sf[iii + 1].initial_line_color[2];
                          z_buffer[j * sr->xr + i] = depth;
                        }
                      }
                    }
                  }
                }
              }
            }
          }

          if ((sf[iii + 1].deform_style == 1) ||
              (sf[iii + 1].deform_style == 4)) {
            normal_flag = 1;
            for (k = 0; k < result[iii].n_patch; k++) {
              for (i = 0; i < 3; i++) {
                if (stat_para[70] == 0) {
                  for (j = 0; j < 3; j++)
                    f_disp[i][j] =
                        result[iii]
                            .vertex[(result[iii].patch[k * 3 + i] - 1) * 3 +
                                    j] +
                        result[iii]
                                .disp[(result[iii].patch[k * 3 + i] - 1) * 3 +
                                      j] *
                            (sf[iii + 1].disp_scale /
                             (sr->deform_num_of_frames) * (jjj + 1));
                } else {
                  for (j = 0; j < 3; j++)
                    f_disp[i][j] =
                        result[iii]
                            .vertex[(result[iii].patch[k * 3 + i] - 1) * 3 +
                                    j] +
                        result[iii]
                                .disp[(result[iii].patch[k * 3 + i] - 1) * 3 +
                                      j] *
                            sf[iii + 1].disp_scale;
                }
                c_value[i] =
                    result[iii].color[result[iii].patch[k * 3 + i] - 1];
              }
              for (mmm = 0; mmm < 3; mmm++) {
                mmm1                = mmm + 1;
                if (mmm1 == 3) mmm1 = 0;
                for (j = 0; j < 3; j++) {
                  iso_p[j]     = f_disp[mmm][j];
                  iso_p[j + 3] = f_disp[mmm1][j];
                }

                transform_frame4(sr->screen_point, iso_p, coff_matrix, n_iso);
                find_projection_range3(view_point, n_iso, pixel_d, iso_p);
                for (j = 0; j < 2; j++) {
                  pixel_i[j][0] =
                      (int)((pixel_d[j][0] - scr_area[0]) / xd - 0.5) + start_x;
                  pixel_i[j][1] =
                      (int)((pixel_d[j][1] - scr_area[2]) / yd + 0.5) + start_y;
                }
                if (pixel_i[1][0] != pixel_i[0][0]) {
                  slip_y = (double)(pixel_i[1][1] - pixel_i[0][1]) /
                           (double)(pixel_i[1][0] - pixel_i[0][0]);
                  last_j = pixel_i[0][1];
                  for (i = pixel_i[0][0]; i <= pixel_i[1][0]; i++) {
                    next_j = (int)(pixel_i[0][1] +
                                   slip_y * (double)(i - pixel_i[0][0]));
                    if (last_j <= next_j) {
                      for (j = last_j; j <= next_j; j++) {
                        x          = xd * (i - start_x + 0.5) + scr_area[0];
                        y          = yd * (j - start_y + 0.5) + scr_area[2];
                        point_s[0] = x;
                        point_s[1] = y;
                        point_s[2] = 0.0;
                        tranverse_transform(sr->screen_point, point_s,
                                            inv_matrix, point_o);

                        intersection =
                            find_point_depth(f_disp, point_o, sr->view_point_d,
                                             f_f, inter_point, normal_flag);
                        if (intersection == 1)
                          depth =
                              sqrt(SQR(inter_point[0] - sr->view_point_d[0]) +
                                   SQR(inter_point[1] - sr->view_point_d[1]) +
                                   SQR(inter_point[2] - sr->view_point_d[2]));
                        else
                          depth = 1.0E17;
                        if (depth <= z_buffer[j * sr->xr + i] + EPSILON) {
                          if (sf[iii + 1].deform_style == 4) {
                            count_dash_num++;
                            if ((count_dash_num % 4) < 1) {
                              image[(j * sr->xr + i) * 3] =
                                  sf[iii + 1].deform_line_color[0];
                              image[(j * sr->xr + i) * 3 + 1] =
                                  sf[iii + 1].deform_line_color[1];
                              image[(j * sr->xr + i) * 3 + 2] =
                                  sf[iii + 1].deform_line_color[2];
                              z_buffer[j * sr->xr + i] = depth;
                            }
                          } else if (sf[iii + 1].deform_style == 1) {
                            image[(j * sr->xr + i) * 3] =
                                sf[iii + 1].deform_line_color[0];
                            image[(j * sr->xr + i) * 3 + 1] =
                                sf[iii + 1].deform_line_color[1];
                            image[(j * sr->xr + i) * 3 + 2] =
                                sf[iii + 1].deform_line_color[2];
                            z_buffer[j * sr->xr + i] = depth;
                          }
                        }
                      }
                    } else {
                      for (j = last_j; j >= next_j; j--) {
                        x          = xd * (i - start_x + 0.5) + scr_area[0];
                        y          = yd * (j - start_y + 0.5) + scr_area[2];
                        point_s[0] = x;
                        point_s[1] = y;
                        point_s[2] = 0.0;
                        tranverse_transform(sr->screen_point, point_s,
                                            inv_matrix, point_o);

                        intersection =
                            find_point_depth(f_disp, point_o, sr->view_point_d,
                                             f_f, inter_point, normal_flag);
                        if (intersection == 1)
                          depth =
                              sqrt(SQR(inter_point[0] - sr->view_point_d[0]) +
                                   SQR(inter_point[1] - sr->view_point_d[1]) +
                                   SQR(inter_point[2] - sr->view_point_d[2]));
                        else
                          depth = 1.0E17;
                        if (depth <= z_buffer[j * sr->xr + i] + EPSILON) {
                          if (sf[iii + 1].deform_style == 4) {
                            count_dash_num++;
                            if ((count_dash_num % 4) < 1) {
                              image[(j * sr->xr + i) * 3] =
                                  sf[iii + 1].deform_line_color[0];
                              image[(j * sr->xr + i) * 3 + 1] =
                                  sf[iii + 1].deform_line_color[1];
                              image[(j * sr->xr + i) * 3 + 2] =
                                  sf[iii + 1].deform_line_color[2];
                              z_buffer[j * sr->xr + i] = depth;
                            }
                          } else if (sf[iii + 1].deform_style == 1) {
                            image[(j * sr->xr + i) * 3] =
                                sf[iii + 1].deform_line_color[0];
                            image[(j * sr->xr + i) * 3 + 1] =
                                sf[iii + 1].deform_line_color[1];
                            image[(j * sr->xr + i) * 3 + 2] =
                                sf[iii + 1].deform_line_color[2];
                            z_buffer[j * sr->xr + i] = depth;
                          }
                        }
                      }
                    }

                    last_j = next_j;
                  }
                } else {
                  i = pixel_i[0][0];
                  if (pixel_i[0][1] <= pixel_i[1][1]) {
                    for (j = pixel_i[0][1]; j <= pixel_i[1][1]; j++) {
                      x = xd * (i - start_x + 0.5) + scr_area[0];
                      y = yd * (j - start_y + 0.5) + scr_area[2];
                      /*transform to the original frame */
                      point_s[0] = x;
                      point_s[1] = y;
                      point_s[2] = 0.0;
                      tranverse_transform(sr->screen_point, point_s, inv_matrix,
                                          point_o);

                      intersection =
                          find_point_depth(f_disp, point_o, sr->view_point_d,
                                           f_f, inter_point, normal_flag);
                      if (intersection == 1)
                        depth = sqrt(SQR(inter_point[0] - sr->view_point_d[0]) +
                                     SQR(inter_point[1] - sr->view_point_d[1]) +
                                     SQR(inter_point[2] - sr->view_point_d[2]));
                      else
                        depth = 1.0E17;
                      if (depth <= z_buffer[j * sr->xr + i] + EPSILON) {
                        if (sf[iii + 1].deform_style == 4) {
                          count_dash_num++;
                          if ((count_dash_num % 4) < 1) {
                            image[(j * sr->xr + i) * 3] =
                                sf[iii + 1].deform_line_color[0];
                            image[(j * sr->xr + i) * 3 + 1] =
                                sf[iii + 1].deform_line_color[1];
                            image[(j * sr->xr + i) * 3 + 2] =
                                sf[iii + 1].deform_line_color[2];
                            z_buffer[j * sr->xr + i] = depth;
                          }
                        } else if (sf[iii + 1].deform_style == 1) {
                          image[(j * sr->xr + i) * 3] =
                              sf[iii + 1].deform_line_color[0];
                          image[(j * sr->xr + i) * 3 + 1] =
                              sf[iii + 1].deform_line_color[1];
                          image[(j * sr->xr + i) * 3 + 2] =
                              sf[iii + 1].deform_line_color[2];
                          z_buffer[j * sr->xr + i] = depth;
                        }
                      }
                    }
                  } else if (pixel_i[0][1] > pixel_i[1][1]) {
                    for (j = pixel_i[1][1]; j <= pixel_i[0][1]; j++) {
                      x = xd * (i - start_x + 0.5) + scr_area[0];
                      y = yd * (j - start_y + 0.5) + scr_area[2];
                      /*transform to the original frame */
                      point_s[0] = x;
                      point_s[1] = y;
                      point_s[2] = 0.0;
                      tranverse_transform(sr->screen_point, point_s, inv_matrix,
                                          point_o);

                      intersection =
                          find_point_depth(f_disp, point_o, sr->view_point_d,
                                           f_f, inter_point, normal_flag);
                      if (intersection == 1)
                        depth = sqrt(SQR(inter_point[0] - sr->view_point_d[0]) +
                                     SQR(inter_point[1] - sr->view_point_d[1]) +
                                     SQR(inter_point[2] - sr->view_point_d[2]));
                      else
                        depth = 1.0E17;
                      if (depth <= z_buffer[j * sr->xr + i] + EPSILON) {
                        if (sf[iii + 1].deform_style == 4) {
                          count_dash_num++;
                          if ((count_dash_num % 4) < 1) {
                            image[(j * sr->xr + i) * 3] =
                                sf[iii + 1].deform_line_color[0];
                            image[(j * sr->xr + i) * 3 + 1] =
                                sf[iii + 1].deform_line_color[1];
                            image[(j * sr->xr + i) * 3 + 2] =
                                sf[iii + 1].deform_line_color[2];
                            z_buffer[j * sr->xr + i] = depth;
                          }
                        } else if (sf[iii + 1].deform_style == 1) {
                          image[(j * sr->xr + i) * 3] =
                              sf[iii + 1].deform_line_color[0];
                          image[(j * sr->xr + i) * 3 + 1] =
                              sf[iii + 1].deform_line_color[1];
                          image[(j * sr->xr + i) * 3 + 2] =
                              sf[iii + 1].deform_line_color[2];
                          z_buffer[j * sr->xr + i] = depth;
                        }
                      }
                    }
                  }
                }
              }
            }
          }

        } /* if deform_display=1 */
        else if (sf[iii + 1].deform_display_on == 0) {
          if ((sf[iii + 1].display_method != 1) &&
              (sf[iii + 1].display_method != 4)) {
            deltac = (maxcolor - mincolor) / (sf[iii + 1].isoline_number + 1);
            iso_value = (double *)HECMW_calloc(sf[iii + 1].isoline_number,
                                               sizeof(double));
            iso_rgb = (double *)HECMW_calloc(sf[iii + 1].isoline_number * 3,
                                             sizeof(double));
            for (m = 0; m < sf[iii + 1].isoline_number; m++) {
              iso_value[m] = mincolor + deltac * (m + 1);
              if (sr->isoline_color[0] == -1.0) {
                value_to_rgb(iso_value[m], color_iso, mincolor, maxcolor,
                             sr->color_mapping_style, sr->interval_point,
                             sr->interval_mapping_num, sr->color_system_type);
                iso_rgb[m * 3]     = color_iso[0];
                iso_rgb[m * 3 + 1] = color_iso[1];
                iso_rgb[m * 3 + 2] = color_iso[2];
              } else {
                iso_rgb[m * 3]     = sr->isoline_color[0];
                iso_rgb[m * 3 + 1] = sr->isoline_color[1];
                iso_rgb[m * 3 + 2] = sr->isoline_color[2];
              }
            }
          }
          if ((sf[iii + 1].surface_style == 3) && (sf[iii + 1].coef[0] == 0) &&
              (sf[iii + 1].coef[1] == 0) && (sf[iii + 1].coef[2] == 0) &&
              (sf[iii + 1].coef[3] == 0) && (sf[iii + 1].coef[4] == 0) &&
              (sf[iii + 1].coef[5] == 0)) {
            f_f[0] = sf[iii + 1].coef[6];
            f_f[1] = sf[iii + 1].coef[7];
            f_f[2] = sf[iii + 1].coef[8];
            f_f[3] = sf[iii + 1].coef[9];
            norm   = sqrt(f_f[0] * f_f[0] + f_f[1] * f_f[1] + f_f[2] * f_f[2]);
            for (m      = 0; m < 4; m++) f_f[m] /= norm;
            normal_flag = 0;
          } else {
            normal_flag = 1;
            if (sr->smooth_shading == 1) {
              v_normal = (double *)HECMW_calloc(result[iii].n_vertex * 3,
                                                sizeof(double));
              p_normal = (double *)HECMW_calloc(result[iii].n_patch * 4,
                                                sizeof(double));
              v_num = (int *)HECMW_calloc(result[iii].n_vertex, sizeof(int));
              if ((v_normal == NULL) || (v_num == NULL) || (p_normal == NULL))
                HECMW_vis_memory_exit("v_normal, p_normal");
              for (k = 0; k < result[iii].n_vertex; k++) v_num[k] = 0;
              for (k = 0; k < result[iii].n_vertex * 3; k++) v_normal[k] = 0.0;
              for (k = 0; k < result[iii].n_patch; k++) {
                for (i = 0; i < 3; i++) {
                  for (j = 0; j < 3; j++)
                    f[i][j] =
                        result[iii]
                            .vertex[(result[iii].patch[k * 3 + i] - 1) * 3 + j];
                }
                p_normal[k * 4 + 0] =
                    (f[1][1] - f[0][1]) * (f[2][2] - f[0][2]) -
                    (f[2][1] - f[0][1]) * (f[1][2] - f[0][2]);
                p_normal[k * 4 + 1] =
                    -(f[1][0] - f[0][0]) * (f[2][2] - f[0][2]) +
                    (f[2][0] - f[0][0]) * (f[1][2] - f[0][2]);
                p_normal[k * 4 + 2] =
                    (f[1][0] - f[0][0]) * (f[2][1] - f[0][1]) -
                    (f[2][0] - f[0][0]) * (f[1][1] - f[0][1]);
                norm = sqrt(p_normal[k * 4 + 0] * p_normal[k * 4 + 0] +
                            p_normal[k * 4 + 1] * p_normal[k * 4 + 1] +
                            p_normal[k * 4 + 2] * p_normal[k * 4 + 2]);
                if (fabs(norm) > EPSILON) {
                  for (j = 0; j < 3; j++) p_normal[k * 4 + j] /= norm;
                }
                p_normal[k * 4 + 3] = -p_normal[k * 4 + 0] * f[0][0] -
                                      p_normal[k * 4 + 1] * f[0][1] -
                                      p_normal[k * 4 + 2] * f[0][2];
                for (i = 0; i < 3; i++) {
                  v_num[result[iii].patch[k * 3 + i] - 1]++;
                  for (j = 0; j < 3; j++)
                    v_normal[(result[iii].patch[k * 3 + i] - 1) * 3 + j] +=
                        p_normal[k * 4 + j];
                }
              }
              for (k = 0; k < result[iii].n_vertex; k++) {
                if (v_num[k] != 0) {
                  for (j = 0; j < 3; j++) v_normal[k * 3 + j] /= v_num[k];
                }
              }
            }
          }

          for (k = 0; k < result[iii].n_patch; k++) {
            for (i = 0; i < 3; i++) {
              for (j = 0; j < 3; j++)
                f[i][j] =
                    result[iii]
                        .vertex[(result[iii].patch[k * 3 + i] - 1) * 3 + j];
              c_value[i] = result[iii].color[result[iii].patch[k * 3 + i] - 1];
            }
            if (sf[iii + 1].display_method != 2) {
              if (sr->smooth_shading == 1) {
                for (i = 0; i < 3; i++)
                  for (j = 0; j < 3; j++)
                    normal[i * 3 + j] =
                        v_normal[(result[iii].patch[k * 3 + i] - 1) * 3 + j];
              }

              transform_frame3(sr->screen_point, f, coff_matrix, n_f);
              find_projection_range2(view_point, n_f, scr_f);
              starti = (int)((scr_f[0] - scr_area[0]) / xd) + start_x - 1;
              endi   = (int)((scr_f[1] - scr_area[0]) / xd) + start_x + 1;
              startj = (int)((scr_f[2] - scr_area[2]) / yd) + start_y - 1;
              endj   = (int)((scr_f[3] - scr_area[2]) / yd) + start_y + 1;
              for (j = startj; j < endj; j++)
                for (i = starti; i < endi; i++) {
                  x = xd * (i - start_x + 0.5) + scr_area[0];
                  y = yd * (j - start_y + 0.5) + scr_area[2];
                  /*transform to the original frame */
                  point_s[0] = x;
                  point_s[1] = y;
                  point_s[2] = 0.0;
                  tranverse_transform(sr->screen_point, point_s, inv_matrix,
                                      point_o);

                  intersection = find_surface_point(
                      f, point_o, sr->view_point_d, inter_point, f_f, v_f,
                      c_value, &value, normal, normal_flag, sr->smooth_shading);
                  if (intersection == 1) {
                    depth = sqrt(SQR(inter_point[0] - sr->view_point_d[0]) +
                                 SQR(inter_point[1] - sr->view_point_d[1]) +
                                 SQR(inter_point[2] - sr->view_point_d[2]));
                    if (depth < z_buffer[j * sr->xr + i]) {
                      z_buffer[j * sr->xr + i] = depth;
                      if (sf[iii + 1].display_method == 5) {
                        tmp_i =
                            (int)((value - mincolor) / (maxcolor - mincolor) *
                                  (sf[iii + 1].isoline_number + 1));
                        value = (maxcolor - mincolor) /
                                    (sf[iii + 1].isoline_number + 1) *
                                    (double)tmp_i +
                                mincolor;
                      }

                      compute_color_sf(
                          inter_point, value, normal, sr->color_mapping_style,
                          sr->interval_point, sr->view_point_d,
                          sr->interval_mapping_num, sr->color_system_type,
                          sr->num_of_lights, sr->light_point, sr->k_ads,
                          accum_rgba, mincolor, maxcolor,
                          sf[iii + 1].display_method);
                      image[(j * sr->xr + i) * 3]     = accum_rgba[0];
                      image[(j * sr->xr + i) * 3 + 1] = accum_rgba[1];
                      image[(j * sr->xr + i) * 3 + 2] = accum_rgba[2];
                    }
                  }
                }
            }

            if ((sf[iii + 1].display_method == 2) ||
                (sf[iii + 1].display_method == 3)) {
              for (m = 0; m < sf[iii + 1].isoline_number; m++) {
                isocolor = iso_value[m];
                flag     = 0;
                flag     = find_line_segment(f, c_value, isocolor, iso_p);
                if (flag == 1) {
                  transform_frame4(sr->screen_point, iso_p, coff_matrix, n_iso);
                  find_projection_range3(view_point, n_iso, pixel_d, iso_p);
                  for (j = 0; j < 2; j++) {
                    pixel_i[j][0] =
                        (int)((pixel_d[j][0] - scr_area[0]) / xd - 0.5) +
                        start_x;
                    pixel_i[j][1] =
                        (int)((pixel_d[j][1] - scr_area[2]) / yd + 0.5) +
                        start_y;
                  }
                  if (pixel_i[1][0] != pixel_i[0][0]) {
                    slip_y = (double)(pixel_i[1][1] - pixel_i[0][1]) /
                             (double)(pixel_i[1][0] - pixel_i[0][0]);
                    last_j = pixel_i[0][1];
                    for (i = pixel_i[0][0]; i <= pixel_i[1][0]; i++) {
                      next_j = (int)(pixel_i[0][1] +
                                     slip_y * (double)(i - pixel_i[0][0]));
                      if (last_j <= next_j) {
                        for (j = last_j; j <= next_j; j++) {
                          x          = xd * (i - start_x + 0.5) + scr_area[0];
                          y          = yd * (j - start_y + 0.5) + scr_area[2];
                          point_s[0] = x;
                          point_s[1] = y;
                          point_s[2] = 0.0;
                          tranverse_transform(sr->screen_point, point_s,
                                              inv_matrix, point_o);

                          intersection =
                              find_point_depth(f, point_o, sr->view_point_d,
                                               f_f, inter_point, normal_flag);
                          if (intersection == 1)
                            depth =
                                sqrt(SQR(inter_point[0] - sr->view_point_d[0]) +
                                     SQR(inter_point[1] - sr->view_point_d[1]) +
                                     SQR(inter_point[2] - sr->view_point_d[2]));
                          else
                            depth = 1.0E17;
                          if (depth <= z_buffer[j * sr->xr + i] + EPSILON) {
                            image[(j * sr->xr + i) * 3] = iso_rgb[m * 3];
                            image[(j * sr->xr + i) * 3 + 1] =
                                iso_rgb[m * 3 + 1];
                            image[(j * sr->xr + i) * 3 + 2] =
                                iso_rgb[m * 3 + 2];
                            z_buffer[j * sr->xr + i] = depth;
                          }
                        }
                      } else {
                        for (j = last_j; j >= next_j; j--) {
                          x          = xd * (i - start_x + 0.5) + scr_area[0];
                          y          = yd * (j - start_y + 0.5) + scr_area[2];
                          point_s[0] = x;
                          point_s[1] = y;
                          point_s[2] = 0.0;
                          tranverse_transform(sr->screen_point, point_s,
                                              inv_matrix, point_o);

                          intersection =
                              find_point_depth(f, point_o, sr->view_point_d,
                                               f_f, inter_point, normal_flag);
                          if (intersection == 1)
                            depth =
                                sqrt(SQR(inter_point[0] - sr->view_point_d[0]) +
                                     SQR(inter_point[1] - sr->view_point_d[1]) +
                                     SQR(inter_point[2] - sr->view_point_d[2]));
                          else
                            depth = 1.0E17;
                          if (depth <= z_buffer[j * sr->xr + i] + EPSILON) {
                            image[(j * sr->xr + i) * 3] = iso_rgb[m * 3];
                            image[(j * sr->xr + i) * 3 + 1] =
                                iso_rgb[m * 3 + 1];
                            image[(j * sr->xr + i) * 3 + 2] =
                                iso_rgb[m * 3 + 2];
                            z_buffer[j * sr->xr + i] = depth;
                          }
                        }
                      }

                      last_j = next_j;
                    }
                  } else {
                    i = pixel_i[0][0];
                    if (pixel_i[0][1] <= pixel_i[1][1]) {
                      for (j = pixel_i[0][1]; j <= pixel_i[1][1]; j++) {
                        x = xd * (i - start_x + 0.5) + scr_area[0];
                        y = yd * (j - start_y + 0.5) + scr_area[2];
                        /*transform to the original frame */
                        point_s[0] = x;
                        point_s[1] = y;
                        point_s[2] = 0.0;
                        tranverse_transform(sr->screen_point, point_s,
                                            inv_matrix, point_o);

                        intersection =
                            find_point_depth(f, point_o, sr->view_point_d, f_f,
                                             inter_point, normal_flag);
                        if (intersection == 1)
                          depth =
                              sqrt(SQR(inter_point[0] - sr->view_point_d[0]) +
                                   SQR(inter_point[1] - sr->view_point_d[1]) +
                                   SQR(inter_point[2] - sr->view_point_d[2]));
                        else
                          depth = 1.0E17;
                        if (depth <= z_buffer[j * sr->xr + i] + EPSILON) {
                          image[(j * sr->xr + i) * 3]     = iso_rgb[m * 3];
                          image[(j * sr->xr + i) * 3 + 1] = iso_rgb[m * 3 + 1];
                          image[(j * sr->xr + i) * 3 + 2] = iso_rgb[m * 3 + 2];
                          z_buffer[j * sr->xr + i]        = depth;
                        }
                      }
                    } else if (pixel_i[0][1] > pixel_i[1][1]) {
                      for (j = pixel_i[1][1]; j <= pixel_i[0][1]; j++) {
                        x = xd * (i - start_x + 0.5) + scr_area[0];
                        y = yd * (j - start_y + 0.5) + scr_area[2];
                        /*transform to the original frame */
                        point_s[0] = x;
                        point_s[1] = y;
                        point_s[2] = 0.0;
                        tranverse_transform(sr->screen_point, point_s,
                                            inv_matrix, point_o);

                        intersection =
                            find_point_depth(f, point_o, sr->view_point_d, f_f,
                                             inter_point, normal_flag);
                        if (intersection == 1)
                          depth =
                              sqrt(SQR(inter_point[0] - sr->view_point_d[0]) +
                                   SQR(inter_point[1] - sr->view_point_d[1]) +
                                   SQR(inter_point[2] - sr->view_point_d[2]));
                        else
                          depth = 1.0E17;
                        if (depth <= z_buffer[j * sr->xr + i] + EPSILON) {
                          image[(j * sr->xr + i) * 3]     = iso_rgb[m * 3];
                          image[(j * sr->xr + i) * 3 + 1] = iso_rgb[m * 3 + 1];
                          image[(j * sr->xr + i) * 3 + 2] = iso_rgb[m * 3 + 2];
                          z_buffer[j * sr->xr + i]        = depth;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          if ((sf[iii + 1].display_method != 1) &&
              (sf[iii + 1].display_method != 4)) {
            HECMW_free(iso_value);
            HECMW_free(iso_rgb);
          }
        } /* if deform_display=0 */
      }

      HECMW_Barrier(VIS_COMM);
      fprintf(stderr, "Finish the computing on each PE\n");
      if (mynode == 0) {
        /*
fprintf(stderr, "Finish the computing on each PE\n");
t3=HECMW_Wtime();
fprintf(stderr, "The  volume rendering on each PE is %lf the total time up to
now is %lf\n", t3-t2, t3-t1);
t2=t3;
         */
      }
      if (pesize > 1) {
        pix_num               = (int)(sr->xr * sr->yr / pesize);
        pix0_num              = sr->xr * sr->yr - pix_num * (pesize - 1);
        if (mynode == 0) pixn = pix0_num;
        if (mynode != 0) pixn = pix_num;
        n_subimage = (double *)HECMW_calloc(pesize * pixn * 3, sizeof(double));
        n_subopa   = (double *)HECMW_calloc(pesize * pixn, sizeof(double));
        if ((n_subimage == NULL) || (n_subopa == NULL))
          HECMW_vis_memory_exit("n_subimage, n_subopa");

        for (j = 0; j < pixn * 3; j++) {
          if (mynode == 0)
            n_subimage[j] = image[j];
          else if (mynode != 0)
            n_subimage[mynode * pixn * 3 + j] =
                image[pix0_num * 3 + (mynode - 1) * pix_num * 3 + j];
        }
        if (mynode == 0)
          for (j = 0; j < pixn; j++) n_subopa[j] = z_buffer[j];
        else if (mynode != 0)
          for (j = 0; j < pixn; j++)
            n_subopa[mynode * pixn + j] =
                z_buffer[pix0_num + (mynode - 1) * pix_num + j];

        subimage = (double *)HECMW_calloc(pixn * 3, sizeof(double));
        subopa   = (double *)HECMW_calloc(pixn, sizeof(double));
        if ((subimage == NULL) || (subopa == NULL))
          HECMW_vis_memory_exit("subimage, subopa");
        HECMW_Barrier(VIS_COMM);

        for (i = 0; i < pesize; i++) {
          if (mynode == i) {
            for (k = 0; k < pesize; k++) {
              if (k != i) {
                if (k == 0) {
                  subimage1 =
                      (double *)HECMW_calloc(pix0_num * 3, sizeof(double));
                  subopa1 = (double *)HECMW_calloc(pix0_num, sizeof(double));
                } else if (k != 0) {
                  subimage1 =
                      (double *)HECMW_calloc(pix_num * 3, sizeof(double));
                  subopa1 = (double *)HECMW_calloc(pix_num, sizeof(double));
                }

                if ((subimage1 == NULL) || (subopa1 == NULL))
                  HECMW_vis_memory_exit("subimage1, subopa1");
                if (k == 0) {
                  for (j = 0; j < pix0_num * 3; j++) subimage1[j] = image[j];
                  for (j = 0; j < pix0_num; j++) subopa1[j] = z_buffer[j];
                } else if (k != 0) {
                  for (j = 0; j < pix_num * 3; j++)
                    subimage1[j] =
                        image[pix0_num * 3 + (k - 1) * pix_num * 3 + j];
                  for (j       = 0; j < pix_num; j++)
                    subopa1[j] = z_buffer[pix0_num + (k - 1) * pix_num + j];
                }
                if (k == 0) {
                  HECMW_Send(subimage1, pix0_num * 3, HECMW_DOUBLE, k, 0,
                             VIS_COMM);
                  HECMW_Send(subopa1, pix0_num, HECMW_DOUBLE, k, 0, VIS_COMM);
                } else if (k != 0) {
                  HECMW_Send(subimage1, pix_num * 3, HECMW_DOUBLE, k, 0,
                             VIS_COMM);
                  HECMW_Send(subopa1, pix_num, HECMW_DOUBLE, k, 0, VIS_COMM);
                }

                HECMW_free(subimage1);
                HECMW_free(subopa1);
              } /*if k!=i */
            }   /*loop k*/
          }
          if (mynode != i) {
            HECMW_Recv(subimage, pixn * 3, HECMW_DOUBLE, i, HECMW_ANY_TAG,
                       VIS_COMM, &stat);
            HECMW_Recv(subopa, pixn, HECMW_DOUBLE, i, HECMW_ANY_TAG, VIS_COMM,
                       &stat);
            for (j                         = 0; j < pixn * 3; j++)
              n_subimage[i * pixn * 3 + j] = subimage[j];
            for (j = 0; j < pixn; j++) n_subopa[i * pixn + j] = subopa[j];
          }
          HECMW_Barrier(VIS_COMM);
        }

        /*	if(time_step==0) {
         */

        composite_subimage_sf(pesize, pixn, n_subimage, n_subopa, subimage,
                              subopa);
        HECMW_free(n_subimage);
        HECMW_free(n_subopa);
        /*	HECMW_free(subopa);
         */
        /*send subimage to master PE */
        if (mynode != 0) {
          HECMW_Send(subimage, pixn * 3, HECMW_DOUBLE, 0, 0, VIS_COMM);
          HECMW_Send(subopa, pixn, HECMW_DOUBLE, 0, 0, VIS_COMM);
        }

        if (mynode == 0) {
          for (j = 0; j < pix0_num * 3; j++) image[j] = subimage[j];
          for (j = 0; j < pix0_num; j++) z_buffer[j] = subopa[j];
          for (i = 1; i < pesize; i++) {
            HECMW_Recv(subimage, pix_num * 3, HECMW_DOUBLE, i, HECMW_ANY_TAG,
                       VIS_COMM, &stat);
            HECMW_Recv(subopa, pix_num, HECMW_DOUBLE, i, HECMW_ANY_TAG,
                       VIS_COMM, &stat);

            for (j = 0; j < pix_num * 3; j++)
              image[pix0_num * 3 + (i - 1) * pix_num * 3 + j] = subimage[j];
            for (j                                       = 0; j < pix_num; j++)
              z_buffer[pix0_num + (i - 1) * pix_num + j] = subopa[j];
          }
          if (sr->boundary_line_on == 1) { /* draw boundary line */
            fprintf(stderr,
                    "Drawing boundary line: trange=%lf %lf %lf %lf %lf %lf\n",
                    trange[0], trange[1], trange[2], trange[3], trange[4],
                    trange[5]);
            for (m1 = 0; m1 < 6; m1++)
              for (m2 = 0; m2 < 4; m2++) {
                if (((m1 == 2) && (m2 == 0)) || ((m1 == 4) && (m2 == 1))) {
                  iso_p[0] = trange[0];
                  iso_p[1] = trange[2];
                  iso_p[2] = trange[4];
                  iso_p[3] = trange[1];
                  iso_p[4] = trange[2];
                  iso_p[5] = trange[4];
                } else if (((m1 == 0) && (m2 == 0)) ||
                           ((m1 == 4) && (m2 == 0))) {
                  iso_p[0] = trange[0];
                  iso_p[1] = trange[2];
                  iso_p[2] = trange[4];
                  iso_p[3] = trange[0];
                  iso_p[4] = trange[3];
                  iso_p[5] = trange[4];
                } else if (((m1 == 0) && (m2 == 1)) ||
                           ((m1 == 2) && (m2 == 3))) {
                  iso_p[0] = trange[0];
                  iso_p[1] = trange[2];
                  iso_p[2] = trange[4];
                  iso_p[3] = trange[0];
                  iso_p[4] = trange[2];
                  iso_p[5] = trange[5];
                } else if (((m1 == 1) && (m2 == 0)) ||
                           ((m1 == 4) && (m2 == 2))) {
                  iso_p[0] = trange[1];
                  iso_p[1] = trange[2];
                  iso_p[2] = trange[4];
                  iso_p[3] = trange[1];
                  iso_p[4] = trange[3];
                  iso_p[5] = trange[4];
                } else if (((m1 == 1) && (m2 == 1)) ||
                           ((m1 == 2) && (m2 == 1))) {
                  iso_p[0] = trange[1];
                  iso_p[1] = trange[2];
                  iso_p[2] = trange[4];
                  iso_p[3] = trange[1];
                  iso_p[4] = trange[2];
                  iso_p[5] = trange[5];
                } else if (((m1 == 2) && (m2 == 2)) ||
                           ((m1 == 5) && (m2 == 1))) {
                  iso_p[0] = trange[1];
                  iso_p[1] = trange[2];
                  iso_p[2] = trange[5];
                  iso_p[3] = trange[0];
                  iso_p[4] = trange[2];
                  iso_p[5] = trange[5];
                } else if (((m1 == 1) && (m2 == 2)) ||
                           ((m1 == 5) && (m2 == 2))) {
                  iso_p[0] = trange[1];
                  iso_p[1] = trange[2];
                  iso_p[2] = trange[5];
                  iso_p[3] = trange[1];
                  iso_p[4] = trange[3];
                  iso_p[5] = trange[5];
                } else if (((m1 == 0) && (m2 == 2)) ||
                           ((m1 == 5) && (m2 == 0))) {
                  iso_p[0] = trange[0];
                  iso_p[1] = trange[2];
                  iso_p[2] = trange[5];
                  iso_p[3] = trange[0];
                  iso_p[4] = trange[3];
                  iso_p[5] = trange[5];
                } else if (((m1 == 0) && (m2 == 3)) ||
                           ((m1 == 3) && (m2 == 0))) {
                  iso_p[0] = trange[0];
                  iso_p[1] = trange[3];
                  iso_p[2] = trange[4];
                  iso_p[3] = trange[0];
                  iso_p[4] = trange[3];
                  iso_p[5] = trange[5];
                } else if (((m1 == 3) && (m2 == 1)) ||
                           ((m1 == 4) && (m2 == 3))) {
                  iso_p[0] = trange[0];
                  iso_p[1] = trange[3];
                  iso_p[2] = trange[4];
                  iso_p[3] = trange[1];
                  iso_p[4] = trange[3];
                  iso_p[5] = trange[4];
                } else if (((m1 == 1) && (m2 == 3)) ||
                           ((m1 == 3) && (m2 == 2))) {
                  iso_p[0] = trange[1];
                  iso_p[1] = trange[3];
                  iso_p[2] = trange[4];
                  iso_p[3] = trange[1];
                  iso_p[4] = trange[3];
                  iso_p[5] = trange[5];
                } else if (((m1 == 3) && (m2 == 3)) ||
                           ((m1 == 5) && (m2 == 3))) {
                  iso_p[0] = trange[1];
                  iso_p[1] = trange[3];
                  iso_p[2] = trange[5];
                  iso_p[3] = trange[0];
                  iso_p[4] = trange[3];
                  iso_p[5] = trange[5];
                }

                if (m1 == 0) {
                  f[0][0] = trange[0];
                  f[0][1] = trange[2];
                  f[0][2] = trange[4];
                  f[1][0] = trange[0];
                  f[1][1] = trange[3];
                  f[1][2] = trange[4];
                  f[2][0] = trange[0];
                  f[2][1] = trange[3];
                  f[2][2] = trange[5];
                } else if (m1 == 1) {
                  f[0][0] = trange[1];
                  f[0][1] = trange[2];
                  f[0][2] = trange[4];
                  f[1][0] = trange[1];
                  f[1][1] = trange[2];
                  f[1][2] = trange[5];
                  f[2][0] = trange[1];
                  f[2][1] = trange[3];
                  f[2][2] = trange[4];
                } else if (m1 == 2) {
                  f[0][0] = trange[0];
                  f[0][1] = trange[2];
                  f[0][2] = trange[4];
                  f[1][0] = trange[0];
                  f[1][1] = trange[2];
                  f[1][2] = trange[5];
                  f[2][0] = trange[1];
                  f[2][1] = trange[2];
                  f[2][2] = trange[5];
                } else if (m1 == 3) {
                  f[0][0] = trange[0];
                  f[0][1] = trange[3];
                  f[0][2] = trange[4];
                  f[1][0] = trange[1];
                  f[1][1] = trange[3];
                  f[1][2] = trange[4];
                  f[2][0] = trange[1];
                  f[2][1] = trange[3];
                  f[2][2] = trange[5];
                }

                else if (m1 == 4) {
                  f[0][0] = trange[0];
                  f[0][1] = trange[2];
                  f[0][2] = trange[4];
                  f[1][0] = trange[1];
                  f[1][1] = trange[2];
                  f[1][2] = trange[4];
                  f[2][0] = trange[0];
                  f[2][1] = trange[3];
                  f[2][2] = trange[4];
                } else if (m1 == 5) {
                  f[0][0] = trange[0];
                  f[0][1] = trange[3];
                  f[0][2] = trange[5];
                  f[1][0] = trange[1];
                  f[1][1] = trange[3];
                  f[1][2] = trange[5];
                  f[2][0] = trange[1];
                  f[2][1] = trange[2];
                  f[2][2] = trange[5];
                }

                transform_frame4(sr->screen_point, iso_p, coff_matrix, n_iso);
                find_projection_range3(view_point, n_iso, pixel_d, iso_p);
                for (j = 0; j < 2; j++) {
                  pixel_i[j][0] =
                      (int)((pixel_d[j][0] - scr_area[0]) / xd - 0.5) + start_x;
                  pixel_i[j][1] =
                      (int)((pixel_d[j][1] - scr_area[2]) / yd + 0.5) + start_y;
                }
                if (pixel_i[1][0] != pixel_i[0][0]) {
                  slip_y = (double)(pixel_i[1][1] - pixel_i[0][1]) /
                           (double)(pixel_i[1][0] - pixel_i[0][0]);
                  last_j = pixel_i[0][1];
                  for (i = pixel_i[0][0]; i <= pixel_i[1][0]; i++) {
                    next_j = (int)(pixel_i[0][1] +
                                   slip_y * (double)(i - pixel_i[0][0]));
                    if (last_j <= next_j) {
                      for (j = last_j; j <= next_j; j++) {
                        x          = xd * (i - start_x + 0.5) + scr_area[0];
                        y          = yd * (j - start_y + 0.5) + scr_area[2];
                        point_s[0] = x;
                        point_s[1] = y;
                        point_s[2] = 0.0;
                        tranverse_transform(sr->screen_point, point_s,
                                            inv_matrix, point_o);
                        intersection = find2_point_depth(
                            f, point_o, sr->view_point_d, inter_point);
                        if (intersection == 1)
                          depth =
                              sqrt(SQR(inter_point[0] - sr->view_point_d[0]) +
                                   SQR(inter_point[1] - sr->view_point_d[1]) +
                                   SQR(inter_point[2] - sr->view_point_d[2]));
                        else
                          depth = 1.0E17 + 1;
                        if (depth <= z_buffer[j * sr->xr + i] + EPSILON) {
                          image[(j * sr->xr + i) * 3]     = sr->font_color[0];
                          image[(j * sr->xr + i) * 3 + 1] = sr->font_color[1];
                          image[(j * sr->xr + i) * 3 + 2] = sr->font_color[2];
                          z_buffer[j * sr->xr + i]        = depth;
                        }
                      }
                    } else {
                      for (j = last_j; j >= next_j; j--) {
                        x          = xd * (i - start_x + 0.5) + scr_area[0];
                        y          = yd * (j - start_y + 0.5) + scr_area[2];
                        point_s[0] = x;
                        point_s[1] = y;
                        point_s[2] = 0.0;
                        tranverse_transform(sr->screen_point, point_s,
                                            inv_matrix, point_o);
                        intersection = find2_point_depth(
                            f, point_o, sr->view_point_d, inter_point);
                        if (intersection == 1)
                          depth =
                              sqrt(SQR(inter_point[0] - sr->view_point_d[0]) +
                                   SQR(inter_point[1] - sr->view_point_d[1]) +
                                   SQR(inter_point[2] - sr->view_point_d[2]));
                        else
                          depth = 1.0E17 + 1;
                        if (depth <= z_buffer[j * sr->xr + i] + EPSILON) {
                          image[(j * sr->xr + i) * 3]     = sr->font_color[0];
                          image[(j * sr->xr + i) * 3 + 1] = sr->font_color[1];
                          image[(j * sr->xr + i) * 3 + 2] = sr->font_color[2];
                          z_buffer[j * sr->xr + i]        = depth;
                        }
                      }
                    }

                    last_j = next_j;
                  }
                } else {
                  i = pixel_i[0][0];
                  if (pixel_i[0][1] <= pixel_i[1][1]) {
                    for (j = pixel_i[0][1]; j <= pixel_i[1][1]; j++) {
                      x = xd * (i - start_x + 0.5) + scr_area[0];
                      y = yd * (j - start_y + 0.5) + scr_area[2];
                      /*transform to the original frame */
                      point_s[0] = x;
                      point_s[1] = y;
                      point_s[2] = 0.0;
                      tranverse_transform(sr->screen_point, point_s, inv_matrix,
                                          point_o);
                      intersection = find2_point_depth(
                          f, point_o, sr->view_point_d, inter_point);
                      if (intersection == 1)
                        depth = sqrt(SQR(inter_point[0] - sr->view_point_d[0]) +
                                     SQR(inter_point[1] - sr->view_point_d[1]) +
                                     SQR(inter_point[2] - sr->view_point_d[2]));
                      else
                        depth = 1.0E17 + 1;
                      if (depth <= z_buffer[j * sr->xr + i] + EPSILON) {
                        image[(j * sr->xr + i) * 3]     = sr->font_color[0];
                        image[(j * sr->xr + i) * 3 + 1] = sr->font_color[1];
                        image[(j * sr->xr + i) * 3 + 2] = sr->font_color[2];
                        z_buffer[j * sr->xr + i]        = depth;
                      }
                    }
                  } else if (pixel_i[0][1] > pixel_i[1][1]) {
                    for (j = pixel_i[1][1]; j <= pixel_i[0][1]; j++) {
                      x = xd * (i - start_x + 0.5) + scr_area[0];
                      y = yd * (j - start_y + 0.5) + scr_area[2];
                      /*transform to the original frame */
                      point_s[0] = x;
                      point_s[1] = y;
                      point_s[2] = 0.0;
                      tranverse_transform(sr->screen_point, point_s, inv_matrix,
                                          point_o);
                      intersection = find2_point_depth(
                          f, point_o, sr->view_point_d, inter_point);
                      if (intersection == 1)
                        depth = sqrt(SQR(inter_point[0] - sr->view_point_d[0]) +
                                     SQR(inter_point[1] - sr->view_point_d[1]) +
                                     SQR(inter_point[2] - sr->view_point_d[2]));
                      else
                        depth = 1.0E17 + 1;
                      if (depth <= z_buffer[j * sr->xr + i] + EPSILON) {
                        image[(j * sr->xr + i) * 3]     = sr->font_color[0];
                        image[(j * sr->xr + i) * 3 + 1] = sr->font_color[1];
                        image[(j * sr->xr + i) * 3 + 2] = sr->font_color[2];
                        z_buffer[j * sr->xr + i]        = depth;
                      }
                    }
                  }
                }
              }
          }
        }
        HECMW_free(subimage);
        HECMW_free(subopa);
      }
      if (mynode == 0) {
        fprintf(stderr, " Finish combining \n");
        /*
t3=HECMW_Wtime();
fprintf(stderr, "The combining time is %lf the total time up to now is %lf\n",
t3-t2, t3-t1);
t2=t3;
         */

        /*output result in image of master pe */
        /*	  sprintf(rname, "%s-%d", outfile,
step_start+time_step*step_interv);

if ((outfp = fopen(rname, "w")) == NULL) {
fprintf(stderr, "There is no such an input data file:\n");
exit (0);
}
fprintf(outfp, "%d %d\n", vr->xr, vr->yr);
for(j=0;j<yr;j++)
for(i=xr-1;i>=0;i--)
  fprintf(outfp, "%lf %lf %lf\n", image[(j*xr+i)*3], image[(j*xr+i)*3+1],
  image[(j*xr+i)*3+2]);
fclose(outfp);
         */
        if (sr->fixed_scale_mark == 0)
          generate_color_bar(
              sr->scale_marking_on, sr->font_size, sr->color_bar_style,
              sr->mark_0_on, sr->color_mapping_bar_on, sr->xr, sr->yr,
              sr->font_color, sr->color_system_type, sr->color_mapping_style,
              sr->interval_point, sr->interval_mapping_num, sr->num_of_scale,
              tmincolor, tmaxcolor, org_mincolor, org_maxcolor, image);
        if (sr->fixed_scale_mark == 1)
          generate_color_bar(
              sr->scale_marking_on, sr->font_size, sr->color_bar_style,
              sr->mark_0_on, sr->color_mapping_bar_on, sr->xr, sr->yr,
              sr->font_color, sr->color_system_type, sr->color_mapping_style,
              sr->interval_point, sr->interval_mapping_num, sr->num_of_scale,
              tmincolor, tmaxcolor, tmincolor, tmaxcolor, image);
        /*        fprintf(stderr, "time_mark_on is %d timestep is %d\n",
         * sr->time_mark_on, timestep);
         */
        if (sr->time_mark_on == 1)
          mark_time_label(sr->font_size, sr->xr, sr->yr, sr->font_color,
                          sr->background_color, start_time, time_interval,
                          timestep, len_step, image);
        /*	(void)GEOFEM_file_translate_filename(geofem_output_file,
                 fname, GEOFEM_INDEPENDENT);
{
int n = strlen(fname);
char *ptr = fname + n;
while(*ptr != '.') {
ptr--;
}
         *ptr = '\0';
}

{
char buf[128];
strcat(fname, ".");
sprintf(buf, "%d.%d", v->count, ii);
strcat(fname, buf);
strcat(fname, ".bmp");
}

         */
        if ((sr->rotate_style > 0) && (sr->deform_num_of_frames < 2)) {
          if (ii >= 10)
            sprintf(fname, "%s.%d.bmp", outfile, ii);
          else
            sprintf(fname, "%s.0%d.bmp", outfile, ii);
        }
        if ((sr->rotate_style > 0) && (sr->deform_num_of_frames > 1)) {
          if (jjj >= 10)
            sprintf(fname, "%s.%d.%d.bmp", outfile, ii, jjj);
          else
            sprintf(fname, "%s.%d.0%d.bmp", outfile, ii, jjj);
        }
        if ((sr->rotate_style == 0) && (sr->deform_num_of_frames > 1)) {
          if (jjj >= 10)
            sprintf(fname, "%s.%d.bmp", outfile, jjj);
          else
            sprintf(fname, "%s.0%d.bmp", outfile, jjj);
        }

        if ((sr->rotate_style == 0) && (sr->deform_num_of_frames < 2))
          sprintf(fname, "%s.bmp", outfile);

        outfp = fopen(fname, "wb");

        if (!outfp)
          HECMW_vis_print_exit(
              "ERROR: HEC-MW-VIS-E0009: Cannot open output file");

        header.bfSize = 54 + 3 * sr->xr * sr->yr;
#ifdef CONVERSE_ORDER
        header.bfSize = change_unsigned_int_order(54 + 3 * sr->xr * sr->yr);
#endif
        header.bfReserved1 = 0;
#ifdef CONVERSE_ORDER
        header.bfReserved1 = change_short_int_order(0);
#endif

        header.bfReserved2 = 0;
#ifdef CONVERSE_ORDER
        header.bfReserved2 = change_short_int_order(0);
#endif

        header.bfOffBits = 54;
#ifdef CONVERSE_ORDER
        header.bfOffBits = change_unsigned_int_order(54);
#endif

        info.biBitCount = 24;
#ifdef CONVERSE_ORDER
        info.biBitCount = change_short_int_order(24);
#endif

        info.biSize = 40;
#ifdef CONVERSE_ORDER
        info.biSize = change_unsigned_int_order(40);
#endif

        info.biWidth = sr->xr;
#ifdef CONVERSE_ORDER
        info.biWidth = change_int_order(sr->xr);
#endif

        info.biHeight = sr->yr;
#ifdef CONVERSE_ORDER
        info.biHeight = change_int_order(sr->yr);
#endif

        info.biSizeImage = 3 * sr->xr * sr->yr;
#ifdef CONVERSE_ORDER
        info.biSizeImage = change_unsigned_int_order(3 * sr->xr * sr->yr);
#endif

        info.biClrImportant = 0;
#ifdef CONVERSE_ORDER
        info.biClrImportant = change_unsigned_int_order(0);
#endif

        info.biClrUsed = 0;
#ifdef CONVERSE_ORDER
        info.biClrUsed = change_unsigned_int_order(0);
#endif

        info.biCompression = 0;
#ifdef CONVERSE_ORDER
        info.biCompression = change_unsigned_int_order(0);
#endif

        info.biPlanes = 1;
#ifdef CONVERSE_ORDER
        info.biPlanes = change_short_int_order(1);
#endif

        info.biXPelsPerMeter = 3780;
#ifdef CONVERSE_ORDER
        info.biXPelsPerMeter = change_int_order(3780);
#endif

        info.biYPelsPerMeter = 3780;
#ifdef CONVERSE_ORDER
        info.biYPelsPerMeter = change_int_order(3780);
#endif

        putc('B', outfp);
        putc('M', outfp);
        fwrite(&(header.bfSize), sizeof(unsigned int), 1, outfp);
        fwrite(&header.bfReserved1, sizeof(unsigned short int), 1, outfp);
        fwrite(&header.bfReserved2, sizeof(unsigned short int), 1, outfp);
        fwrite(&header.bfOffBits, sizeof(unsigned int), 1, outfp);
        fwrite(&info, 40, 1, outfp);
        for (j = 0; j < sr->yr; j++)
          for (i = sr->xr - 1; i >= 0; i--) {
            if ((image[(j * sr->xr + i) * 3] < 0.1) &&
                (image[(j * sr->xr + i) * 3 + 1] < 0.1) &&
                (image[(j * sr->xr + i) * 3 + 2] < 0.1)) {
              ri = (int)(sr->background_color[0] * 256);
              gi = (int)(sr->background_color[1] * 256);
              bi = (int)(sr->background_color[2] * 256);
            } else {
              ri = (int)(image[(j * sr->xr + i) * 3] * 256);
              gi = (int)(image[(j * sr->xr + i) * 3 + 1] * 256);
              bi = (int)(image[(j * sr->xr + i) * 3 + 2] * 256);
            }
            if (ri < 0) ri   = 0;
            if (ri > 255) ri = 255;
            if (gi < 0) gi   = 0;
            if (gi > 255) gi = 255;
            if (bi < 0) bi   = 0;
            if (bi > 255) bi = 255;
            r                = ri;
            g                = gi;
            b                = bi;
            putc(b, outfp);
            putc(g, outfp);
            putc(r, outfp);
          }

        fclose(outfp);
      }
      HECMW_free(image);
      HECMW_free(z_buffer);
      if (sr->smooth_shading == 1) {
        HECMW_free(v_normal);
        HECMW_free(p_normal);
        HECMW_free(v_num);
      }
    }
  } /*for jjj loop */

  /*
if(pvr->surface_on==1) {
  HECMW_free(surface->color);
  HECMW_free(surface->vertex);
  HECMW_free(surface->patch);
  HECMW_free(surface);
}
   */
  /*  fclose(fp2);
   */
  /*  if(mynode==0) {
fprintf(stderr, " Finish rendering\n");
t3=HECMW_Wtime();
fprintf(stderr, " the time for rendering and output is %lf the total time of the
module is %lf\n", t3-t2, t3-t1);
}
   */

  return;
}
