/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include "hecmw_vis_combine.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hecmw_vis_comm_util.h"
#include "hecmw_vis_mem_util.h"
#include "hecmw_vis_generate_histogram_sf.h"
#include "hecmw_malloc.h"

static void find_minmax_disp(struct surface_module *sf, Result *result,
                             HECMW_Comm VIS_COMM, double disp_min[5],
                             double disp_max[5], int pesize) {
  int i, j, ii;
  double l_disp_min[5], l_disp_max[5], tmp;

  for (i = 0; i < 5; i++) {
    l_disp_min[i] = 1.0E17;
    l_disp_max[i] = -1.0E17;
  }
  for (ii = 0; ii < sf[0].surface_style; ii++) {
    for (i = 0; i < result[ii].n_vertex; i++) {
      for (j = 0; j < 3; j++) {
        if (l_disp_min[j + 1] > result[ii].disp[i * 3 + j])
          l_disp_min[j + 1] = result[ii].disp[i * 3 + j];
        if (l_disp_max[j + 1] < result[ii].disp[i * 3 + j])
          l_disp_max[j + 1] = result[ii].disp[i * 3 + j];
      }
      tmp = sqrt(result[ii].disp[i * 3] * result[ii].disp[i * 3] +
                 result[ii].disp[i * 3 + 1] * result[ii].disp[i * 3 + 1] +
                 result[ii].disp[i * 3 + 2] * result[ii].disp[i * 3 + 2]);
      if (l_disp_min[0] > tmp) l_disp_min[0] = tmp;
      if (l_disp_max[0] < tmp) l_disp_max[0] = tmp;
      if (l_disp_min[4] > result[ii].color[i])
        l_disp_min[4] = result[ii].color[i];
      if (l_disp_max[4] < result[ii].color[i])
        l_disp_max[4] = result[ii].color[i];
    }
  }
  if (pesize > 1) {
    HECMW_Allreduce(l_disp_min, disp_min, 5, HECMW_DOUBLE, HECMW_MIN, VIS_COMM);
    HECMW_Allreduce(l_disp_max, disp_max, 5, HECMW_DOUBLE, HECMW_MAX, VIS_COMM);
  } else {
    for (i = 0; i < 5; i++) {
      disp_min[i] = l_disp_min[i];
      disp_max[i] = l_disp_max[i];
    }
  }
  return;
}

void HECMW_vis_combine(struct surface_module *sf,
                       struct hecmwST_local_mesh *mesh,
                       struct hecmwST_result_data *data, int tvertex,
                       int tpatch, int *color_list, double *minvalue,
                       double *maxvalue, Result *result, char *outfile,
                       HECMW_Comm VIS_COMM) {
  int i, j, ii;
  double color, *mivalue, *mavalue;
  int mynode, pesize;
  HECMW_Status stat;
  double *vcoord, *coord, *vcolor, *ccolor;
  int s_vertex, s_patch;
  int *c_colorid, *v_colorid;

  int *c_style, *vc_style;
  int *n_vertex, *n_patch;
  int *plist, *pplist;
  int nbase, pbase;
  int n_node, n_element;
  double r, g, b;
  FILE *outfp;
  double rgbrange[3], value;
  double range[6], minx, miny, minz, maxx, maxy, maxz, tminx, tminy, tminz,
      tmaxx, tmaxy, tmaxz, trange[6];
  double *cdisp, *vdisp, *tdisp, *tcolor, disp_min[5], disp_max[5], tmp;
  int icol, isid, isop, istyp, itopo;

  HECMW_Comm_size(VIS_COMM, &pesize);
  HECMW_Comm_rank(VIS_COMM, &mynode);
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
  rgbrange[0] = rgbrange[1] = rgbrange[2] = 1.0;
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

  if (mynode != MASTER_PE) {
    HECMW_Send(&tvertex, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);
    HECMW_Send(&tpatch, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);
  }
  if (mynode == MASTER_PE) {
    n_vertex = (int *)HECMW_calloc(pesize, sizeof(int));
    n_patch  = (int *)HECMW_calloc(pesize, sizeof(int));

    n_vertex[MASTER_PE] = tvertex;
    n_patch[MASTER_PE]  = tpatch;
    for (i = 1; i < pesize; i++) {
      HECMW_Recv(&tvertex, 1, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM, &stat);
      n_vertex[i] = tvertex;

      HECMW_Recv(&tpatch, 1, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM, &stat);
      n_patch[i] = tpatch;
    }
    n_node    = 0;
    n_element = 0;
    for (i = 0; i < pesize; i++) {
      n_node += n_vertex[i];
      n_element += n_patch[i];
    }
    for (i = 0; i < data->nn_component; i++) {
      if (color_list[i] == 1) {
        fprintf(stderr,
                "for the color component %s ,the current minimum color=%e  "
                "maximum color=%e\n",
                data->node_label[i], mivalue[i], mavalue[i]);
        /*  fprintf(stderr, "Please input the mincolor and maxcolor you would
like:\n");
scanf("%lf  %lf",&mivalue[i], &mavalue[i]);
         */
      }
    }

    if (sf[1].output_type == 1)
      strcat(outfile, ".inp");
    else if (sf[1].output_type == 2)
      strcat(outfile, ".fec");
    else if (sf[1].output_type == 4)
      strcat(outfile, ".neu");

    outfp = fopen(outfile, "w");

    if (!outfp)
      HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0009: Cannot open output file");

    if (n_node == 0) {
      fprintf(stderr, "There is no cross-section.\n");
      return;
    }

    if (sf[1].output_type == 1) {
      if (sf[1].deform_display_on == 0)
        fprintf(outfp, "%d %d 1 0 0\n", n_node, n_element);
      else
        fprintf(outfp, "%d %d 2 0 0\n", n_node, n_element);

    } else if (sf[1].output_type == 2)
      fprintf(outfp, "%d %d %d\n", n_node, n_element, 0);
    else if (sf[1].output_type == 4) {
      put_neutral_head(outfp);
      /* start writing material data */
      /* put_neutral_601(outfp, mesh); */
      put_neutral_402(outfp, mesh);
      fprintf(outfp, "   -1\n");
      fprintf(outfp, "   403\n");
    }
  }

  if (mynode != 0) {
    if (tvertex > 0) {
      vcoord = (double *)HECMW_calloc(tvertex * 3, sizeof(double));
      if (vcoord == NULL) HECMW_vis_memory_exit("vcoord");

      s_vertex = 0;
      for (ii = 0; ii < sf[0].surface_style; ii++) {
        if (result[ii].n_vertex > 0) {
          for (i = 0; i < result[ii].n_vertex; i++) {
            for (j                           = 0; j < 3; j++)
              vcoord[(i + s_vertex) * 3 + j] = result[ii].vertex[i * 3 + j];
          }
          s_vertex += result[ii].n_vertex;
        }
      }
      HECMW_Send(vcoord, tvertex * 3, HECMW_DOUBLE, MASTER_PE, 0, VIS_COMM);
      /*
outfp = fopen(outfile2, "w");
if(outfp==NULL) {
fprintf(stderr, "Couldn't create the tmp files for results: %s---\n", outfile2);
exit(0);
}
for(i=0;i<tvertex;i++)
fprintf(outfp, "%lf %lf %lf\n", vcoord[i*3], vcoord[i*3+1], vcoord[i*3+2]);
fclose(outfp);
       */
    }
  }

  if (mynode == 0) {
    if (n_vertex[0] > 0) {
      vcoord = (double *)HECMW_calloc(n_vertex[0] * 3, sizeof(double));
      if (vcoord == NULL) HECMW_vis_memory_exit("vcoord");
      s_vertex = 0;
      for (ii = 0; ii < sf[0].surface_style; ii++) {
        if (result[ii].n_vertex > 0) {
          for (i = 0; i < result[ii].n_vertex; i++) {
            if (sf[1].output_type == 1)
              fprintf(outfp, "%d %e %e %e\n", i + s_vertex + 1,
                      result[ii].vertex[i * 3], result[ii].vertex[i * 3 + 1],
                      result[ii].vertex[i * 3 + 2]);
            else if (sf[1].output_type == 2)
              fprintf(outfp, "%e %e %e\n", result[ii].vertex[i * 3],
                      result[ii].vertex[i * 3 + 1],
                      result[ii].vertex[i * 3 + 2]);
            else if (sf[1].output_type == 4)
              fprintf(outfp, "%d,0,0,1,46,0,0,0,0,0,0,%e,%e,%e,\n",
                      i + s_vertex + 1, result[ii].vertex[i * 3],
                      result[ii].vertex[i * 3 + 1],
                      result[ii].vertex[i * 3 + 2]);
            for (j                           = 0; j < 3; j++)
              vcoord[(i + s_vertex) * 3 + j] = result[ii].vertex[i * 3 + j];
          }
          s_vertex += result[ii].n_vertex;
        }
      }
    }
    for (i = 1; i < pesize; i++) {
      if (n_vertex[i] > 0) {
        coord = (double *)HECMW_calloc(n_vertex[i] * 3, sizeof(double));
        HECMW_Recv(coord, n_vertex[i] * 3, HECMW_DOUBLE, i, HECMW_ANY_TAG,
                   VIS_COMM, &stat);
        /*  sprintf(outfile3, "%s.%d.%d.tmp", outfile1, *timestep, i);
tmpfp=fopen(outfile3, "r");
if(tmpfp==NULL) {
fprintf(stderr, "cannot open the tmpfile %s\n");
exit(0);
}
for(j=0;j<n_vertex[i];j++)
fscanf(tmpfp, "%lf %lf %lf", &coord[j*3], &coord[j*3+1],&coord[j*3+2]);
fclose(tmpfp);
         */

        nbase = 0;
        for (j = 0; j < i; j++) nbase += n_vertex[j];
        for (j = 0; j < n_vertex[i]; j++) {
          if (sf[1].output_type == 1)
            fprintf(outfp, "%d %e %e %lf\n", nbase + j + 1, coord[j * 3],
                    coord[j * 3 + 1], coord[j * 3 + 2]);
          else if (sf[1].output_type == 2)
            fprintf(outfp, "%e %e %e\n", coord[j * 3], coord[j * 3 + 1],
                    coord[j * 3 + 2]);
          else if (sf[1].output_type == 4)
            fprintf(outfp, "%d,0,0,1,46,0,0,0,0,0,0,%e,%e,%e,\n", nbase + j + 1,
                    coord[j * 3], coord[j * 3 + 1], coord[j * 3 + 2]);
        }
        HECMW_free(coord);
      }
    }
  }
  /* send patches----- */
  if (mynode != 0) {
    if (tpatch > 0) {
      plist = (int *)HECMW_calloc(tpatch * 3, sizeof(int));
      if (plist == NULL) HECMW_vis_memory_exit("plist");
      s_patch = 0;
      for (ii = 0; ii < sf[0].surface_style; ii++) {
        if (result[ii].n_patch > 0) {
          for (i = 0; i < result[ii].n_patch; i++) {
            for (j                         = 0; j < 3; j++)
              plist[(i + s_patch) * 3 + j] = result[ii].patch[i * 3 + j];
          }
          s_patch += result[ii].n_patch;
        }
      }

      HECMW_Send(plist, tpatch * 3, HECMW_INT, MASTER_PE, 0, VIS_COMM);
    }
  }

  if (mynode == 0) {
    if (sf[1].output_type == 4) {
      fprintf(outfp, "   -1\n");
      fprintf(outfp, "   -1\n");
      fprintf(outfp, "   404\n");
      isid                 = mesh->section_ID[0];
      isop                 = mesh->section->sect_opt[isid];
      istyp                = 25;
      if (isop == 1) istyp = 19;
      if (isop == 2) istyp = 35;
      itopo                = 2;
      icol                 = 124;
    }

    if (n_patch[0] > 0) {
      plist = (int *)HECMW_calloc(n_patch[0] * 3, sizeof(int));
      if (plist == NULL) HECMW_vis_memory_exit("plist");
      /*	for(i=0;i<n_patch[0];i++) {
fscanf(pfile, "%d %d %d %d", &id,&v1, &v2, &v3);
if(sff->output_type==1)
fprintf(FP, "%d 1 tri %d %d %d\n", id,v1, v2, v3);
else if(sff->output_type==2)
fprintf(FP, "%d %d %d\n", v1, v2, v3);

plist[i*3]=v1;
plist[i*3+1]=v2;
plist[i*3+2]=v3;
}
       */
      s_patch = 0;
      for (ii = 0; ii < sf[0].surface_style; ii++) {
        if (result[ii].n_patch > 0) {
          for (i = 0; i < result[ii].n_patch; i++) {
            if (sf[1].output_type == 1)
              fprintf(outfp, "%d 1 tri %d %d %d\n", i + s_patch + 1,
                      result[ii].patch[i * 3], result[ii].patch[i * 3 + 1],
                      result[ii].patch[i * 3 + 2]);
            else if (sf[1].output_type == 2)
              fprintf(outfp, "%d %d %d\n", result[ii].patch[i * 3],
                      result[ii].patch[i * 3 + 1], result[ii].patch[i * 3 + 2]);
            else if (sf[1].output_type == 4) {
              fprintf(outfp, "%d, %d,%d,%d,%d,1,0,0,0,0,0,0,0,\n",
                      i + s_patch + 1, icol, isid, istyp, itopo);
              fprintf(outfp, "%8d, %8d, %8d, 0,0,0,0,0,0,0,\n",
                      result[ii].patch[i * 3], result[ii].patch[i * 3 + 1],
                      result[ii].patch[i * 3 + 2]);
              fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
              fprintf(outfp, "0,0,0,\n");
              fprintf(outfp, "0,0,0,\n");
              fprintf(outfp, "0,0,0,\n");
              fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,\n");
            }

            for (j                         = 0; j < 3; j++)
              plist[(i + s_patch) * 3 + j] = result[ii].patch[i * 3 + j];
          }
          s_patch += result[ii].n_patch;
        }
      }
    }

    for (i = 1; i < pesize; i++) {
      if (n_patch[i] > 0) {
        pplist = (int *)HECMW_calloc(n_patch[i] * 3, sizeof(int));
        HECMW_Recv(pplist, n_patch[i] * 3, HECMW_INT, i, HECMW_ANY_TAG,
                   VIS_COMM, &stat);
        nbase = 0;
        pbase = 0;
        for (j = 0; j < i; j++) {
          nbase += n_vertex[j];
          pbase += n_patch[j];
        }

        for (j = 0; j < n_patch[i]; j++) {
          if (sf[1].output_type == 1)
            fprintf(outfp, "%d 1 tri %d %d %d\n", pbase + j + 1,
                    pplist[j * 3] + nbase, pplist[j * 3 + 1] + nbase,
                    pplist[j * 3 + 2] + nbase);
          else if (sf[1].output_type == 2)
            fprintf(outfp, "%d %d %d\n", pplist[j * 3] + nbase,
                    pplist[j * 3 + 1] + nbase, pplist[j * 3 + 2] + nbase);
          else if (sf[1].output_type == 4) {
            fprintf(outfp, "%d, icol,%d,%d,%d,%d,1,0,0,0,0,0,0,0,\n",
                    pbase + j + 1, icol, isid, istyp, itopo);
            fprintf(outfp, "%d, %d, %d, 0,0,0,0,0,0,0,0,0,0,0,0,\n",
                    pplist[j * 3] + nbase, pplist[j * 3 + 1] + nbase,
                    pplist[j * 3 + 2] + nbase);
            fprintf(outfp, "0,0,0,\n");
            fprintf(outfp, "0,0,0,\n");
            fprintf(outfp, "0,0,0,\n");
            fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,\n");
          }
        }

        HECMW_free(pplist);
      }
    }
  }

  if (mynode == 0) {
    if (sf[1].output_type == 4) {
      put_neutral_middle(outfp);
    }
  }

  if (mynode != 0) {
    if (tvertex > 0) {
      vcolor = (double *)HECMW_calloc(tvertex, sizeof(double));
      if (sf[1].deform_display_on == 1)
        vdisp = (double *)HECMW_calloc(tvertex * 3, sizeof(double));

      c_style   = (int *)HECMW_calloc(tvertex, sizeof(int));
      c_colorid = (int *)HECMW_calloc(tvertex, sizeof(int));
      if ((vcolor == NULL) || (c_style == NULL) || (c_colorid == NULL))
        HECMW_vis_memory_exit("vcolor, c_style and c_colorid");
      s_vertex = 0;
      for (ii = 0; ii < sf[0].surface_style; ii++) {
        if (result[ii].n_vertex > 0) {
          for (i = 0; i < result[ii].n_vertex; i++) {
            vcolor[i + s_vertex] = result[ii].color[i];
            if (sf[ii + 1].deform_display_on == 1) {
              vdisp[(i + s_vertex) * 3]     = result[ii].disp[i * 3];
              vdisp[(i + s_vertex) * 3 + 1] = result[ii].disp[i * 3 + 1];
              vdisp[(i + s_vertex) * 3 + 2] = result[ii].disp[i * 3 + 2];
            }
            c_style[i + s_vertex]   = sf[ii + 1].display_method;
            c_colorid[i + s_vertex] = sf[ii + 1].color_comp;
          }
          s_vertex += result[ii].n_vertex;
        }
      }
      HECMW_Send(vcolor, tvertex, HECMW_DOUBLE, MASTER_PE, 0, VIS_COMM);
      HECMW_Send(c_style, tvertex, HECMW_INT, MASTER_PE, 0, VIS_COMM);
      HECMW_Send(c_colorid, tvertex, HECMW_INT, MASTER_PE, 0, VIS_COMM);
      if (sf[1].deform_display_on == 1)
        HECMW_Send(vdisp, tvertex * 3, HECMW_DOUBLE, MASTER_PE, 0, VIS_COMM);
    }

    if (tvertex > 0) {
      HECMW_free(vcoord);
      HECMW_free(vcolor);
      HECMW_free(c_style);
      HECMW_free(c_colorid);
      if (sf[1].deform_display_on == 1) HECMW_free(vdisp);
    }
    if (tpatch > 0) HECMW_free(plist);
  }
  if (mynode == 0) {
    if (sf[1].output_type == 1) {
      if (sf[1].deform_display_on == 0) {
        fprintf(outfp, "1 1\n");
        fprintf(outfp, "data, unit_unknown\n");
      } else {
        fprintf(outfp, "2 1 3\n");
        fprintf(outfp, "data, unit_unknown\n");
        fprintf(outfp, "disp, unit_unknown\n");
      }
    }
    if (sf[1].output_type == 4) {
#ifdef old_version
      fprintf(outfp, "   -1\n");
      fprintf(outfp, "   -1\n");
      fprintf(outfp, "   409\n");
      fprintf(outfp, "1\n");
      fprintf(outfp, "Default XY View\n");
      fprintf(outfp, "2,0,1,\n");
      fprintf(outfp, "35.2644,-45.,0.,\n");
      fprintf(outfp, "%e, %e, %e,\n", (tminx + tmaxx) / 2.0,
              (tminy + tmaxy) / 2.0, (tminz + tmaxz) / 2.0);
      fprintf(outfp, "0.90909,1.,0,%e,%e,%e,0.,0.,0.,\n", (tminx + tmaxx) / 2.0,
              (tminy + tmaxy) / 2.0, (tminz + tmaxz) / 2.0);
      fprintf(outfp, "%e,%e,%e,\n", (tminx + tmaxx) / 2.0,
              (tminy + tmaxy) / 2.0, (tminz + tmaxz) / 2.0);
      fprintf(outfp, "1.,\n");
      fprintf(outfp, "0.,0.,1.,1.,\n");
      fprintf(outfp, "2,0,1,1,1,\n");
      fprintf(outfp, "0,0,0,1,0,0,2,1,0,4000000,\n");
      fprintf(outfp, "9,\n");
      fprintf(outfp, "0,0,0,\n");
      fprintf(outfp, "0,0,0,\n");
      fprintf(outfp, "0,0,0,\n");
      fprintf(outfp, "0,0,0,\n");
      fprintf(outfp, "0,0,0,\n");
      fprintf(outfp, "0,0,0,\n");
      fprintf(outfp, "0,0,0,\n");
      fprintf(outfp, "0,0,0,\n");
      fprintf(outfp, "0,0,0,\n");
      fprintf(outfp, "100.,100.,1,7,\n");
      fprintf(outfp, "0,1,1.,\n");
      fprintf(outfp, "0.,0.,0.,\n");
      fprintf(outfp, "0.,0.,1.,\n"); /* sec_nor*/
      fprintf(outfp, "2.,1.,70.,0.5,\n");
      fprintf(outfp, "0.,25.,0.,0.,0.,1,100.,1000.,0.,0.,0.,\n");
      fprintf(outfp, "5.,90.,10.,10.,1.,\n");
      fprintf(outfp, "4,276,0,0,0,0,0,0,0,0,0.,0.,1.,\n");
      fprintf(outfp, "0,0,0,0,0,0,14,110,\n");
      fprintf(outfp, "0,1,1,1,0,1,0,1,1,0,1,1,0,1,1,1,0,0,1,\n");
      fprintf(outfp, "0,0,0.00000001,25.,100.,0.,0.,0.,20,\n");
      fprintf(outfp, "0,1,1,0,0,1,20.,0,\n");
      fprintf(outfp, "12,\n"); /*max_lev*/
      fprintf(outfp, "0,0.,\n");
      fprintf(outfp, "0,0.,\n");
      fprintf(outfp, "0,0.,\n");
      fprintf(outfp, "0,0.,\n");
      fprintf(outfp, "0,0.,\n");
      fprintf(outfp, "0,0.,\n");
      fprintf(outfp, "0,0.,\n");
      fprintf(outfp, "0,0.,\n");
      fprintf(outfp, "0,0.,\n");
      fprintf(outfp, "0,0.,\n");
      fprintf(outfp, "0,0.,\n");
      fprintf(outfp, "0,0.,\n");
      fprintf(outfp, "0,5,0,0,0,0.,25.,\n");
      fprintf(outfp, "4,16408,20,16504,100,16488,\n");
      fprintf(outfp, "0.,0.,\n");
      fprintf(outfp, "0.,0.,0.,0.,\n");
      fprintf(outfp, "9,\n"); /*max_xy*/
      fprintf(outfp, "1.,\n");
      fprintf(outfp, "1.,\n");
      fprintf(outfp, "1.,\n");
      fprintf(outfp, "1.,\n");
      fprintf(outfp, "1.,\n");
      fprintf(outfp, "1.,\n");
      fprintf(outfp, "1.,\n");
      fprintf(outfp, "1.,\n");
      fprintf(outfp, "1.,\n");
      fprintf(outfp, "2,\n"); /*max_xyt*/
      fprintf(outfp, "<NULL>\n");
      fprintf(outfp, "<NULL>\n");
      fprintf(outfp, "0,0,0,0,\n");
      fprintf(outfp, "0.,0.,0.,0.,\n");
      fprintf(outfp, "0.,0.,0.,0.,\n");
      fprintf(outfp, "90,\n");
      fprintf(outfp, "1,124,1,0,\n");
      fprintf(outfp, "0,60,1,1,\n");
      fprintf(outfp, "0,24,0,0,\n");
      fprintf(outfp, "0,120,1,0,\n");
      fprintf(outfp, "0,60,0,1,\n");
      fprintf(outfp, "0,24642,0,1,\n");
      fprintf(outfp, "0,124,0,1,\n");
      fprintf(outfp, "0,46,0,1,\n");
      fprintf(outfp, "0,120,1,1,\n");
      fprintf(outfp, "0,124,0,1,\n");
      fprintf(outfp, "0,124,0,0,\n");
      fprintf(outfp, "1,14,0,1,\n");
      fprintf(outfp, "0,62,0,0,\n");
      fprintf(outfp, "0,62,0,0,\n");
      fprintf(outfp, "0,10,1,1,\n");
      fprintf(outfp, "0,52,1,1,\n");
      fprintf(outfp, "0,4,1,1,\n");
      fprintf(outfp, "0,120,1,1,\n");
      fprintf(outfp, "0,12,1,1,\n");
      fprintf(outfp, "0,2,1,1,\n");
      fprintf(outfp, "0,120,1,1,\n");
      fprintf(outfp, "0,8312,1,1,\n");
      fprintf(outfp, "0,24600,0,0,\n");
      fprintf(outfp, "0,0,0,0,\n");
      fprintf(outfp, "1,123,0,1,\n");
      fprintf(outfp, "0,0,0,0,\n");
      fprintf(outfp, "2,124,0,1,\n");
      fprintf(outfp, "0,24636,0,0,\n");
      fprintf(outfp, "0,124,0,0,\n");
      fprintf(outfp, "0,4,0,0,\n");
      fprintf(outfp, "0,100,0,0,\n");
      fprintf(outfp, "0,124,0,0,\n");
      fprintf(outfp, "0,60,0,0,\n");
      fprintf(outfp, "0,56,0,0,\n");
      fprintf(outfp, "0,24,0,0,\n");
      fprintf(outfp, "0,8216,0,1,\n");
      fprintf(outfp, "0,4,0,0,\n");
      fprintf(outfp, "0,124,2,0,\n");
      fprintf(outfp, "0,0,1,1,\n");
      fprintf(outfp, "0,0,0,1,\n");
      fprintf(outfp, "1,124,5,1,\n");
      fprintf(outfp, "0,0,0,1,\n");
      fprintf(outfp, "0,24,0,1,\n");
      fprintf(outfp, "0,124,0,0,\n");
      fprintf(outfp, "0,100,0,1,\n");
      fprintf(outfp, "1,100,0,1,\n");
      fprintf(outfp, "0,0,0,1,\n");
      fprintf(outfp, "0,16,0,0,\n");
      fprintf(outfp, "0,124,4,1,\n");
      fprintf(outfp, "0,62,0,0,\n");
      fprintf(outfp, "2,124,1,1,\n");
      fprintf(outfp, "1,8254,0,0,\n");
      fprintf(outfp, "0,124,1,1,\n");
      fprintf(outfp, "1,0,5,1,\n");
      fprintf(outfp, "0,124,0,1,\n");
      fprintf(outfp, "0,100,0,1,\n");
      fprintf(outfp, "0,100,0,1,\n");
      fprintf(outfp, "1,46,0,1,\n");
      fprintf(outfp, "1,120,0,1,\n");
      fprintf(outfp, "1,4,0,1,\n");
      fprintf(outfp, "1,52,0,1,\n");
      fprintf(outfp, "1,24,0,1,\n");
      fprintf(outfp, "1,93,0,1,\n");
      fprintf(outfp, "1,12,0,1,\n");
      fprintf(outfp, "1,10,0,1,\n");
      fprintf(outfp, "1,104,0,1,\n");
      fprintf(outfp, "0,100,1,1,\n");
      fprintf(outfp, "0,24,1,1,\n");
      fprintf(outfp, "0,60,1,1,\n");
      fprintf(outfp, "0,104,1,1,\n");
      fprintf(outfp, "0,0,0,0,\n");
      fprintf(outfp, "0,90,1,1,\n");
      fprintf(outfp, "0,14,1,1,\n");
      fprintf(outfp, "0,8261,1,1,\n");
      fprintf(outfp, "0,19,0,0,\n");
      fprintf(outfp, "1,24,1,1,\n");
      fprintf(outfp, "0,4,0,0,\n");
      fprintf(outfp, "0,0,1,0,\n");
      fprintf(outfp, "0,0,0,0,\n");
      fprintf(outfp, "0,8201,0,1,\n");
      fprintf(outfp, "0,0,0,0,\n");
      fprintf(outfp, "0,0,1,1,\n");
      fprintf(outfp, "0,0,1,1,\n");
      fprintf(outfp, "0,0,1,1,\n");
      fprintf(outfp, "0,0,1,1,\n");
      fprintf(outfp, "0,0,1,1,\n");
      fprintf(outfp, "0,0,1,1,\n");
      fprintf(outfp, "0,0,1,1,\n");
      fprintf(outfp, "0,0,1,1,\n");
      fprintf(outfp, "0,0,1,1,\n");
      fprintf(outfp, "-1,\n");
      fprintf(outfp, "   -1\n");
#endif
    }
    if (n_vertex[0] > 0) {
      if ((sf[1].output_type == 4) && (sf[1].deform_display_on == 1)) {
        fprintf(outfp, "   -1\n");
        fprintf(outfp, "   451\n");
        /* first find maximum data value */
        find_minmax_disp(sf, result, VIS_COMM, disp_min, disp_max, pesize);
        tcolor = (double *)HECMW_calloc(n_node, sizeof(double));
        tdisp  = (double *)HECMW_calloc(n_node * 4, sizeof(double));
      }
      if ((sf[1].output_type == 4) && (sf[1].deform_display_on != 1)) {
        fprintf(outfp, "   -1\n");
        fprintf(outfp, "   451\n");
        fprintf(outfp, "1,1,1,\n");
        fprintf(outfp, "%s\n", sf[ii].color_comp_name);
        tmp = fabs(mavalue[sf[1].color_comp]);
        if (fabs(mivalue[sf[1].color_comp]) > tmp)
          tmp = fabs(mivalue[sf[1].color_comp]);
        fprintf(outfp, "%e %e %e,\n", mivalue[sf[1].color_comp],
                mavalue[sf[1].color_comp], tmp);
        fprintf(outfp, "1,0,0,0,0,0,0,0,0,0,\n");
        fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
        fprintf(outfp, "0,0,1,7,\n");
        fprintf(outfp, "1,1,1,\n");
      }
      vcolor    = (double *)HECMW_calloc(n_vertex[0], sizeof(double));
      c_style   = (int *)HECMW_calloc(n_vertex[0], sizeof(int));
      c_colorid = (int *)HECMW_calloc(n_vertex[0], sizeof(int));
      if (sf[1].deform_display_on == 1)
        vdisp = (double *)HECMW_calloc(n_vertex[0] * 3, sizeof(double));

      /*  for(i=0;i<n_vertex[0];i++) {
fscanf(cfile, "%d %lf %d %d ", &id,&value, &cstyle, &colorid);
       */
      s_vertex = 0;
      for (ii = 0; ii < sf[0].surface_style; ii++) {
        if (result[ii].n_vertex > 0) {
          for (i = 0; i < result[ii].n_vertex; i++) {
            vcolor[i + s_vertex]    = result[ii].color[i];
            c_style[i + s_vertex]   = sf[ii + 1].display_method;
            c_colorid[i + s_vertex] = sf[ii + 1].color_comp;
            if (sf[ii + 1].deform_display_on == 1) {
              vdisp[(i + s_vertex) * 3]     = result[ii].disp[i];
              vdisp[(i + s_vertex) * 3 + 1] = result[ii].disp[i * 3 + 1];
              vdisp[(i + s_vertex) * 3 + 2] = result[ii].disp[i * 3 + 2];
            }

            if (sf[ii + 1].display_method != 4) {
              if (sf[1].normalize_flag == 1) {
                color = (result[ii].color[i] - mivalue[sf[ii + 1].color_comp]) /
                        (mavalue[sf[ii + 1].color_comp] -
                         mivalue[sf[ii + 1].color_comp]);
                if (color < 0.0) color = 0.0;
                if (color > 1.0) color = 1.0;
              } else {
                color = result[ii].color[i];
                if (color < mivalue[sf[ii + 1].color_comp])
                  color = mivalue[sf[ii + 1].color_comp];
                if (color > mavalue[sf[ii + 1].color_comp])
                  color = mavalue[sf[ii + 1].color_comp];
              }
            } else if (sf[ii + 1].display_method == 4)
              color = result[ii].color[i];
            if (sf[1].output_type == 1) {
              if (sf[ii + 1].deform_display_on == 0)
                fprintf(outfp, "%d %e\n", i + s_vertex + 1, color);
              else if (sf[ii + 1].deform_display_on == 1)
                fprintf(outfp, "%d %e %e %e %e\n", i + s_vertex + 1, color,
                        result[ii].disp[i * 3], result[ii].disp[i * 3 + 1],
                        result[ii].disp[i * 3 + 2]);

            } else if (sf[ii + 1].output_type == 4) {
              if (sf[ii + 1].deform_display_on == 0)
                fprintf(outfp, "%d, %e,\n", i + s_vertex + 1,
                        result[ii].color[i]);
              else if (sf[ii + 1].deform_display_on == 1) {
                tcolor[i + s_vertex]          = result[ii].color[i];
                tdisp[(i + s_vertex) * 4 + 1] = result[ii].disp[i * 3];
                tdisp[(i + s_vertex) * 4 + 2] = result[ii].disp[i * 3 + 1];
                tdisp[(i + s_vertex) * 4 + 3] = result[ii].disp[i * 3 + 2];
                tdisp[(i + s_vertex) * 4] =
                    sqrt(tdisp[(i + s_vertex) * 4 + 1] *
                             tdisp[(i + s_vertex) * 4 + 1] +
                         tdisp[(i + s_vertex) * 4 + 2] *
                             tdisp[(i + s_vertex) * 4 + 2] +
                         tdisp[(i + s_vertex) * 4 + 3] *
                             tdisp[(i + s_vertex) * 4 + 3]);
              }
            } else if (sf[1].output_type == 2) {
              if (color <= 0.25 * rgbrange[0]) {
                r = g = 0.0;
                b     = (0.5 * rgbrange[0] - color) * 2 / rgbrange[0];
              } else if ((color > 0.25 * rgbrange[0]) &&
                         (color <= 0.5 * rgbrange[0])) {
                r = 0.0;
                g = (color - 0.25 * rgbrange[0]) * 2 / rgbrange[0];
                b = (0.5 * rgbrange[0] - color) * 2 / rgbrange[0];
              } else if ((color > 0.5 * rgbrange[0]) &&
                         (color <
                          ((1 - rgbrange[0] * 0.5) / 2 + 0.5 * rgbrange[0]))) {
                r = (color - 0.5 * rgbrange[0]) /
                    ((1 - rgbrange[0] * 0.5) / 2) * 0.5;
                g = (((1 - rgbrange[0] * 0.5) / 2 + 0.5 * rgbrange[0]) -
                     color) /
                    ((1 - rgbrange[0] * 0.5) / 2) * 0.5;
                b = 0.0;
              } else if (color >=
                         ((1 - rgbrange[0] * 0.5) / 2 + 0.5 * rgbrange[0])) {
                r = (color - 0.5 * rgbrange[0]) /
                    ((1 - rgbrange[0] * 0.5) / 2) * 0.5;
                g = b = 0.0;
              }

              /*       fprintf(FP, "%d %lf\n",i+1, color);
							 */ fprintf(
                  outfp, "%e %e %e\n", r, g, b);
            }
          }
          s_vertex += result[ii].n_vertex;
        }
      }
    }
    for (i = 1; i < pesize; i++) {
      if (n_vertex[i] > 0) {
        ccolor    = (double *)HECMW_calloc(n_vertex[i], sizeof(double));
        vc_style  = (int *)HECMW_calloc(n_vertex[i], sizeof(int));
        v_colorid = (int *)HECMW_calloc(n_vertex[i], sizeof(int));
        if (sf[1].deform_display_on == 1)
          cdisp = (double *)HECMW_calloc(n_vertex[i] * 3, sizeof(double));

        HECMW_Recv(ccolor, n_vertex[i], HECMW_DOUBLE, i, HECMW_ANY_TAG,
                   VIS_COMM, &stat);
        HECMW_Recv(vc_style, n_vertex[i], HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM,
                   &stat);
        HECMW_Recv(v_colorid, n_vertex[i], HECMW_INT, i, HECMW_ANY_TAG,
                   VIS_COMM, &stat);
        if (sf[1].deform_display_on == 1)
          HECMW_Recv(cdisp, n_vertex[i] * 3, HECMW_DOUBLE, i, HECMW_ANY_TAG,
                     VIS_COMM, &stat);

        nbase = 0;
        for (j = 0; j < i; j++) nbase += n_vertex[j];
        for (j = 0; j < n_vertex[i]; j++) {
          value = ccolor[j];
          if (vc_style[j] != 4) {
            if (sf[1].output_type == 2) {
              color = (value - mivalue[v_colorid[j]]) /
                      (mavalue[v_colorid[j]] - mivalue[v_colorid[j]]);
              if (color > 1.0) color = 1.0;
              if (color < 0.0) color = 0.0;
            }
            if ((sf[1].output_type == 1) || (sf[1].output_type == 4)) {
              color                                    = value;
              if (color < mivalue[v_colorid[j]]) color = mivalue[v_colorid[j]];
              if (color > mavalue[v_colorid[j]]) color = mavalue[v_colorid[j]];
              if (sf[1].normalize_flag == 1)
                color = (color - mivalue[v_colorid[j]]) /
                        (mavalue[v_colorid[j]] - mivalue[v_colorid[j]]);
            }
          } else if (vc_style[j] == 4)
            color = value;
          /*	  if(color>1.0) color=1.0;
if(color<0.0) color=0.0;
           */

          if (sf[1].output_type == 2) {
            if (color <= 0.25 * rgbrange[0]) {
              r = g = 0.0;
              b     = (0.5 * rgbrange[0] - color) * 2 / rgbrange[0];
            } else if ((color > 0.25 * rgbrange[0]) &&
                       (color <= 0.5 * rgbrange[0])) {
              r = 0.0;
              g = (color - 0.25 * rgbrange[0]) * 2 / rgbrange[0];
              b = (0.5 * rgbrange[0] - color) * 2 / rgbrange[0];
            } else if ((color > 0.5 * rgbrange[0]) &&
                       (color <
                        ((1 - rgbrange[0] * 0.5) / 2 + 0.5 * rgbrange[0]))) {
              r = (color - 0.5 * rgbrange[0]) / ((1 - rgbrange[0] * 0.5) / 2) *
                  0.5;
              g = (((1 - rgbrange[0] * 0.5) / 2 + 0.5 * rgbrange[0]) - color) /
                  ((1 - rgbrange[0] * 0.5) / 2) * 0.5;
              b = 0.0;
            } else if (color >=
                       ((1 - rgbrange[0] * 0.5) / 2 + 0.5 * rgbrange[0])) {
              r = (color - 0.5 * rgbrange[0]) / ((1 - rgbrange[0] * 0.5) / 2) *
                  0.5;
              g = b = 0.0;
            }
          }
          if (sf[1].output_type == 1) {
            if (sf[1].deform_display_on == 0)
              fprintf(outfp, "%d %e\n", nbase + j + 1, color);
            else
              fprintf(outfp, "%d %e %e %e %e\n", nbase + j + 1, color,
                      cdisp[j * 3], cdisp[j * 3 + 1], cdisp[j * 3 + 2]);
          } else if (sf[1].output_type == 4) {
            if (sf[1].deform_display_on == 0)
              fprintf(outfp, "%d, %e,\n", nbase + j + 1, color);
            else if (sf[1].deform_display_on == 1) {
              tcolor[nbase + j]          = ccolor[j];
              tdisp[(nbase + j) * 4 + 1] = cdisp[j * 3];
              tdisp[(nbase + j) * 4 + 2] = cdisp[j * 3 + 1];
              tdisp[(nbase + j) * 4 + 3] = cdisp[j * 3 + 2];
              tdisp[(nbase + j) * 4] =
                  sqrt(cdisp[j * 3] * cdisp[j * 3] +
                       cdisp[j * 3 + 1] * cdisp[j * 3 + 1] +
                       cdisp[j * 3 + 2] * cdisp[j * 3 + 2]);
            }
          } else if (sf[1].output_type == 2)
            fprintf(outfp, "%e %e %e\n", r, g, b);
        }
        HECMW_free(ccolor);
        HECMW_free(vc_style);
        HECMW_free(v_colorid);
        if (sf[1].deform_display_on == 1) HECMW_free(cdisp);
      }
    }

    if (n_vertex[0] > 0) {
      HECMW_free(vcoord);

      HECMW_free(vcolor);
      HECMW_free(c_style);
      HECMW_free(c_colorid);
      if (sf[1].deform_display_on == 1) HECMW_free(vdisp);
    }
    if (n_patch[0] > 0) HECMW_free(plist);
  }
  if (mynode == 0) {
    if ((sf[1].output_type == 4) && (sf[1].deform_display_on == 0))
      fprintf(outfp, "   -1\n");
    else if ((sf[1].output_type == 4) && (sf[1].deform_display_on == 1)) {
      for (j = 0; j < 5; j++) {
      }
    }
  }
  if (mynode == 0) fclose(outfp);
  return;
}

void put_neutral_head(FILE *outfp) {
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   100\n");
  fprintf(outfp, "<NULL>\n");
  fprintf(outfp, "8.2,\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   405\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   475\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   410\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   413\n");
  fprintf(outfp, "1,124,\n");
  fprintf(outfp, "default layer\n");
  fprintf(outfp, "9999,124,\n");
  fprintf(outfp, "<NULL>\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   570\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   571\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   572\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   573\n");
  fprintf(outfp, "   -1\n");
}

void put_neutral_601(FILE *outfp, struct hecmwST_local_mesh *mesh) {
  int i;
  int im;
  double ee, pp, rh0, alfa, gg;
  double rdum[10];

  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   601\n");
  for (im = 0; im < mesh->material->n_mat; im++) {
    fprintf(outfp, "%6d", im + 1);
    fprintf(outfp, ",-601,55,0,0,1,0,\n");
    fprintf(outfp, "<NULL>\n");
    fprintf(outfp, "10,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "25,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,\n");
    fprintf(outfp, "200,\n");
    ee      = mesh->material->mat_val[im * 4];
    pp      = mesh->material->mat_val[im * 4 + 1];
    rh0     = mesh->material->mat_val[im * 4 + 2];
    alfa    = mesh->material->mat_val[im * 4 + 3];
    gg      = ee * 0.5 / (1.0 + pp);
    rdum[0] = ee;
    rdum[1] = ee;
    rdum[2] = ee;
    rdum[3] = rdum[4] = rdum[5] = gg;
    rdum[6] = rdum[7] = rdum[8] = pp;
    rdum[9]                     = 0.0;
    for (i = 0; i < 10; i++) fprintf(outfp, "%9.2e,", rdum[i]);
    fprintf(outfp, "\n");
    fprintf(outfp, "0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,\n");
    fprintf(outfp, "0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,\n");
    rdum[0] = 0.0;
    rdum[1] = 0.0;
    rdum[2] = 0.0;
    rdum[3] = rdum[4] = rdum[5] = 0.0;
    rdum[6]                     = alfa;
    rdum[7] = rdum[8] = 0.0;
    rdum[9]           = alfa;
    for (i = 0; i < 10; i++) fprintf(outfp, "%9.2e,", rdum[i]);
    fprintf(outfp, "\n");
    rdum[0] = 0.0;
    rdum[1] = alfa;
    rdum[2] = 0.0;
    rdum[3] = rdum[4] = rdum[5] = 0.0;
    rdum[6]                     = 0.0;
    rdum[7] = rdum[8] = 0.0;
    rdum[9]           = rh0;
    for (i = 0; i < 10; i++) fprintf(outfp, "%9.2e,", rdum[i]);
    fprintf(outfp, "\n");

    fprintf(outfp, "0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,\n");
    fprintf(outfp, "0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,\n");
    fprintf(outfp, "0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,\n");
    fprintf(outfp, "0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,\n");
    fprintf(outfp, "0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,\n");
    fprintf(outfp, "0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,\n");
    fprintf(outfp, "0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,\n");
    fprintf(outfp, "0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,\n");
    fprintf(outfp, "0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,\n");
    fprintf(outfp, "0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,\n");
    fprintf(outfp, "0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,\n");
    fprintf(outfp, "0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,\n");
    fprintf(outfp, "0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,\n");
    fprintf(outfp, "0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,\n");
    fprintf(outfp, "0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,\n");
    fprintf(outfp, "50,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "70,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
  }
  fprintf(outfp, "   -1\n");
  return;
}

void put_neutral_402(FILE *outfp, struct hecmwST_local_mesh *mesh) {
  int im, idum[6], i;

  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   402\n");
  for (im = 0; im < mesh->section->n_sect; im++) {
    idum[0] = im + 1;
    idum[1] = 110;
    idum[2] = mesh->section->sect_mat_ID_item[im];
    if (mesh->section->sect_type[im] == 1)
      idum[3] = 25;
    else if (mesh->section->sect_type[im] == 2)
      idum[3] = 17;
    else if (mesh->section->sect_type[im] == 3)
      idum[3] = 5;
    else if (mesh->section->sect_type[im] == 4)
      idum[3] = 9;
    idum[4]   = 1;
    idum[5]   = 0;
    for (i = 0; i < 6; i++) fprintf(outfp, "%6i,", idum[i]);
    fprintf(outfp, "\n");
    fprintf(outfp, "<NULL>\n");
    fprintf(outfp, "0,0,0,0,\n");
    fprintf(outfp, "90,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,\n");
    fprintf(outfp, "0,0,\n");
    fprintf(outfp, "190,\n");
    for (i = 0; i < 38; i++) fprintf(outfp, "0.,0.,0.,0.,0.,\n");
    fprintf(outfp, "0,\n");
    fprintf(outfp, "0,\n");
  }
  fprintf(outfp, "   -1\n");

  return;
}

void put_neutral_middle(FILE *outfp) {
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "  615\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "  514\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "  506\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "  507\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "  408\n");
  fprintf(outfp, "   -1\n");

  /* put View : BLOCK NO. = 409 */
  /*      put_neutral_409(outfp);
   */

  fprintf(outfp, "   -1\n");
  fprintf(outfp, "  411\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "  420\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "  822\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "  412\n");
  fprintf(outfp, "1,1,0,0,0,\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "   -1\n");
  fprintf(outfp, "  450\n");
  fprintf(outfp, "1,\n");
  fprintf(outfp, "hecmw_FSTR_result\n");
  fprintf(outfp, "0,0,\n");
  fprintf(outfp, "0.,\n");
  fprintf(outfp, "1,\n");
  fprintf(outfp, "<NULL>\n");
  fprintf(outfp, "   -1\n");
}

void put_neutral_409(FILE *outfp) {
  /*
!C==put View : BLOCK NO. = 409 )
   */
  fprintf(outfp, "   -1\n");
  fprintf(outfp, " 409\n");
  fprintf(outfp, "1,\n");
  fprintf(outfp, "Default XY View\n");
  fprintf(outfp, "2,0,1,\n");
  /*      write(INEU,*) '35.2644,-45.,0.,'
   */
  fprintf(outfp, "0.,0.,0.,\n");
  fprintf(outfp, "2.5,1.25,1.5,\n");
  fprintf(outfp, "1.,1.,0,0.,0.,0.,0.,0.,0.,\n");
  fprintf(outfp, "1.03572,0.51035,0.,\n");
  fprintf(outfp, "1.2574,\n");
  fprintf(outfp, "0.,0.,1.,1.,\n");
  fprintf(outfp, "2,0,1,1,0,\n");
  fprintf(outfp, "-1,-1,0,1,0,1,1,60031,0,4000000,\n");
  fprintf(outfp, "9,\n");
  fprintf(outfp, "0,0,0,\n");
  fprintf(outfp, "0,0,0,\n");
  fprintf(outfp, "0,0,0,\n");
  fprintf(outfp, "0,0,0,\n");
  fprintf(outfp, "0,0,0,\n");
  fprintf(outfp, "0,0,0,\n");
  fprintf(outfp, "0,0,0,\n");
  fprintf(outfp, "0,0,0,\n");
  fprintf(outfp, "0,0,0,\n");
  fprintf(outfp, "100.,100.,1,7,\n");
  fprintf(outfp, "0,1,1.,\n");
  fprintf(outfp, "0.,0.,0.,\n");
  fprintf(outfp, "0.,0.,1.,\n");
  fprintf(outfp, "0,1,0,0,\n");
  fprintf(outfp, "0.,0.,0.,\n");
  fprintf(outfp, "1.,0.,0.,\n");
  fprintf(outfp, "0.,0.,0.,\n");
  fprintf(outfp, "0.,1.,0.,\n");
  fprintf(outfp, "2.,1.,70.,0.5,\n");
  fprintf(outfp, "0.,0.,0.,0.,0.,0,100.,1000.,0.,0.,0.,100.,-100.,1,\n");
  fprintf(outfp, "5.,90.,10.,10.,1.,\n");
  fprintf(outfp, "4,176,0,0,0,0,0,0,0,0,0.,0.,0.,\n");
  fprintf(outfp, "0,0,0,0,0,0,14,110,\n");
  fprintf(outfp, "0,1,1,1,0,1,0,1,1,0,1,1,0,1,1,1,0,0,1,\n");
  fprintf(outfp, "0,0,0.00000001,25.,100.,0.,0.,0.,20,\n");
  fprintf(outfp, "0,1,1,0,0,1,20.,0,\n");
  fprintf(outfp, "12,\n");
  fprintf(outfp, "0,0.,\n");
  fprintf(outfp, "0,0.,\n");
  fprintf(outfp, "0,0.,\n");
  fprintf(outfp, "0,0.,\n");
  fprintf(outfp, "0,0.,\n");
  fprintf(outfp, "0,0.,\n");
  fprintf(outfp, "0,0.,\n");
  fprintf(outfp, "0,0.,\n");
  fprintf(outfp, "0,0.,\n");
  fprintf(outfp, "0,0.,\n");
  fprintf(outfp, "0,0.,\n");
  fprintf(outfp, "0,0.,\n");
  fprintf(outfp, "0,5,0,0,0,0.,25.,\n");
  fprintf(outfp, "4,16408,20,16504,100,16488,\n");
  fprintf(outfp, "0.,0.,\n");
  fprintf(outfp, "0.,0.,0.,0.,\n");
  fprintf(outfp, "9,\n");
  fprintf(outfp, "1.,\n");
  fprintf(outfp, "1.,\n");
  fprintf(outfp, "1.,\n");
  fprintf(outfp, "1.,\n");
  fprintf(outfp, "1.,\n");
  fprintf(outfp, "1.,\n");
  fprintf(outfp, "1.,\n");
  fprintf(outfp, "1.,\n");
  fprintf(outfp, "1.,\n");
  fprintf(outfp, "2,\n");
  fprintf(outfp, "<NULL>\n");
  fprintf(outfp, "<NULL>\n");
  fprintf(outfp, "0,0,0,0,\n");
  fprintf(outfp, "0.,0.,0.,0.,\n");
  fprintf(outfp, "0.,0.,0.,0.,\n");
  fprintf(outfp, "90,1,124,1,0,\n");
  fprintf(outfp, "0,60,0,0,\n");
  fprintf(outfp, "0,24,0,0,\n");
  fprintf(outfp, "0,100,0,0,\n");
  fprintf(outfp, "0,2,0,0,\n");
  fprintf(outfp, "0,24580,0,0,\n");
  fprintf(outfp, "0,124,0,0,\n");
  fprintf(outfp, "0,46,0,0,\n");
  fprintf(outfp, "0,120,0,0,\n");
  fprintf(outfp, "0,124,0,1,\n");
  fprintf(outfp, "0,124,0,0,\n");
  fprintf(outfp, "0,12,0,1,\n");
  fprintf(outfp, "0,62,0,0,\n");
  fprintf(outfp, "0,62,0,0,\n");
  fprintf(outfp, "0,10,0,0,\n");
  fprintf(outfp, "0,52,0,0,\n");
  fprintf(outfp, "0,4,0,0,\n");
  fprintf(outfp, "0,120,0,0,\n");
  fprintf(outfp, "0,12,0,0,\n");
  fprintf(outfp, "0,2,0,0,\n");
  fprintf(outfp, "0,120,0,0,\n");
  fprintf(outfp, "0,8312,0,0,\n");
  fprintf(outfp, "0,24600,0,0,\n");
  fprintf(outfp, "0,0,0,0,\n");
  fprintf(outfp, "1,74,0,1,\n");
  fprintf(outfp, "0,0,0,0,\n");
  fprintf(outfp, "3,124,0,1,\n");
  fprintf(outfp, "0,24636,0,0,\n");
  fprintf(outfp, "0,0,0,0,\n");
  fprintf(outfp, "0,4,0,0,\n");
  fprintf(outfp, "0,100,0,0,\n");
  fprintf(outfp, "0,124,0,1,\n");
  fprintf(outfp, "0,60,0,1,\n");
  fprintf(outfp, "0,56,0,1,\n");
  fprintf(outfp, "0,24,0,0,\n");
  fprintf(outfp, "0,8216,0,1,\n");
  fprintf(outfp, "0,4,0,0,\n");
  fprintf(outfp, "0,124,2,0,\n");
  fprintf(outfp, "0,0,1,1,\n");
  fprintf(outfp, "0,0,0,1,\n");
  fprintf(outfp, "1,124,5,1,\n");
  fprintf(outfp, "0,0,0,1,\n");
  fprintf(outfp, "0,24,0,1,\n");
  fprintf(outfp, "0,124,0,0,\n");
  fprintf(outfp, "0,100,0,1,\n");
  fprintf(outfp, "1,100,0,1,\n");
  fprintf(outfp, "0,0,0,1,\n");
  fprintf(outfp, "0,16,0,0,\n");
  fprintf(outfp, "0,124,4,1,\n");
  fprintf(outfp, "0,62,0,0,\n");
  fprintf(outfp, "2,124,1,1,\n");
  fprintf(outfp, "1,8254,0,0,\n");
  fprintf(outfp, "0,124,1,1,\n");
  fprintf(outfp, "1,0,5,1,\n");
  fprintf(outfp, "0,124,0,1,\n");
  fprintf(outfp, "0,100,0,1,\n");
  fprintf(outfp, "0,100,0,1,\n");
  fprintf(outfp, "1,46,0,1,\n");
  fprintf(outfp, "1,120,0,1,\n");
  fprintf(outfp, "1,4,0,1,\n");
  fprintf(outfp, "1,52,0,1,\n");
  fprintf(outfp, "1,24,0,1,\n");
  fprintf(outfp, "1,93,0,1,\n");
  fprintf(outfp, "1,12,0,1,\n");
  fprintf(outfp, "1,10,0,1,\n");
  fprintf(outfp, "1,104,0,1,\n");
  fprintf(outfp, "0,100,0,0,\n");
  fprintf(outfp, "0,24,0,0,\n");
  fprintf(outfp, "0,60,0,0,\n");
  fprintf(outfp, "0,104,0,0,\n");
  fprintf(outfp, "0,0,0,0,\n");
  fprintf(outfp, "0,0,1,1,\n");
  fprintf(outfp, "0,0,1,1,\n");
  fprintf(outfp, "0,0,1,1,\n");
  fprintf(outfp, "0,0,1,1,\n");
  fprintf(outfp, "0,0,1,1,\n");
  fprintf(outfp, "0,4,0,0,\n");
  fprintf(outfp, "0,0,1,0,\n");
  fprintf(outfp, "0,0,0,0,\n");
  fprintf(outfp, "0,0,1,1,\n");
  fprintf(outfp, "0,0,1,1,\n");
  fprintf(outfp, "0,0,1,1,\n");
  fprintf(outfp, "0,0,1,1,\n");
  fprintf(outfp, "0,0,1,1,\n");
  fprintf(outfp, "0,0,1,1,\n");
  fprintf(outfp, "0,0,1,1,\n");
  fprintf(outfp, "0,62,1,1,\n");
  fprintf(outfp, "0,60,4,0,\n");
  fprintf(outfp, "0,0,1,1,\n");
  fprintf(outfp, "0,0,1,1,\n");
  fprintf(outfp, "-1,\n");
  fprintf(outfp, "   -1\n");

  return;
}
