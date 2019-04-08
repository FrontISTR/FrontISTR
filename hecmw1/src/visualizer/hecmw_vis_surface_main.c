/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include "hecmw_vis_surface_main.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hecmw_vis_mem_util.h"
#include "hecmw_vis_surface_compute.h"
#include "hecmw_vis_connectivity_build.h"
#include "hecmw_vis_rendering.h"
#include "hecmw_vis_combine.h"
#include "hecmw_fstr_output_femap.h"
#include "hecmw_fstr_output_vtk.h"
#include "hecmw_malloc.h"

static void conv_compname_integer(struct surface_module *sf,
                                  struct hecmwST_result_data *data);
static void chk_sf_param(struct surface_module *sf,
                         struct hecmwST_result_data *data, int ii);
static void copy_control_para(Surface *sff, struct surface_module *sf, int ii);
static void read_surface_equation(Surface *sff, struct surface_module *sf,
                                  int ii);
static void shell2hexa(struct hecmwST_local_mesh *mesh,
                       struct hecmwST_result_data *data, HECMW_Comm VIS_COMM);

Connect_inf *global_connect;

void HECMW_vis_psf_rendering(struct hecmwST_local_mesh *mesh,
                             struct hecmwST_result_data *data, int *timestep,
                             struct surface_module *sf, Parameter_rendering *sr,
                             int stat_para[NUM_CONTROL_PSF], char *outfile1, char *body,
                             HECMW_Comm VIS_COMM) {
  int pesize, mynode;
  Surface *sff;

  int *bdflag;
  double *minvalue, *maxvalue, *mincolor, *maxcolor, mic, mac;
  int sum_v, sum_t, tvertex, tpatch;
  int ii, j, find;
  int *color_list;
  double t1, t2, t3;
  Result *result;
  int name_len;
  static int first_boundary = 1;

  HECMW_Comm_size(VIS_COMM, &pesize);
  HECMW_Comm_rank(VIS_COMM, &mynode);

  if (mynode == 0) {
    t1 = HECMW_Wtime();
    t2 = t1;
  }
  if (sf[1].output_type == 5) {
    find = 0;
    for (j = 0; j < data->nn_component; j++) {
      name_len = strlen(data->node_label[j]);
      if (strncmp("NodalSTRESS", data->node_label[j], name_len) == 0 ||
          strncmp("NodalSTRESSplus", data->node_label[j], name_len) == 0 ||
          strncmp("TEMPERATURE", data->node_label[j], name_len) == 0) {
        find = 1;
        break;
      }
    }
    if (find == 0) {
      if (mynode == 0) {
        fprintf(stderr,
                "The current library only support the FEMAP output for "
                "structure/heat analysis result.\n");
      }
      return;
    } else {
      ii = 0;
      while ((outfile1[ii] != ' ') && (ii < 100) && (outfile1[ii] != '\0')) {
        if (outfile1[ii] == '.') outfile1[ii] = '_';
        ii++;
      }
      HECMW_fstr_output_femap(mesh, data, outfile1, VIS_COMM);
      return;
    }
  }
  if ((sf[1].output_type == 6) || (sf[1].output_type == 7) ||
      (sf[1].output_type == 8) || (sf[1].output_type == 9)) {
    strcat(outfile1, ".inp");
    if (sf[1].output_type == 6) {
      HECMW_avs_output(mesh, data, outfile1, VIS_COMM);
    } else if (sf[1].output_type == 8) {
      HECMW_reorder_avs_output(mesh, data, outfile1, VIS_COMM);
    } else if (sf[1].output_type == 7) {
      HECMW_bin_avs_output(mesh, data, outfile1, VIS_COMM);
    } else if (sf[1].output_type == 9) {
      HECMW_microavs_output(mesh, data, outfile1, VIS_COMM);
    }

    return;
  }
  if (sf[1].output_type == 10) {
    char buf[16];
    sprintf(buf, "_%d.inp", mynode);
    strcat(outfile1, buf);
    HECMW_separate_avs_output(mesh, data, outfile1);
    return;
  } else if(sf[1].output_type==15) {
    HECMW_vtk_output(mesh, data, body, outfile1, VIS_COMM);
    return;
  } else if(sf[1].output_type==16) {
    HECMW_bin_vtk_output(mesh, data, body, outfile1, VIS_COMM);
    return;
  }

  if (mesh->elem_type[0] > 700) shell2hexa(mesh, data, VIS_COMM);
  conv_compname_integer(sf, data);

  bdflag = (int *)HECMW_calloc(mesh->n_elem, sizeof(int));
  if (bdflag == NULL) HECMW_vis_memory_exit("bdflag");

  tvertex    = 0;
  tpatch     = 0;
  color_list = (int *)HECMW_calloc(data->nn_component, sizeof(int));
  if (color_list == NULL) HECMW_vis_memory_exit("color_list");
  for (ii = 0; ii < data->nn_component; ii++) color_list[ii] = 0;
  for (ii = 1; ii < sf[0].surface_style + 1; ii++) {
    if (sf[ii].display_method != 4) {
      if ((sf[ii].color_comp >= 0) && (sf[ii].color_comp < data->nn_component))
        color_list[sf[ii].color_comp] = 1;
    }
  }
  minvalue = (double *)HECMW_calloc(data->nn_component, sizeof(double));
  mincolor = (double *)HECMW_calloc(data->nn_component, sizeof(double));
  maxvalue = (double *)HECMW_calloc(data->nn_component, sizeof(double));
  maxcolor = (double *)HECMW_calloc(data->nn_component, sizeof(double));
  if ((minvalue == NULL) || (mincolor == NULL) || (maxvalue == NULL) ||
      (maxcolor == NULL))
    HECMW_vis_memory_exit("minmax_color");
  for (ii = 0; ii < data->nn_component; ii++) {
    minvalue[ii] = mincolor[ii] = 1.0E17;
    maxvalue[ii] = maxcolor[ii] = -1.0E17;
  }
  mic = 1.0E17;
  mac = -1.0E17;

  sff    = (Surface *)HECMW_malloc(sizeof(Surface));
  result = (Result *)HECMW_calloc(sf[0].surface_style, sizeof(Result));
  if ((sff == NULL) || (result == NULL))
    HECMW_vis_memory_exit("sff and result");
  for (ii = 0; ii < sf[0].surface_style; ii++) {
    result[ii].n_patch  = 0;
    result[ii].n_vertex = 0;
    result[ii].vertex   = NULL;
    result[ii].color    = NULL;
    result[ii].patch    = NULL;
  }

  for (ii = 1; ii < sf[0].surface_style + 1; ii++) {
    /* check parameters */

    chk_sf_param(sf, data, ii);

    copy_control_para(sff, sf, ii);
    /*   fprintf(stderr, " coef=%lf %lf %lf %lf surface_style %d  display_way:
%d  color_comp=%d\n",sff->cont_equ[6], sff->cont_equ[7], sff->cont_equ[8],
sff->cont_equ[9],
sff->surface_style, sff->display_way, sff->color_comp);
     */
    if (mesh->elem_type[0] < 300) {
      if (sff->display_way != 4) {
        mic = mincolor[sff->color_comp];
        mac = maxcolor[sff->color_comp];
        HECMW_vis_surface_compute(sff, mesh, data, bdflag, &sum_v, &sum_t,
                                  &tvertex, &tpatch, &mic, &mac, result, ii - 1,
                                  mynode, VIS_COMM);

        if (mic < minvalue[sff->color_comp]) minvalue[sff->color_comp] = mic;
        if (mac > maxvalue[sff->color_comp]) maxvalue[sff->color_comp] = mac;
      } else if (sff->display_way == 4)
        HECMW_vis_surface_compute(sff, mesh, data, bdflag, &sum_v, &sum_t,
                                  &tvertex, &tpatch, &mic, &mac, result, ii - 1,
                                  mynode, VIS_COMM);
    } else if (mesh->elem_type[0] > 300) {
      if (sff->surface_style == 1) {
        if (first_boundary == 1) {
          global_connect = (Connect_inf *)HECMW_malloc(sizeof(Connect_inf));
          if (global_connect == NULL) HECMW_vis_memory_exit("global_connect");
        }
        if (sff->display_way != 4) {
          mic = mincolor[sff->color_comp];
          mac = maxcolor[sff->color_comp];
          /*    find_fault_surface(sff, v, &tvertex, &tpatch, &mic,&mac, result,
           * ii-1);
           */
          /*  if(strncmp(sff->group_name, "boundary", 8)==0)
           */
          HECMW_vis_find_boundary_surface(sff, mesh, data, &tvertex, &tpatch,
                                          &mic, &mac, result, ii - 1, VIS_COMM,
                                          first_boundary, global_connect);
          if (mic < minvalue[sff->color_comp]) minvalue[sff->color_comp] = mic;
          if (mac > maxvalue[sff->color_comp]) maxvalue[sff->color_comp] = mac;
        } else if (sff->display_way == 4) {
          /*  if(strncmp(sff->group_name, "boundary", 8)==0)
           */
          HECMW_vis_find_boundary_surface(sff, mesh, data, &tvertex, &tpatch,
                                          &mic, &mac, result, ii - 1, VIS_COMM,
                                          first_boundary, global_connect);
        }
        first_boundary = 0;
      } else if (sff->surface_style == 2) {
        if (mynode == 0)
          fprintf(stderr, "iso surface parameters: %d %lf\n", sff->data_comp,
                  sff->iso_value);
        HECMW_vis_chk_bounds(mesh, bdflag);
        sum_v = 0;
        sum_t = 0;
        if (sff->display_way != 4) {
          mic = mincolor[sff->color_comp];
          mac = maxcolor[sff->color_comp];
          HECMW_vis_surface_compute(sff, mesh, data, bdflag, &sum_v, &sum_t,
                                    &tvertex, &tpatch, &mic, &mac, result,
                                    ii - 1, mynode, VIS_COMM);

          if (mic < minvalue[sff->color_comp]) minvalue[sff->color_comp] = mic;
          if (mac > maxvalue[sff->color_comp]) maxvalue[sff->color_comp] = mac;
        } else if (sff->display_way == 4)
          HECMW_vis_surface_compute(sff, mesh, data, bdflag, &sum_v, &sum_t,
                                    &tvertex, &tpatch, &mic, &mac, result,
                                    ii - 1, mynode, VIS_COMM);

      } else if (sff->surface_style == 3) {
        /* get_onepoint_on_section(css, p);

for(i=0;i<3;i++)
pp[i]=p[i];

loopflag=1;
loopnum=0;
loopadd=1;  loopmin=0;
while(loopflag==1) {
loopnum++;*/
        HECMW_vis_chk_bounds(mesh, bdflag);
        sum_v = 0;
        sum_t = 0;
        if (sff->display_way != 4) {
          mic = mincolor[sff->color_comp];
          mac = maxcolor[sff->color_comp];

          HECMW_vis_surface_compute(sff, mesh, data, bdflag, &sum_v, &sum_t,
                                    &tvertex, &tpatch, &mic, &mac, result,
                                    ii - 1, mynode, VIS_COMM);

          if (mic < minvalue[sff->color_comp]) minvalue[sff->color_comp] = mic;
          if (mac > maxvalue[sff->color_comp]) maxvalue[sff->color_comp] = mac;
        } else if (sff->display_way == 4)
          HECMW_vis_surface_compute(sff, mesh, data, bdflag, &sum_v, &sum_t,
                                    &tvertex, &tpatch, &mic, &mac, result,
                                    ii - 1, mynode, VIS_COMM);
      }
    }
  }
  /*  fclose(vfile); fclose(pfile);  fclose(cfile);
   */
  HECMW_free(bdflag);
  HECMW_Barrier(VIS_COMM);
  if (mynode == MASTER_PE) {
    t3 = HECMW_Wtime();
    fprintf(stderr, "Finish visual computing, the time is: %lf\n", t3 - t2);
    t2 = t3;
  }
  if (sf[1].output_type == 4) {
    ii = 0;
    while ((outfile1[ii] != ' ') && (ii < 100) && (outfile1[ii] != '\0')) {
      if (outfile1[ii] == '.') outfile1[ii] = '_';
      ii++;
    }
  }
  if ((sf[1].output_type == 1) || (sf[1].output_type == 2) ||
      (sf[1].output_type == 4))
    HECMW_vis_combine(sf, mesh, data, tvertex, tpatch, color_list, minvalue,
                      maxvalue, result, outfile1, VIS_COMM);
  else if (sf[1].output_type == 3)
    HECMW_vis_rendering_surface(sf, sr, mesh, data, tvertex, tpatch, color_list,
                                minvalue, maxvalue, result, outfile1, stat_para,
                                VIS_COMM, *timestep);

  if (mynode == 0) {
    fprintf(stderr, "surface module finish\n");
    t3 = HECMW_Wtime();
    fprintf(stderr, "the rendering and output file time is %lf\n", t3 - t2);
  }
  for (ii = 0; ii < sf[0].surface_style; ii++) {
    if (result[ii].n_vertex > 0) {
      HECMW_free(result[ii].vertex);
      HECMW_free(result[ii].color);
    }
    if (result[ii].n_patch > 0) HECMW_free(result[ii].patch);
  }
  HECMW_free(result);
  /*
HECMW_free(sff);
   */
  HECMW_free(color_list);
  HECMW_free(mincolor);
  HECMW_free(maxcolor);
  HECMW_free(minvalue);
  HECMW_free(maxvalue);
  HECMW_Barrier(VIS_COMM);
  /*}
if(mynode==0) {
fprintf(stderr,"surface module finish\n");
t3=HECMW_Wtime();
fprintf(stderr, "the rendering and output file time is %lf\n", t3-t2);
}

HECMW_Finalize();
   */

  return;
}

static void conv_compname_integer(struct surface_module *sf,
                                  struct hecmwST_result_data *data) {
  int i, j, name_len, find;

  for (i = 1; i < sf[0].surface_style + 1; i++) {
    if (sf[i].surface_style == 2) {
      if (sf[i].data_comp < 0) {
        find = 0;
        if (strncmp(sf[i].data_comp_name, "NULL", 4) != 0) {
          for (j = 0; j < data->nn_component; j++) {
            name_len = strlen(data->node_label[j]);
            if (strncmp(sf[i].data_comp_name, data->node_label[j], name_len) ==
                0) {
              sf[i].data_comp = j;
              find            = 1;
              break;
            }
          }
        }
        if (find == 0) {
          fprintf(stderr,
                  "ERROR: HEC-MW-VIS-E1051:the name for data component is not "
                  "correct:%s\n",
                  sf[i].data_comp_name);
          HECMW_vis_print_exit("Please check it again");
        }
      }

      if (data->nn_dof[sf[i].data_comp] > 1) {
        if (sf[i].data_subcomp < 0) {
          find = 0;
          if (strncmp(sf[i].data_subcomp_name, "norm", 4) == 0) {
            sf[i].data_subcomp = 0;
            find               = 1;
          } else if (strncmp(sf[i].data_subcomp_name, "x", 1) == 0) {
            sf[i].data_subcomp = 1;
            find               = 1;
          } else if (strncmp(sf[i].data_subcomp_name, "y", 1) == 0) {
            sf[i].data_subcomp = 2;
            find               = 1;
          } else if (strncmp(sf[i].data_subcomp_name, "z", 1) == 0) {
            sf[i].data_subcomp = 3;
            find               = 1;
          } else {
            fprintf(stderr,
                    "ERROR: HEC-MW-VIS-E1052:The subcompnent name is not "
                    "correct, it must be norm, x, y, or z. \n");
            HECMW_vis_print_exit("Please check it again");
          }
        }
      }
    }

    if (sf[i].deform_display_on == 1) {
      if (sf[i].disp_comp < 0) {
        find = 0;
        if (strncmp(sf[i].disp_comp_name, "NULL", 4) != 0) {
          for (j = 0; j < data->nn_component; j++) {
            name_len = strlen(data->node_label[j]);
            if (strncmp(sf[i].disp_comp_name, data->node_label[j], name_len) ==
                0) {
              sf[i].disp_comp = j;
              find            = 1;
              break;
            }
          }
        }
        if (find == 0) {
          fprintf(stderr,
                  "ERROR: HEC-MW-VIS-E1051:the name for deformation component "
                  "is not correct:%s\n",
                  sf[i].disp_comp_name);
          HECMW_vis_print_exit("Please check it again");
        }
      }
    }

    if (sf[i].color_comp < 0) {
      find = 0;
      if (strncmp(sf[i].color_comp_name, "NULL", 4) != 0) {
        for (j = 0; j < data->nn_component; j++) {
          name_len = strlen(data->node_label[j]);
          if (strncmp(sf[i].color_comp_name, data->node_label[j], name_len) ==
              0) {
            sf[i].color_comp = j;
            find             = 1;
            break;
          }
        }
      }
      if (find == 0) {
        fprintf(stderr,
                "ERROR: HEC-MW-VIS-E1053:the name for color component is not "
                "correct:%s\n",
                sf[i].color_comp_name);
        HECMW_vis_print_exit("Please check it again");
      }
    }
    if (data->nn_dof[sf[i].color_comp] > 1) {
      if (sf[i].color_subcomp < 0) {
        find = 0;
        if (strncmp(sf[i].color_subcomp_name, "norm", 4) == 0) {
          sf[i].color_subcomp = 0;
          find                = 1;
        } else if (strncmp(sf[i].color_subcomp_name, "x", 1) == 0) {
          sf[i].color_subcomp = 1;
          find                = 1;
        } else if (strncmp(sf[i].color_subcomp_name, "y", 1) == 0) {
          sf[i].color_subcomp = 2;
          find                = 1;
        } else if (strncmp(sf[i].color_subcomp_name, "z", 1) == 0) {
          sf[i].color_subcomp = 3;
          find                = 1;
        } else {
          fprintf(stderr,
                  "ERROR: HEC-MW-VIS-E1052:The subcompnent name is not "
                  "correct, it must be norm, x, y, or z. \n");
          HECMW_vis_print_exit("Please check it again");
        }
      }
    }
  }
  return;
}

static void chk_sf_param(struct surface_module *sf,
                         struct hecmwST_result_data *data, int ii) {
  if (data->nn_component < sf[ii].color_comp) {
    HECMW_vis_print_exit(
        "ERROR: HEC-MW-VIS-E1054: color_comp is wrong: >nn_component");
  }
  if (data->nn_dof[sf[ii].color_comp] > 1) {
    if (sf[ii].color_subcomp > data->nn_dof[sf[ii].color_comp]) {
      HECMW_vis_print_exit(
          "ERROR: HEC-MW-VIS-E1055: color_subcomp is wrong: >dof");
    }
  }
  /*	if(sf[ii].surface_style==1) {
          if((sf[ii].defined_style<1) || (sf[ii].defined_style>2)){
  HECMW_vis_print_exit("Error: the style of group defined: 1 or 2");
  }
  }
   */
  if (sf[ii].surface_style == 2) {
    if (data->nn_component < sf[ii].data_comp)
      HECMW_vis_print_exit(
          "ERROR: HEC-MW-VIS-E1056: data component number is wrong: "
          ">nn_component");
    if (data->nn_dof[sf[ii].data_comp] > 1) {
      if (sf[ii].data_subcomp > data->nn_dof[sf[ii].data_comp])
        HECMW_vis_print_exit(
            "ERROR: HEC-MW-VIS-E1057:  data_subcomp is wrong: >dof");
    }
  }
  if (sf[ii].surface_style == 3) {
    if ((sf[ii].method < 1) || (sf[ii].method > 5)) {
      HECMW_vis_print_exit(
          "ERROR: HEC-MW-VIS-E1058:  the number of method, it must be between "
          "1 and 5");
    }
    if ((sf[ii].display_method < 1) || (sf[ii].display_method > 5)) {
      HECMW_vis_print_exit(
          "ERROR: HEC-MW-VIS-E1004: the number of display_method,it should be "
          "between 1 and 5");
    }
    if (sf[ii].isoline_number < 0) {
      HECMW_vis_print_exit(
          "ERROR: HEC-MW-VIS-E1059:  the number of isolines,it should be >=0");
    }
  }
  return;
}

static void copy_control_para(Surface *sff, struct surface_module *sf, int ii) {
  int j;

  sff->surface_style = sf[ii].surface_style;
  if (sff->surface_style == 1) {
    /*		sff->group_name=(char *)HECMW_calloc(100, sizeof(char));
    if(sff->group_name==NULL)
            HECMW_vis_memory_exit("group_name");
    sff->group_name=sf[ii].group_name;
    sff->defined_style=sf[ii].defined_style;
     */
  } else if (sff->surface_style == 2) {
    sff->data_comp    = sf[ii].data_comp;
    sff->data_subcomp = sf[ii].data_subcomp;
    sff->iso_value    = sf[ii].iso_value;
  } else if (sff->surface_style == 3) {
    /*  fscanf(fp1, "%d", &(css->color_comp));
     */
    /*  fscanf(fp1, "%d", &(css->cross_type));
		 */ sff->cross_type = 2;
    read_surface_equation(sff, sf, ii);
  }
  /*  fscanf(fp1,"%d",&(css->display_way));
	 */ sff->display_way = sf[ii].display_method;
  if (sff->display_way == 1) {
    sff->color_comp    = sf[ii].color_comp;
    sff->color_subcomp = sf[ii].color_subcomp;
    sff->rgbrange[0] = sff->rgbrange[1] = sff->rgbrange[2] = 1.0;
    sff->isonumber                                         = 0;
  } else if (sff->display_way == 2) {
    sff->color_comp    = sf[ii].color_comp;
    sff->color_subcomp = sf[ii].color_subcomp;
    sff->rgbrange[0] = sff->rgbrange[1] = sff->rgbrange[2] = 0.0;
    /*	fscanf(fp1,"%d", &(css->isonumber));
		  */ sff->isonumber = sf[ii].isoline_number;
  } else if (sff->display_way == 3) {
    /*	  fscanf(fp1, "%lf %lf %lf", &(css->rgbrange[0]),&(css->rgbrange[1]),&(css->rgbrange[2]));
		  */ sff->color_comp = sf[ii].color_comp;
    sff->color_subcomp               = sf[ii].color_subcomp;

    sff->rgbrange[0] = sff->rgbrange[1] = sff->rgbrange[2] = 1.0;
    sff->isonumber = sf[ii].isoline_number;

                  /*	fscanf(fp1,"%d", &(css->isonumber));
		   */  }
	 else if(sff->display_way==4) {
                    /*	  fscanf(fp1, "%lf %lf %lf", &(css->rgbrange[0]),&(css->rgbrange[1]),&(css->rgbrange[2]));
		  */ sff->specified_color = sf[ii].specified_color;
                    sff->rgbrange[0] = sff->rgbrange[1] = sff->rgbrange[2] =
                        1.0;

                  /*	fscanf(fp1,"%d", &(css->isonumber));
		   */  }
	 else if(sff->display_way==5) {
                    /*	  fscanf(fp1, "%lf %lf %lf", &(css->rgbrange[0]),&(css->rgbrange[1]),&(css->rgbrange[2]));
		  */ sff->color_comp   = sf[ii].color_comp;
                    sff->color_subcomp = sf[ii].color_subcomp;

                    sff->rgbrange[0] = sff->rgbrange[1] = sff->rgbrange[2] =
                        1.0;
                    sff->isonumber = sf[ii].isoline_number;

                  /*	fscanf(fp1,"%d", &(css->isonumber));
		   */  }
                  sff->deform_display_on = sf[ii].deform_display_on;
                  sff->output_type       = sf[ii].output_type;
                  if (sff->deform_display_on == 1) {
                    sff->disp_scale    = sf[ii].disp_scale;
                    sff->initial_style = sf[ii].initial_style;
                    sff->deform_style  = sf[ii].deform_style;
                    sff->disp_comp     = sf[ii].disp_comp;
                    for (j = 0; j < 3; j++) {
                      sff->initial_line_color[j] = sf[ii].initial_line_color[j];
                      sff->deform_line_color[j]  = sf[ii].deform_line_color[j];
                    }
                  }
}

static void read_surface_equation(Surface *sff, struct surface_module *sf,
                                  int ii) {
  int input_way;
  int i;
  double point[3], coff[10];
  double laxis[3], radius;
  for (i = 0; i < 10; i++) coff[i] = 0.0;

  input_way = sf[ii].method;
  switch (input_way) {
    case 1: /*the surface is a sphere*/
      for (i = 0; i < 3; i++) point[i] = sf[ii].point[i];
      radius                           = sf[ii].radius;
      coff[0] = coff[1] = coff[2] = 1;
      coff[6]                     = -2 * point[0];
      coff[7]                     = -2 * point[1];
      coff[8]                     = -2 * point[2];
      coff[9]                     = point[0] * point[0] + point[1] * point[1] +
                point[2] * point[2] - radius * radius;
      break;

    case 2: /*the surface is an ellipsoidal surface*/
      for (i = 0; i < 3; i++) {
        point[i] = sf[ii].point[i];
        laxis[i] = sf[ii].length[i];
      }
      for (i = 0; i < 3; i++) coff[i] = 1 / (laxis[i] * laxis[i]);
      for (i        = 0; i < 3; i++)
        coff[i + 6] = -2.0 * point[i] / (laxis[i] * laxis[i]);
      coff[9]       = point[0] * point[0] / (laxis[0] * laxis[0]) +
                (point[1] * point[1]) / (laxis[1] * laxis[1]) +
                (point[2] * point[2]) / (laxis[2] * laxis[2]) - 1;

      break;

    case 3: /*the surface is a hyperboloid surface*/
      for (i = 0; i < 3; i++) {
        point[i] = sf[ii].point[i];
        laxis[i] = sf[ii].length[i];
      }

      for (i = 0; i < 3; i++) coff[i] = 1 / (laxis[i] * laxis[i]);
      coff[2]                         = -coff[2];
      for (i        = 0; i < 3; i++)
        coff[i + 6] = -2.0 * point[i] / (laxis[i] * laxis[i]);
      coff[8]       = -coff[8];
      coff[9]       = point[0] * point[0] / (laxis[0] * laxis[0]) +
                (point[1] * point[1]) / (laxis[1] * laxis[1]) -
                (point[2] * point[2]) / (laxis[2] * laxis[2]) - 1;
      break;

    case 4: /*the surface is parabolic*/
      for (i = 0; i < 3; i++) {
        point[i] = sf[ii].point[i];
        laxis[i] = sf[ii].length[i];
      }
      for (i = 0; i < 2; i++) coff[i] = 1 / (laxis[i] * laxis[i]);
      coff[1]                         = -coff[1];
      coff[2]                         = 0;
      coff[6]                         = -2.0 * point[0] / (laxis[0] * laxis[0]);
      coff[7]                         = 2.0 * point[1] / (laxis[1] * laxis[1]);
      coff[8]                         = -1;
      coff[9] = point[0] * point[0] / (laxis[0] * laxis[0]) -
                (point[1] * point[1]) / (laxis[1] * laxis[1]) + point[2];
      break;

    case 5:
      for (i = 0; i < 10; i++) coff[i] = sf[ii].coef[i];
      break;
  }

  for (i = 0; i < 10; i++) sff->cont_equ[i] = coff[i];
  return;
}

#if 0
double cal_matrix(double a[3][3])
{
	double value;
	value=a[0][0]*a[1][1]*a[2][2]+a[0][1]*a[1][2]*a[2][0]+a[0][2]*a[1][0]*a[2][1]
	                                                                           -a[0][0]*a[1][2]*a[2][1]-a[0][1]*a[1][0]*a[2][2]-a[0][2]*a[1][1]*a[2][0];
	return(value);
}
#endif

#ifdef old
int chk_gid(struct visual_buf *vol, int *bdflag, int *gid, int *over_elem_id) {
  int i;
  int count;

  count = 0;
  for (i = 0; i < vol->mesh->n_elem; i++) {
    gid[i] = vol->mesh->global_elem_id[i];
    if ((bdflag[i] % 1024) != HEX_FACE_INDEX) {
      over_elem_id[count] = i;
      count++;
    }
  }

  return count;
}
#endif

#if 0
int judge_box(double box[HEX_N_NODE*3], double cont_equ[10])
{
	int loopflag1,i,sum;
	double x,y,z,f[HEX_N_NODE];
	for(i=0;i<HEX_N_NODE;i++) {
		x=box[i*3];
		y=box[i*3+1];
		z=box[i*3+2];
		f[i]=x*cont_equ[0]+y*cont_equ[1]+z*cont_equ[2]+cont_equ[3];
	}
	sum=0;
	for(i=0;i<HEX_N_NODE;i++) {
		if(f[i]>0) sum++;
		if(f[i]<0) sum--;
	}
	if((sum==-8) || (sum==8))
		loopflag1=0;
	else loopflag1=1;
	return(loopflag1);
}
#endif

#ifdef old_version
void shell2hexa(struct hecmwST_local_mesh *mesh,
                struct hecmwST_result_data *data) {
  int i, j, k;
  double *coord;
  int *connect;
  double thickness;
  double *value;
  int tn_component;
  fprintf(stderr, "Start transforming...\n");
  /*	thickness=1.0;
   */
  thickness = mesh->section->sect_R_item[0];
  coord     = (double *)HECMW_calloc(2 * mesh->n_node * 3, sizeof(double));
  if (coord == NULL) HECMW_vis_memory_exit("coord");

  for (i = 0; i < mesh->n_node; i++) {
    for (j = 0; j < 3; j++) coord[i * 2 * 3 + j] = mesh->node[i * 3 + j];
    for (j = 0; j < 2; j++) coord[(i * 2 + 1) * 3 + j] = mesh->node[i * 3 + j];
    coord[(i * 2 + 1) * 3 + 2] = mesh->node[i * 3 + 2] + thickness;
  }
  HECMW_free(mesh->node);

  mesh->node = coord;
  if ((mesh->elem_type[0] == 731) || (mesh->elem_type[0] == 732)) {
    connect = (int *)HECMW_calloc(mesh->n_elem * 6, sizeof(int));
    if (connect == NULL) HECMW_vis_memory_exit("connect");
    for (i = 0; i < mesh->n_elem; i++) {
      for (j = 0; j < 3; j++)
        connect[i * 6 + j] =
            (mesh->elem_node_item[mesh->elem_node_index[i] + j] - 1) * 2 + 1;
      for (j = 0; j < 3; j++)
        connect[i * 6 + 3 + j] =
            (mesh->elem_node_item[mesh->elem_node_index[i] + j] - 1) * 2 + 1 +
            1;
    }
    for (i = 0; i < mesh->n_elem; i++) mesh->elem_type[i] = 351;
    for (i = 0; i < mesh->n_elem + 1; i++) mesh->elem_node_index[i] = 6 * i;
    HECMW_free(mesh->elem_node_item);
    mesh->elem_node_item = connect;
  } else if ((mesh->elem_type[0] == 741) || (mesh->elem_type[0] == 742) ||
             (mesh->elem_type[0] == 743)) {
    connect = (int *)HECMW_calloc(mesh->n_elem * 8, sizeof(int));
    if (connect == NULL) HECMW_vis_memory_exit("connect");
    for (i = 0; i < mesh->n_elem; i++) {
      for (j = 0; j < 4; j++)
        connect[i * 8 + j] =
            (mesh->elem_node_item[mesh->elem_node_index[i] + j] - 1) * 2 + 1;
      for (j = 0; j < 4; j++)
        connect[i * 8 + 4 + j] =
            (mesh->elem_node_item[mesh->elem_node_index[i] + j] - 1) * 2 + 1 +
            1;
      /*			fprintf(stderr, "%d ", i);
      for(j=0;j<8;j++)
              fprintf(stderr, "%d ", connect[i*8+j]);
      fprintf(stderr, "\n");
       */
    }
    for (i = 0; i < mesh->n_elem; i++) mesh->elem_type[i] = 361;
    for (i = 0; i < mesh->n_elem + 1; i++) mesh->elem_node_index[i] = 8 * i;
    HECMW_free(mesh->elem_node_item);
    mesh->elem_node_item = connect;
  }
  tn_component = 0;
  for (i = 0; i < data->nn_component; i++) tn_component += data->nn_dof[i];
  value =
      (double *)HECMW_calloc(tn_component * mesh->n_node * 2, sizeof(double));
  if (value == NULL) HECMW_vis_memory_exit("value");
  for (i = 0; i < mesh->n_node; i++) {
    for (j = 0; j < tn_component; j++)
      value[i * 2 * tn_component + j] =
          data->node_val_item[i * tn_component + j];
    for (j = 0; j < tn_component; j++)
      value[(i * 2 + 1) * tn_component + j] =
          data->node_val_item[i * tn_component + j];
  }
  HECMW_free(data->node_val_item);
  data->node_val_item = value;

  fprintf(stderr, "It is ok to transform shell to solid\n");
  mesh->n_node *= 2;
  mesh->nn_internal *= 2;

  return;
}

#endif

static void shell2hexa(struct hecmwST_local_mesh *mesh,
                       struct hecmwST_result_data *data, HECMW_Comm VIS_COMM) {
  int i, j;
  double *coord;
  int *connect;
  double thickness;
  double *value;
  int tn_component, pesize;
  double range[6], trange[6];
  fprintf(stderr, "Start transforming...\n");
  /*	thickness=1.0;
   */

  HECMW_Comm_size(VIS_COMM, &pesize);
  range[0] = range[2] = range[4] = 1.0E17;
  range[1] = range[3] = range[5] = -1.0E17;

  for (i = 0; i < mesh->n_node; i++) {
    for (j = 0; j < 3; j++) {
      if (mesh->node[i * 3 + j] < range[j * 2])
        range[j * 2] = mesh->node[i * 3 + j];
      if (mesh->node[i * 3 + j] > range[j * 2 + 1])
        range[j * 2 + 1] = mesh->node[i * 3 + j];
    }
  }
  if (pesize > 1) {
    HECMW_Allreduce(&range[0], &trange[0], 1, HECMW_DOUBLE, HECMW_MIN,
                    VIS_COMM);
    HECMW_Allreduce(&range[1], &trange[1], 1, HECMW_DOUBLE, HECMW_MAX,
                    VIS_COMM);
    HECMW_Allreduce(&range[2], &trange[2], 1, HECMW_DOUBLE, HECMW_MIN,
                    VIS_COMM);
    HECMW_Allreduce(&range[3], &trange[3], 1, HECMW_DOUBLE, HECMW_MAX,
                    VIS_COMM);
    HECMW_Allreduce(&range[4], &trange[4], 1, HECMW_DOUBLE, HECMW_MIN,
                    VIS_COMM);
    HECMW_Allreduce(&range[5], &trange[5], 1, HECMW_DOUBLE, HECMW_MAX,
                    VIS_COMM);
  } else {
    for (i = 0; i < 6; i++) trange[i] = range[i];
  }
  thickness = 0.01 * sqrt((trange[1] - trange[0]) * (trange[1] - trange[0]) +
                          (trange[3] - trange[2]) * (trange[3] - trange[2]) +
                          (trange[5] - trange[4]) * (trange[5] - trange[4]));
  coord = (double *)HECMW_calloc(2 * mesh->n_node * 3, sizeof(double));
  if (coord == NULL) HECMW_vis_memory_exit("coord");

  for (i = 0; i < mesh->n_node; i++) {
    for (j = 0; j < 3; j++) coord[i * 2 * 3 + j] = mesh->node[i * 3 + j];
    for (j = 0; j < 2; j++) coord[(i * 2 + 1) * 3 + j] = mesh->node[i * 3 + j];
    coord[(i * 2 + 1) * 3 + 2] = mesh->node[i * 3 + 2] + thickness;
  }
  HECMW_free(mesh->node);

  mesh->node = coord;
  if ((mesh->elem_type[0] == 731) || (mesh->elem_type[0] == 732)) {
    connect = (int *)HECMW_calloc(mesh->n_elem * 6, sizeof(int));
    if (connect == NULL) HECMW_vis_memory_exit("connect");
    for (i = 0; i < mesh->n_elem; i++) {
      for (j = 0; j < 3; j++)
        connect[i * 6 + j] =
            (mesh->elem_node_item[mesh->elem_node_index[i] + j] - 1) * 2 + 1;
      for (j = 0; j < 3; j++)
        connect[i * 6 + 3 + j] =
            (mesh->elem_node_item[mesh->elem_node_index[i] + j] - 1) * 2 + 1 +
            1;
    }
    for (i = 0; i < mesh->n_elem; i++) mesh->elem_type[i] = 351;
    for (i = 0; i < mesh->n_elem + 1; i++) mesh->elem_node_index[i] = 6 * i;
    HECMW_free(mesh->elem_node_item);
    mesh->elem_node_item = connect;
  } else if ((mesh->elem_type[0] == 741) || (mesh->elem_type[0] == 742) ||
             (mesh->elem_type[0] == 743)) {
    connect = (int *)HECMW_calloc(mesh->n_elem * 8, sizeof(int));
    if (connect == NULL) HECMW_vis_memory_exit("connect");
    for (i = 0; i < mesh->n_elem; i++) {
      for (j = 0; j < 4; j++)
        connect[i * 8 + j] =
            (mesh->elem_node_item[mesh->elem_node_index[i] + j] - 1) * 2 + 1;
      for (j = 0; j < 4; j++)
        connect[i * 8 + 4 + j] =
            (mesh->elem_node_item[mesh->elem_node_index[i] + j] - 1) * 2 + 1 +
            1;
      /*			fprintf(stderr, "%d ", i);
      for(j=0;j<8;j++)
              fprintf(stderr, "%d ", connect[i*8+j]);
      fprintf(stderr, "\n");
       */
    }
    for (i = 0; i < mesh->n_elem; i++) mesh->elem_type[i] = 361;
    for (i = 0; i < mesh->n_elem + 1; i++) mesh->elem_node_index[i] = 8 * i;
    HECMW_free(mesh->elem_node_item);
    mesh->elem_node_item = connect;
  }
  tn_component = 0;
  for (i = 0; i < data->nn_component; i++) tn_component += data->nn_dof[i];
  value =
      (double *)HECMW_calloc(tn_component * mesh->n_node * 2, sizeof(double));
  if (value == NULL) HECMW_vis_memory_exit("value");
  for (i = 0; i < mesh->n_node; i++) {
    for (j = 0; j < tn_component; j++)
      value[i * 2 * tn_component + j] =
          data->node_val_item[i * tn_component + j];
    for (j = 0; j < tn_component; j++)
      value[(i * 2 + 1) * tn_component + j] =
          data->node_val_item[i * tn_component + j];
  }
  HECMW_free(data->node_val_item);
  data->node_val_item = value;

  fprintf(stderr, "It is ok to transform shell to solid\n");
  mesh->n_node *= 2;
  mesh->nn_internal *= 2;

  return;
}
