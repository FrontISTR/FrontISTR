/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include "hecmw_visualizer.h"

#include <stdio.h>
#include <stdlib.h>
#include "hecmw_vis_mem_util.h"
#include "hecmw_vis_read_control.h"
#include "hecmw_vis_surface_main.h"
#include "hecmw_vis_pvr_main.h"
#include "hecmw_malloc.h"

PSF_link *psf;
PVR_link *pvr;

int HECMW_visualize_init(void) {
  return HECMW_visualize_init_by_comm(HECMW_comm_get_comm());
}

int HECMW_visualize_init_by_comm(HECMW_Comm VIS_COMM) {
  FILE *contfp;
  int pesize, mynode;
  char *contfile, buf[HECMW_FILENAME_LEN];

  HECMW_Comm_size(VIS_COMM, &pesize);
  HECMW_Comm_rank(VIS_COMM, &mynode);

  if ((contfp = fopen("hecmw_vis.ini", "r")) == NULL) {
    contfile = HECMW_ctrl_get_control_file("vis_ctrl");
    if ((contfp = fopen(contfile, "r")) == NULL)
      HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0011: Cannot open control file");
  }

  psf = (PSF_link *)HECMW_malloc(sizeof(PSF_link));
  if (psf == NULL) HECMW_vis_memory_exit("psf");
  psf->next_psf   = NULL;
  psf->num_of_psf = 0;
  pvr             = (PVR_link *)HECMW_malloc(sizeof(PVR_link));
  if (pvr == NULL) HECMW_vis_memory_exit("pvr");
  pvr->next_pvr   = NULL;
  pvr->num_of_pvr = 0;

  HECMW_vis_read_control(contfp, pesize, mynode, psf, pvr);
  fclose(contfp);

  return 0;
}

int HECMW_visualize(struct hecmwST_local_mesh *mesh,
                    struct hecmwST_result_data *data, int timestep ) {
  int ii;
  char *outfile, *buf1, outfile1[HECMW_FILENAME_LEN];
  char body[HECMW_FILENAME_LEN];
  PSF_link *tp1;
  PVR_link *tv1;
  int visual_id, init_flag, fg_text;
  Parameter_rendering *sr;
  struct surface_module *sf;
  Parameter_vr *vr;
  int stat_para_sf[NUM_CONTROL_PSF], stat_para_vr[NUM_CONTROL_PVR];
  HECMW_Comm VIS_COMM;
  int pesize, mynode;

  HECMW_Comm_dup(mesh->HECMW_COMM, &VIS_COMM);
  HECMW_Comm_size(VIS_COMM, &pesize);
  HECMW_Comm_rank(VIS_COMM, &mynode);

  outfile = HECMW_ctrl_get_result_fileheader("vis_out", timestep, &fg_text);
  buf1 = HECMW_ctrl_get_result_filebody("vis_out");

  if (HECMW_ctrl_is_subdir()) {
    if (HECMW_ctrl_make_subdir(outfile)) {
      HECMW_vis_print_exit(
          "ERROR: HEC-MW-VIS-E0009: Cannot open output directory");
    }
  }
  if (psf->num_of_psf > 0) {
    init_flag = 1;
    tp1       = psf->next_psf;
    for (visual_id = 0; visual_id < psf->num_of_psf; visual_id++) {
      if (mynode == 0)
        fprintf(stderr, " Start visualize PSF %d at timestep %d\n",
                visual_id + 1, timestep);
      sf = tp1->sf;
      sr = tp1->sr;
      for (ii            = 0; ii < NUM_CONTROL_PSF; ii++)
        stat_para_sf[ii] = tp1->stat_para[ii];
      tp1                = tp1->next_psf;
      if (psf->num_of_psf > 1) {
        if (timestep >= 1000) {
          sprintf(outfile1, "%s_psf%d.%d", outfile, visual_id + 1, timestep);
          sprintf(body, "%s_psf%d.%d", buf1, visual_id + 1, timestep);
        } else if ((timestep >= 100) && (timestep <= 999)) {
          sprintf(outfile1, "%s_psf%d.0%d", outfile, visual_id + 1, timestep);
          sprintf(body, "%s_psf%d.0%d", buf1, visual_id + 1, timestep);
        } else if ((timestep >= 10) && (timestep <= 99)) {
          sprintf(outfile1, "%s_psf%d.00%d", outfile, visual_id + 1, timestep);
          sprintf(body, "%s_psf%d.00%d", buf1, visual_id + 1, timestep);
        } else if (timestep <= 9) {
          sprintf(outfile1, "%s_psf%d.000%d", outfile, visual_id + 1, timestep);
          sprintf(body, "%s_psf%d.000%d", buf1, visual_id + 1, timestep);
        }
      } else {
        if (timestep >= 1000) {
          sprintf(outfile1, "%s_psf.%d", outfile, timestep);
          sprintf(body, "%s_psf.%d", buf1, timestep);
        } else if ((timestep >= 100) && (timestep <= 999)) {
          sprintf(outfile1, "%s_psf.0%d", outfile, timestep);
          sprintf(body, "%s_psf.0%d", buf1, timestep);
        } else if ((timestep >= 10) && (timestep <= 99)) {
          sprintf(outfile1, "%s_psf.00%d", outfile, timestep);
          sprintf(body, "%s_psf.00%d", buf1, timestep);
        } else if (timestep <= 9) {
          sprintf(outfile1, "%s_psf.000%d", outfile, timestep);
          sprintf(body, "%s_psf.000%d", buf1, timestep);
        }
      }
      HECMW_vis_psf_rendering(mesh, data, &timestep, sf, sr, stat_para_sf,
                              outfile1, body, VIS_COMM);
      init_flag = 0;
    }
  }
  if (pvr->num_of_pvr > 0) {
    tv1       = pvr->next_pvr;
    init_flag = 1;
    for (visual_id = 0; visual_id < pvr->num_of_pvr; visual_id++) {
      if (mynode == 0)
        fprintf(stderr, " Start visualize PVR %d at timestep %d\n",
                visual_id + 1, timestep);
      vr = tv1->vr;
      for (ii            = 0; ii < NUM_CONTROL_PVR; ii++)
        stat_para_vr[ii] = tv1->stat_para[ii];
      tv1                = tv1->next_pvr;
      if (pvr->num_of_pvr > 1) {
        if (timestep >= 1000)
          sprintf(outfile1, "%s_pvr%d.%d", outfile, visual_id + 1, timestep);
        else if ((timestep >= 100) && (timestep <= 999))
          sprintf(outfile1, "%s_pvr%d.0%d", outfile, visual_id + 1, timestep);
        else if ((timestep >= 10) && (timestep <= 99))
          sprintf(outfile1, "%s_pvr%d.00%d", outfile, visual_id + 1,
                  timestep);
        else if (timestep <= 9)
          sprintf(outfile1, "%s_pvr%d.000%d", outfile, visual_id + 1,
                  timestep);
      } else {
        if (timestep >= 1000)
          sprintf(outfile1, "%s_pvr.%d", outfile, timestep);
        else if ((timestep >= 100) && (timestep <= 999))
          sprintf(outfile1, "%s_pvr.0%d", outfile, timestep);
        else if ((timestep >= 10) && (timestep <= 99))
          sprintf(outfile1, "%s_pvr.00%d", outfile, timestep);
        else if (timestep <= 9)
          sprintf(outfile1, "%s_pvr.000%d", outfile, timestep);
      }
      HECMW_vis_pvr_rendering(mesh, data, &timestep, &init_flag,
                              pvr->num_of_pvr, vr, stat_para_vr, outfile1,
                              VIS_COMM);
      init_flag = 0;
    }
  }
  HECMW_free(buf1);
  HECMW_free(outfile);
  HECMW_Comm_free(&VIS_COMM);

  return 0;
}

int HECMW_visualize_finalize(void) {
  PSF_link *tp1, *tp2;
  PVR_link *tv1, *tv2;
  int i;

  if (psf->num_of_psf > 0) {
    tp1 = psf->next_psf;
    for (i = 0; i < psf->num_of_psf; i++) {
      tp2 = tp1;
      tp1 = tp1->next_psf;
      HECMW_free(tp2->sf);
      if (tp2->sr->light_point) HECMW_free(tp2->sr->light_point);
      HECMW_free(tp2->sr);
      HECMW_free(tp2);
    }
  }
  HECMW_free(psf);
  if (pvr->num_of_pvr > 0) {
    tv1 = pvr->next_pvr;
    for (i = 0; i < pvr->num_of_pvr; i++) {
      tv2 = tv1;
      tv1 = tv1->next_pvr;
      HECMW_free(tv2->vr);
      HECMW_free(tv2);
    }
  }
  HECMW_free(pvr);

  return 0;
}
