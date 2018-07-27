/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include "hecmw_fstr_output_femap.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "hecmw_malloc.h"
#include "hecmw_etype.h"
#include "hecmw_vis_mem_util.h"
#include "hecmw_vis_comm_util.h"
#include "hecmw_vis_endian.h"
#include "hecmw_vis_combine.h"

static void count_data_components(struct hecmwST_result_data *data,
                                  int *tn_component, int *te_component) {
  int i;
  *tn_component = 0;
  for (i = 0; i < data->nn_component; i++) *tn_component += data->nn_dof[i];

  *te_component = 0;
  for (i = 0; i < data->ne_component; i++) *te_component += data->ne_dof[i];
}

static void femap_write_elem(FILE *outfp, int mynode, int n_elem, int n_node,
                             int *global_elem_ID, int *elem_ID, int *elem_type,
                             int *section_ID, int *elem_node_index,
                             int *elem_node_item, int *global_node_ID,
                             int *sect_opt) {
  int i, j, m;
  int ielm, icol, isid, isop, istyp, ietyp, itopo;
  int nna[10], nnb[10], nn[20];

  for (i = 0; i < n_elem; i++) {
    if (elem_ID[i * 2 + 1] != mynode) continue;

    for (m = 0; m < 20; m++) {
      nn[m] = 0;
    }
    for (m = 0; m < 10; m++) {
      nna[m] = 0;
      nnb[m] = 0;
    }

    ielm = global_elem_ID[i];
    icol = 124;
    isid = section_ID[i];
    isop = sect_opt[isid];
    for (j = 0; j < elem_node_index[i + 1] - elem_node_index[i]; j++) {
      nn[j] = global_node_ID[elem_node_item[elem_node_index[i] + j] - 1];
    }
    ietyp = elem_type[i];
    if (ietyp == 231) {
      istyp                = 25;
      if (isop == 1) istyp = 19;
      if (isop == 2) istyp = 35;
      itopo                = 2;
      nna[0]               = nn[0];
      nna[1]               = nn[1];
      nna[2]               = nn[2];
    } else if (ietyp == 731) {
      istyp  = 17;
      itopo  = 2;
      nna[0] = nn[0];
      nna[1] = nn[1];
      nna[2] = nn[2];
    } else if (ietyp == 232) {
      istyp                = 26;
      if (isop == 1) istyp = 20;
      if (isop == 2) istyp = 36;
      itopo                = 3;
      nna[0]               = nn[0];
      nna[1]               = nn[1];
      nna[2]               = nn[2];
      nna[4]               = nn[5];
      nna[5]               = nn[3];
      nna[6]               = nn[4];
    } else if (ietyp == 732) {
      istyp  = 18;
      itopo  = 3;
      nna[0] = nn[0];
      nna[1] = nn[1];
      nna[2] = nn[2];
      nna[4] = nn[5];
      nna[5] = nn[3];
      nna[6] = nn[4];
    } else if (ietyp == 241) {
      istyp                = 25;
      if (isop == 1) istyp = 19;
      if (isop == 2) istyp = 35;
      itopo                = 4;
      nna[0]               = nn[0];
      nna[1]               = nn[1];
      nna[2]               = nn[2];
      nna[3]               = nn[3];
    } else if (ietyp == 741) {
      istyp  = 17;
      itopo  = 4;
      nna[0] = nn[0];
      nna[1] = nn[1];
      nna[2] = nn[2];
      nna[3] = nn[3];
    } else if (ietyp == 242) {
      istyp                = 26;
      if (isop == 1) istyp = 20;
      if (isop == 2) istyp = 36;
      itopo                = 5;
      for (m = 0; m < 8; m++) nna[m] = nn[m];
    } else if (ietyp == 742) {
      istyp = 18;
      itopo = 5;
      for (m = 0; m < 8; m++) nna[m] = nn[m];
    } else if (ietyp == 743) {
      istyp = 18;
      itopo = 5;
      for (m = 0; m < 8; m++) nna[m] = nn[m];
    } else if (ietyp == 341 || ietyp == 3414) {
      istyp  = 25;
      itopo  = 6;
      nna[0] = nn[0];
      nna[1] = nn[1];
      nna[2] = nn[2];
      nna[4] = nn[3];
    } else if (ietyp == 351) {
      istyp  = 25;
      itopo  = 7;
      nna[0] = nn[0];
      nna[1] = nn[1];
      nna[2] = nn[2];
      nna[4] = nn[3];
      nna[5] = nn[4];
      nna[6] = nn[5];
    } else if (ietyp == 361) {
      istyp = 25;
      itopo = 8;
      for (m = 0; m < 8; m++) nna[m] = nn[m];
    } else if (ietyp == 342 || ietyp == 3422) {
      istyp  = 26;
      itopo  = 10;
      nna[0] = nn[0];
      nna[1] = nn[1];
      nna[2] = nn[2];
      nna[4] = nn[3];
      nna[8] = nn[6];
      nna[9] = nn[4];
      nnb[0] = nn[5];
      nnb[2] = nn[7];
      nnb[3] = nn[8];
      nnb[4] = nn[9];
    } else if (ietyp == 352) {
      istyp  = 26;
      itopo  = 11;
      nna[0] = nn[0];
      nna[1] = nn[1];
      nna[2] = nn[2];
      nna[4] = nn[3];
      nna[5] = nn[4];
      nna[6] = nn[5];
      nna[8] = nn[8];
      nna[9] = nn[6];
      nnb[0] = nn[7];
      nnb[2] = nn[12];
      nnb[3] = nn[13];
      nnb[4] = nn[14];
      nnb[6] = nn[11];
      nnb[7] = nn[9];
      nnb[8] = nn[10];
    } else if (ietyp == 362) {
      istyp = 26;
      itopo = 12;
      for (m = 0; m < 10; m++) nna[m] = nn[m];
      nnb[0]                          = nn[10];
      nnb[1]                          = nn[11];
      nnb[2]                          = nn[16];
      nnb[3]                          = nn[17];
      nnb[4]                          = nn[18];
      nnb[5]                          = nn[19];
      nnb[6]                          = nn[12];
      nnb[7]                          = nn[13];
      nnb[8]                          = nn[14];
      nnb[9]                          = nn[15];
    }

    fprintf(outfp, "%8d,%8d,%8d,%8d,%8d,1,0,0,0,0,0,0,0,\n", ielm, icol, isid,
            istyp, itopo);
    for (m = 0; m < 10; m++) fprintf(outfp, "%8d,", nna[m]);
    fprintf(outfp, "\n");
    for (m = 0; m < 10; m++) fprintf(outfp, "%8d,", nnb[m]);
    fprintf(outfp, "\n");
    fprintf(outfp, "0,0,0,\n");
    fprintf(outfp, "0,0,0,\n");
    fprintf(outfp, "0,0,0,\n");
    fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,\n");
  }
}

void HECMW_fstr_output_femap(struct hecmwST_local_mesh *mesh,
                             struct hecmwST_result_data *data, char *outfile,
                             HECMW_Comm VIS_COMM) {
  int i, j, k, m;
  int mynode, pesize;
  HECMW_Status stat;
  double tmp;
  int tmp_int, tn_component, tmp_int2, te_component;
  double *tmp_recv_d, *tmp_send_d;
  int *tmp_recv_i, *tmp_elem_ID, *tmp_elem_type, *tmp_elem_global_ID,
      *tmp_elem_node_index, *tmp_elem_node_item, *tmp_section_ID,
      *tmp_node_global_ID, *tmp_send_i;
  int disp_comp, stress_comp, disp_base, stress_base, name_len;
  int temp_comp, temp_base;
  FILE *outfp;

  HECMW_Comm_rank(VIS_COMM, &mynode);
  HECMW_Comm_size(VIS_COMM, &pesize);

  /* open file */
  if (mynode == 0) {
    strcat(outfile, ".neu");
    outfp = fopen(outfile, "w");

    if (!outfp)
      HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0009: Cannot open output file");

    put_neutral_head(outfp);
    /* start writing material data */
    put_neutral_601(outfp, mesh);
    put_neutral_402(outfp, mesh);
  }

  /* output nodes */
  if (mynode != 0) {
    HECMW_Send(&mesh->nn_internal, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);
    if (mesh->nn_internal > 0) {
      HECMW_Send(mesh->global_node_ID, mesh->nn_internal, HECMW_INT, MASTER_PE,
                 0, VIS_COMM);
      HECMW_Send(mesh->node, mesh->nn_internal * 3, HECMW_DOUBLE, MASTER_PE, 0,
                 VIS_COMM);
    }
  } else {
    fprintf(outfp, "   -1\n");
    fprintf(outfp, "   403\n");

    for (i = 0; i < mesh->nn_internal; i++)
      fprintf(outfp, "%8d,0,0,1,46,0,0,0,0,0,0,%15.7e,%15.7e,%15.7e\n",
              mesh->global_node_ID[i], mesh->node[i * 3], mesh->node[i * 3 + 1],
              mesh->node[i * 3 + 2]);

    for (i = 1; i < pesize; i++) {
      HECMW_Recv(&tmp_int, 1, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM, &stat);
      if (tmp_int > 0) {
        tmp_recv_i = (int *)HECMW_calloc(tmp_int, sizeof(int));
        tmp_recv_d = (double *)HECMW_calloc(tmp_int * 3, sizeof(double));
        if ((tmp_recv_i == NULL) || (tmp_recv_d == NULL))
          HECMW_vis_memory_exit("tmp_recv");

        HECMW_Recv(tmp_recv_i, tmp_int, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM,
                   &stat);
        HECMW_Recv(tmp_recv_d, tmp_int * 3, HECMW_DOUBLE, i, HECMW_ANY_TAG,
                   VIS_COMM, &stat);

        for (j = 0; j < tmp_int; j++)
          fprintf(outfp, "%8d,0,0,1,46,0,0,0,0,0,0,%15.7e,%15.7e,%15.7e\n",
                  tmp_recv_i[j], tmp_recv_d[j * 3], tmp_recv_d[j * 3 + 1],
                  tmp_recv_d[j * 3 + 2]);

        HECMW_free(tmp_recv_i);
        HECMW_free(tmp_recv_d);
      }
    }
    fprintf(outfp, "   -1\n");
  }

  /*  output element */
  if (mynode != 0) {
    HECMW_Send(&mesh->n_elem, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);
    if (mesh->n_elem > 0) {
      HECMW_Send(&mesh->n_node, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);

      HECMW_Send(mesh->global_elem_ID, mesh->n_elem, HECMW_INT, MASTER_PE, 0,
                 VIS_COMM);
      HECMW_Send(mesh->elem_ID, mesh->n_elem * 2, HECMW_INT, MASTER_PE, 0,
                 VIS_COMM);
      HECMW_Send(mesh->elem_type, mesh->n_elem, HECMW_INT, MASTER_PE, 0,
                 VIS_COMM);
      HECMW_Send(mesh->section_ID, mesh->n_elem, HECMW_INT, MASTER_PE, 0,
                 VIS_COMM);
      HECMW_Send(mesh->elem_node_index, mesh->n_elem + 1, HECMW_INT, MASTER_PE,
                 0, VIS_COMM);
      HECMW_Send(mesh->elem_node_item, mesh->elem_node_index[mesh->n_elem],
                 HECMW_INT, MASTER_PE, 0, VIS_COMM);
      HECMW_Send(mesh->global_node_ID, mesh->n_node, HECMW_INT, MASTER_PE, 0,
                 VIS_COMM);
    }
  } else {
    fprintf(outfp, "   -1\n");
    fprintf(outfp, "   404\n");

    femap_write_elem(
        outfp, mynode, mesh->n_elem, mesh->n_node, mesh->global_elem_ID,
        mesh->elem_ID, mesh->elem_type, mesh->section_ID, mesh->elem_node_index,
        mesh->elem_node_item, mesh->global_node_ID, mesh->section->sect_opt);

    for (j = 1; j < pesize; j++) {
      HECMW_Recv(&tmp_int, 1, HECMW_INT, j, HECMW_ANY_TAG, VIS_COMM, &stat);

      if (tmp_int > 0) {
        HECMW_Recv(&tmp_int2, 1, HECMW_INT, j, HECMW_ANY_TAG, VIS_COMM, &stat);

        tmp_elem_global_ID  = (int *)HECMW_calloc(tmp_int, sizeof(int));
        tmp_elem_ID         = (int *)HECMW_calloc(tmp_int * 2, sizeof(int));
        tmp_elem_type       = (int *)HECMW_calloc(tmp_int, sizeof(int));
        tmp_section_ID      = (int *)HECMW_calloc(tmp_int, sizeof(int));
        tmp_elem_node_index = (int *)HECMW_calloc(tmp_int + 1, sizeof(int));
        if ((tmp_elem_global_ID == NULL) || (tmp_elem_ID == NULL) ||
            (tmp_elem_type == NULL) || (tmp_section_ID == NULL))
          HECMW_vis_memory_exit("tmp recv");

        HECMW_Recv(tmp_elem_global_ID, tmp_int, HECMW_INT, j, HECMW_ANY_TAG,
                   VIS_COMM, &stat);
        HECMW_Recv(tmp_elem_ID, tmp_int * 2, HECMW_INT, j, HECMW_ANY_TAG,
                   VIS_COMM, &stat);
        HECMW_Recv(tmp_elem_type, tmp_int, HECMW_INT, j, HECMW_ANY_TAG,
                   VIS_COMM, &stat);
        HECMW_Recv(tmp_section_ID, tmp_int, HECMW_INT, j, HECMW_ANY_TAG,
                   VIS_COMM, &stat);
        HECMW_Recv(tmp_elem_node_index, tmp_int + 1, HECMW_INT, j,
                   HECMW_ANY_TAG, VIS_COMM, &stat);

        tmp_elem_node_item =
            (int *)HECMW_calloc(tmp_elem_node_index[tmp_int], sizeof(int));
        tmp_node_global_ID = (int *)HECMW_calloc(tmp_int2, sizeof(int));

        HECMW_Recv(tmp_elem_node_item, tmp_elem_node_index[tmp_int], HECMW_INT,
                   j, HECMW_ANY_TAG, VIS_COMM, &stat);
        HECMW_Recv(tmp_node_global_ID, tmp_int2, HECMW_INT, j, HECMW_ANY_TAG,
                   VIS_COMM, &stat);

        femap_write_elem(outfp, j, tmp_int, tmp_int2, tmp_elem_global_ID,
                         tmp_elem_ID, tmp_elem_type, tmp_section_ID,
                         tmp_elem_node_index, tmp_elem_node_item,
                         tmp_node_global_ID, mesh->section->sect_opt);

        HECMW_free(tmp_elem_global_ID);
        HECMW_free(tmp_elem_ID);
        HECMW_free(tmp_elem_type);
        HECMW_free(tmp_section_ID);
        HECMW_free(tmp_elem_node_index);
        HECMW_free(tmp_elem_node_item);
        HECMW_free(tmp_node_global_ID);
      }
    }
    fprintf(outfp, "   -1\n");
  }

  /* output data */

  if (mynode == 0) {
    put_neutral_middle(outfp);
  }
  count_data_components(data, &tn_component, &te_component);

  /* output node-based data */

  /* temperature */

  temp_comp = -1;
  for (j = 0; j < data->nn_component; j++) {
    name_len = strlen(data->node_label[j]);
    if (strncmp("TEMPERATURE", data->node_label[j], name_len) == 0) {
      temp_comp = j;
    }
  }
  if (temp_comp >= 0) {
    temp_base = 0;
    for (i = 0; i < temp_comp; i++) temp_base += data->nn_dof[i];
    if (mynode != 0) {
      HECMW_Send(&mesh->nn_internal, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);
      if (mesh->nn_internal > 0) {
        HECMW_Send(mesh->global_node_ID, mesh->nn_internal, HECMW_INT,
                   MASTER_PE, 0, VIS_COMM);

        tmp_send_d = (double *)HECMW_calloc(mesh->nn_internal, sizeof(double));

        for (i          = 0; i < mesh->nn_internal; i++)
          tmp_send_d[i] = data->node_val_item[i * tn_component + temp_base];

        HECMW_Send(tmp_send_d, mesh->nn_internal, HECMW_DOUBLE, MASTER_PE, 0,
                   VIS_COMM);

        HECMW_free(tmp_send_d);
      }
    } else {
      fprintf(outfp, "   -1\n");
      fprintf(outfp, "  451\n");
      fprintf(outfp, "1,1,1,\n");
      fprintf(outfp, "Temperature\n");
      fprintf(outfp, "1.0, 0.0, 0.0\n");
      fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
      fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
      fprintf(outfp, "0,0,6,7,\n");
      fprintf(outfp, "0,0,1,\n");

      for (i = 0; i < mesh->nn_internal; i++) {
        tmp = data->node_val_item[i * tn_component + temp_base];
        fprintf(outfp, "%8d,%15.7e,\n", mesh->global_node_ID[i], tmp);
      }
      for (i = 1; i < pesize; i++) {
        HECMW_Recv(&tmp_int, 1, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM, &stat);
        if (tmp_int > 0) {
          tmp_recv_i = (int *)HECMW_calloc(tmp_int, sizeof(int));
          tmp_recv_d = (double *)HECMW_calloc(tmp_int, sizeof(double));
          if ((tmp_recv_i == NULL) || (tmp_recv_d == NULL))
            HECMW_vis_memory_exit("tmp_recv");
          HECMW_Recv(tmp_recv_i, tmp_int, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM,
                     &stat);
          HECMW_Recv(tmp_recv_d, tmp_int, HECMW_DOUBLE, i, HECMW_ANY_TAG,
                     VIS_COMM, &stat);
          for (j = 0; j < tmp_int; j++)
            fprintf(outfp, "%8d,%15.7e,\n", tmp_recv_i[j], tmp_recv_d[j]);
          HECMW_free(tmp_recv_i);
          HECMW_free(tmp_recv_d);
        }
      }
      fprintf(outfp, "-1,0.\n");
    }
    return;
  }

  /* displacement and stress */

  disp_comp = 0;
  for (j = 0; j < data->nn_component; j++) {
    name_len = strlen(data->node_label[j]);
    if (strncmp("DISPLACEMENT", data->node_label[j], name_len) == 0) {
      disp_comp = j;
    }
  }
  stress_comp = 0;
  for (j = 0; j < data->nn_component; j++) {
    name_len = strlen(data->node_label[j]);
    if (strncmp("NodalSTRESS", data->node_label[j], name_len) == 0 ||
        strncmp("NodalSTRESSplus", data->node_label[j], name_len) == 0) {
      stress_comp = j;
    }
  }

  disp_base   = 0;
  stress_base = 0;
  for (i = 0; i < disp_comp; i++) disp_base += data->nn_dof[i];
  for (i = 0; i < stress_comp; i++) stress_base += data->nn_dof[i];
  if (mynode != 0) {
    HECMW_Send(&mesh->nn_internal, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);
    if (mesh->nn_internal > 0) {
      HECMW_Send(mesh->global_node_ID, mesh->nn_internal, HECMW_INT, MASTER_PE,
                 0, VIS_COMM);

      tmp_send_d = (double *)HECMW_calloc(mesh->nn_internal, sizeof(double));

      for (i = 0; i < mesh->nn_internal; i++) {
        tmp_send_d[i] =
            sqrt(data->node_val_item[i * tn_component + disp_base] *
                     data->node_val_item[i * tn_component + disp_base] +
                 data->node_val_item[i * tn_component + disp_base + 1] *
                     data->node_val_item[i * tn_component + disp_base + 1] +
                 data->node_val_item[i * tn_component + disp_base + 2] *
                     data->node_val_item[i * tn_component + disp_base + 2]);
      }

      HECMW_Send(tmp_send_d, mesh->nn_internal, HECMW_DOUBLE, MASTER_PE, 0,
                 VIS_COMM);
      HECMW_Send(&mesh->nn_internal, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);
      HECMW_Send(mesh->global_node_ID, mesh->nn_internal, HECMW_INT, MASTER_PE,
                 0, VIS_COMM);

      for (i          = 0; i < mesh->nn_internal; i++)
        tmp_send_d[i] = data->node_val_item[i * tn_component + disp_base];

      HECMW_Send(tmp_send_d, mesh->nn_internal, HECMW_DOUBLE, MASTER_PE, 0,
                 VIS_COMM);
      HECMW_Send(&mesh->nn_internal, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);
      HECMW_Send(mesh->global_node_ID, mesh->nn_internal, HECMW_INT, MASTER_PE,
                 0, VIS_COMM);

      for (i          = 0; i < mesh->nn_internal; i++)
        tmp_send_d[i] = data->node_val_item[i * tn_component + disp_base + 1];

      HECMW_Send(tmp_send_d, mesh->nn_internal, HECMW_DOUBLE, MASTER_PE, 0,
                 VIS_COMM);
      HECMW_Send(&mesh->nn_internal, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);
      HECMW_Send(mesh->global_node_ID, mesh->nn_internal, HECMW_INT, MASTER_PE,
                 0, VIS_COMM);

      for (i          = 0; i < mesh->nn_internal; i++)
        tmp_send_d[i] = data->node_val_item[i * tn_component + disp_base + 2];

      HECMW_Send(tmp_send_d, mesh->nn_internal, HECMW_DOUBLE, MASTER_PE, 0,
                 VIS_COMM);

      if (mesh->n_dof != 6) {
        for (j = 0; j < 7; j++) {
          HECMW_Send(&mesh->nn_internal, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);
          HECMW_Send(mesh->global_node_ID, mesh->nn_internal, HECMW_INT,
                     MASTER_PE, 0, VIS_COMM);
          for (i = 0; i < mesh->nn_internal; i++)
            tmp_send_d[i] =
                data->node_val_item[i * tn_component + stress_base + j];
          HECMW_Send(tmp_send_d, mesh->nn_internal, HECMW_DOUBLE, MASTER_PE, 0,
                     VIS_COMM);
        }
      } else {
        for (j = 0; j < 14; j++) {
          HECMW_Send(&mesh->nn_internal, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);
          HECMW_Send(mesh->global_node_ID, mesh->nn_internal, HECMW_INT,
                     MASTER_PE, 0, VIS_COMM);
          for (i = 0; i < mesh->nn_internal; i++)
            tmp_send_d[i] =
                data->node_val_item[i * tn_component + stress_base + j];
          HECMW_Send(tmp_send_d, mesh->nn_internal, HECMW_DOUBLE, MASTER_PE, 0,
                     VIS_COMM);
        }
      }
      HECMW_free(tmp_send_d);
    }
  } else {
    fprintf(outfp, "   -1\n");
    fprintf(outfp, "  451\n");
    for (k = 0; k < 4; k++) {
      fprintf(outfp, "1,%1d,1,\n", k + 1);
      if (k == 0)
        fprintf(outfp, "Total Translation\n");
      else
        fprintf(outfp, "T%1d Translation\n", k);
      fprintf(outfp, "1.0, 0.0, 0.0\n");
      if (k == 0)
        fprintf(outfp, "2,3,4,0,0,0,0,0,0,0,\n");
      else if (k == 1)
        fprintf(outfp, "2,0,0,0,0,0,0,0,0,0,\n");
      else if (k == 2)
        fprintf(outfp, "0,3,0,0,0,0,0,0,0,0,\n");
      else if (k == 3)
        fprintf(outfp, "0,0,4,0,0,0,0,0,0,0,\n");
      fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
      fprintf(outfp, "0,0,1,7,\n");
      if (k == 0)
        fprintf(outfp, "1,1,1,\n");
      else
        fprintf(outfp, "0,1,1,\n");

      for (i = 0; i < mesh->nn_internal; i++) {
        if (k == 0)
          tmp = sqrt(data->node_val_item[i * tn_component + disp_base] *
                         data->node_val_item[i * tn_component + disp_base] +
                     data->node_val_item[i * tn_component + disp_base + 1] *
                         data->node_val_item[i * tn_component + disp_base + 1] +
                     data->node_val_item[i * tn_component + disp_base + 2] *
                         data->node_val_item[i * tn_component + disp_base + 2]);
        else
          tmp = data->node_val_item[i * tn_component + disp_base + k - 1];
        fprintf(outfp, "%8d,%15.7e,\n", mesh->global_node_ID[i], tmp);
      }
      for (i = 1; i < pesize; i++) {
        HECMW_Recv(&tmp_int, 1, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM, &stat);
        if (tmp_int > 0) {
          tmp_recv_i = (int *)HECMW_calloc(tmp_int, sizeof(int));
          tmp_recv_d = (double *)HECMW_calloc(tmp_int, sizeof(double));
          if ((tmp_recv_i == NULL) || (tmp_recv_d == NULL))
            HECMW_vis_memory_exit("tmp_recv");
          HECMW_Recv(tmp_recv_i, tmp_int, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM,
                     &stat);
          HECMW_Recv(tmp_recv_d, tmp_int, HECMW_DOUBLE, i, HECMW_ANY_TAG,
                     VIS_COMM, &stat);
          for (j = 0; j < tmp_int; j++)
            fprintf(outfp, "%8d,%15.7e,\n", tmp_recv_i[j], tmp_recv_d[j]);
          HECMW_free(tmp_recv_i);
          HECMW_free(tmp_recv_d);
        }
      }
      fprintf(outfp, "-1,0.\n");
    }
    if (mesh->n_dof != 6) {
      char *nodal_stress_titles_solid[] = {
          "Solid X Normal Stress", "Solid Y Normal Stress",
          "Solid Z Normal Stress", "Solid XY Shear Stress",
          "Solid YZ Shear Stress", "Solid XZ Shear Stress",
          "Solid Von Mises Stress"};
      for (k = 0; k < 7; k++) {
        fprintf(outfp, "1,6001%1d,1,\n", k + 1);
        fprintf(outfp, "%s\n", nodal_stress_titles_solid[k]);
        fprintf(outfp, "1.0, 0.0, 0.0\n");
        fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
        fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
        fprintf(outfp, "0,0,4,7,\n");
        fprintf(outfp, "0,0,1,\n");

        for (i = 0; i < mesh->nn_internal; i++) {
          tmp = data->node_val_item[i * tn_component + stress_base + k];
          fprintf(outfp, "%8d,%15.7e,\n", mesh->global_node_ID[i], tmp);
        }
        for (i = 1; i < pesize; i++) {
          HECMW_Recv(&tmp_int, 1, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM, &stat);
          if (tmp_int > 0) {
            tmp_recv_i = (int *)HECMW_calloc(tmp_int, sizeof(int));
            tmp_recv_d = (double *)HECMW_calloc(tmp_int, sizeof(double));
            if ((tmp_recv_i == NULL) || (tmp_recv_d == NULL))
              HECMW_vis_memory_exit("tmp_recv");
            HECMW_Recv(tmp_recv_i, tmp_int, HECMW_INT, i, HECMW_ANY_TAG,
                       VIS_COMM, &stat);
            HECMW_Recv(tmp_recv_d, tmp_int, HECMW_DOUBLE, i, HECMW_ANY_TAG,
                       VIS_COMM, &stat);
            for (j = 0; j < tmp_int; j++)
              fprintf(outfp, "%8d,%15.7e,\n", tmp_recv_i[j], tmp_recv_d[j]);
            HECMW_free(tmp_recv_i);
            HECMW_free(tmp_recv_d);
          }
        }
        fprintf(outfp, "-1,0.\n");
      }
    } else if (mesh->n_dof == 6) {
      char *nodal_stress_ids_shell[] = {
          "1,70011,1,", "1,70012,1,", "1,70013,1,", "1,70014,1,", "1,70015,1,",
          "1,70016,1,", "1,71011,1,", "1,71012,1,", "1,71013,1,", "1,71014,1,",
          "1,71015,1,", "1,71016,1,", "1,70017,1,", "1,71017,1,"};
      char *nodal_stress_titles_shell[] = {
          "Plate Top X Normal Stress",  "Plate Top Y Normal Stress",
          "Plate Top Z Normal Stress",  "Plate Top XY Shear Stress",
          "Plate Top YZ Shear Stress",  "Plate Top XZ Shear Stress",
          "Plate Bot X Normal Stress",  "Plate Bot Y Normal Stress",
          "Plate Bot Z Normal Stress",  "Plate Bot XY Shear Stress",
          "Plate Bot YZ Shear Stress",  "Plate Bot XZ Shear Stress",
          "Plate Top Von Mises Stress", "Plate Bot Von Mises Stress"};
      for (k = 0; k < 14; k++) {
        fprintf(outfp, "%s\n", nodal_stress_ids_shell[k]);
        fprintf(outfp, "%s\n", nodal_stress_titles_shell[k]);
        fprintf(outfp, "1.0, 0.0, 0.0\n");
        fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
        fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
        fprintf(outfp, "0,0,4,7,\n");
        fprintf(outfp, "0,0,1,\n");

        for (i = 0; i < mesh->nn_internal; i++) {
          tmp = data->node_val_item[i * tn_component + stress_base + k];
          fprintf(outfp, "%8i,%15.7e,\n", mesh->global_node_ID[i], tmp);
        }
        for (i = 1; i < pesize; i++) {
          HECMW_Recv(&tmp_int, 1, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM, &stat);
          if (tmp_int > 0) {
            tmp_recv_i = (int *)HECMW_calloc(tmp_int, sizeof(int));
            tmp_recv_d = (double *)HECMW_calloc(tmp_int, sizeof(double));
            if ((tmp_recv_i == NULL) || (tmp_recv_d == NULL))
              HECMW_vis_memory_exit("tmp_recv");

            HECMW_Recv(tmp_recv_i, tmp_int, HECMW_INT, i, HECMW_ANY_TAG,
                       VIS_COMM, &stat);
            HECMW_Recv(tmp_recv_d, tmp_int, HECMW_DOUBLE, i, HECMW_ANY_TAG,
                       VIS_COMM, &stat);

            for (j = 0; j < tmp_int; j++)
              fprintf(outfp, "%8d,%15.7e,\n", tmp_recv_i[j], tmp_recv_d[j]);

            HECMW_free(tmp_recv_i);
            HECMW_free(tmp_recv_d);
          }
        }
        fprintf(outfp, "-1,0.\n");
      }
    }
  }

  /* starting output element-based data */

  if (data->ne_component > 0) {
    stress_comp = 0;
    for (j = 0; j < data->ne_component; j++) {
      name_len = strlen(data->elem_label[j]);
      if (strncmp("ElementalSTRESS", data->elem_label[j], name_len) == 0) {
        stress_comp = j;
      }
    }

    stress_base = 0;
    for (i = 0; i < stress_comp; i++) stress_base += data->ne_dof[i];

    if (mynode != 0) {
      tmp_send_i = (int *)HECMW_calloc(mesh->nn_internal, sizeof(int));
      tmp_send_d = (double *)HECMW_calloc(mesh->ne_internal, sizeof(double));

      for (j = 0; j < 7; j++) {
        HECMW_Send(&mesh->ne_internal, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);

        for (i = 0; i < mesh->ne_internal; i++) {
          tmp_int       = mesh->elem_internal_list[i];
          tmp_send_i[i] = mesh->global_elem_ID[tmp_int];
          tmp_send_d[i] =
              data->elem_val_item[tmp_int * te_component + stress_base + j];
        }

        HECMW_Send(tmp_send_i, mesh->ne_internal, HECMW_INT, MASTER_PE, 0,
                   VIS_COMM);
        HECMW_Send(tmp_send_d, mesh->ne_internal, HECMW_DOUBLE, MASTER_PE, 0,
                   VIS_COMM);
      }

    } else {
      char *elem_stress_titles_solid[] = {
          "Solid X Normal Stress", "Solid Y Normal Stress",
          "Solid Z Normal Stress", "Solid XY Shear Stress",
          "Solid YZ Shear Stress", "Solid XZ Shear Stress",
          "Solid Von Mises Stress"};

      for (k = 0; k < 7; k++) {
        fprintf(outfp, "1,7001%1d,1,\n", k + 1);
        fprintf(outfp, "%s\n", elem_stress_titles_solid[k]);
        fprintf(outfp, "1.0, 0.0, 0.0\n");
        fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
        fprintf(outfp, "0,0,0,0,0,0,0,0,0,0,\n");
        fprintf(outfp, "0,0,4,8,\n");
        fprintf(outfp, "0,0,1,\n");

        for (i = 0; i < mesh->ne_internal; i++) {
          tmp_int = mesh->elem_internal_list[i];
          tmp = data->elem_val_item[tmp_int * te_component + stress_base + k];
          fprintf(outfp, "%8d,%15.7e,\n", mesh->global_elem_ID[tmp_int], tmp);
        }
        for (i = 1; i < pesize; i++) {
          HECMW_Recv(&tmp_int, 1, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM, &stat);
          if (tmp_int > 0) {
            tmp_recv_i = (int *)HECMW_calloc(tmp_int, sizeof(int));
            tmp_recv_d = (double *)HECMW_calloc(tmp_int, sizeof(double));
            if ((tmp_recv_i == NULL) || (tmp_recv_d == NULL))
              HECMW_vis_memory_exit("tmp_recv");

            HECMW_Recv(tmp_recv_i, tmp_int, HECMW_INT, i, HECMW_ANY_TAG,
                       VIS_COMM, &stat);
            HECMW_Recv(tmp_recv_d, tmp_int, HECMW_DOUBLE, i, HECMW_ANY_TAG,
                       VIS_COMM, &stat);

            for (j = 0; j < tmp_int; j++)
              fprintf(outfp, "%8d,%15.7e,\n", tmp_recv_i[j], tmp_recv_d[j]);

            HECMW_free(tmp_recv_i);
            HECMW_free(tmp_recv_d);
          }
        }
        fprintf(outfp, "-1,0.\n");
      }
    }

  } /* end of if ne>0) */

  if (mynode == 0) {
    fprintf(outfp, "   -1\n");
    fprintf(outfp, "   -1\n");

    fclose(outfp);
  }

  return;
}

static void avs_write_header(FILE *outfp, int total_n_node, int total_n_elem,
                             int tn_component, int te_component,
                             int flag_oldUCD) {
  if (flag_oldUCD) {
    fprintf(outfp, "%d %d %d %d 0\n", total_n_node, total_n_elem, tn_component,
            te_component);
  } else {
    fprintf(outfp, "1\n");
    fprintf(outfp, "data\n");
    fprintf(outfp, "step1\n");
    fprintf(outfp, "%d %d\n", total_n_node, total_n_elem);
  }
}

static void avs_write_node_coord(FILE *outfp, int n_node, int *global_node_ID,
                                 double *node, int flag_global_ID,
                                 int nid_offset) {
  int i;

  HECMW_assert(!flag_global_ID || global_node_ID);

  for (i = 0; i < n_node; i++) {
    fprintf(outfp, "%8d %15.7e %15.7e %15.7e\n",
            flag_global_ID ? global_node_ID[i] : i + 1 + nid_offset,
            node[i * 3], node[i * 3 + 1], node[i * 3 + 2]);
  }
}

static const int *avs_elem_node_order(int elem_type) {
  static const int i232[] = {0, 1, 2, 5, 3, 4};
  static const int i342[] = {0, 1, 3, 2, 6, 7, 5, 8, 9, 4};
  static const int i352[] = {0, 2, 1, 3, 5, 4, 7, 6, 8, 10, 9, 11, 12, 14, 13};
  static const int i362[] = {0, 3, 2,  1,  4,  7,  6,  5,  11, 10,
                             9, 8, 15, 14, 13, 12, 16, 19, 18, 17};
  static const int idefault[] = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
  const int *ii;

  if (elem_type == 232 || elem_type == 732)
    ii = i232;
  else if (elem_type == 341 || elem_type == 342 || elem_type == 3414)
    ii = i342;
  else if (elem_type == 351 || elem_type == 352)
    ii = i352;
  else if (elem_type == 361 || elem_type == 362)
    ii = i362;
  else
    ii = idefault;

  return ii;
}

static void avs_write_elem_conn(FILE *outfp, int mynode, int n_elem,
                                int *elem_ID, int *elem_type,
                                int *elem_node_index, int *elem_node_item,
                                int *global_elem_ID, int *global_node_ID,
                                int *node_ID, int flag_global_ID,
                                int flag_skip_external, int flag_oldUCD,
                                int *nid_offsets, int eid_offset) {
  int eid = 0;
  int i, j, etype, node_num;
  int nn[20];
  const int *ii;

  HECMW_assert((flag_global_ID && global_elem_ID && global_node_ID) ||
               (!flag_global_ID && (!nid_offsets || node_ID)));

  for (i = 0; i < n_elem; i++) {
    if (flag_skip_external && elem_ID[i * 2 + 1] != mynode) continue;

    eid++;

    etype = elem_type[i];
    if (flag_oldUCD && etype % 2 == 0) {
      /* old UCD format does not support second order elements */
      if (etype != 3414 && etype != 3614) {
        etype--;
      }
    }
    node_num = HECMW_get_max_node(etype);
    if (HECMW_is_etype_33struct(etype)) node_num /= 2;

    for (j = 0; j < node_num; j++) {
      int nid = elem_node_item[elem_node_index[i] + j] - 1;
      if (flag_global_ID) {
        nn[j] = global_node_ID[nid];
      } else if (nid_offsets != NULL) {
        int rank        = node_ID[nid * 2 + 1];
        int nid_in_rank = node_ID[nid * 2];
        nn[j]           = nid_offsets[rank] + nid_in_rank;
      } else {
        nn[j] = nid + 1;
      }
    }

    fprintf(outfp, "%8d %d %s  ",
            flag_global_ID ? global_elem_ID[i] : eid + eid_offset,
            HECMW_get_etype_class(etype), HECMW_get_ucd_label(etype));

    ii = avs_elem_node_order(etype);
    for (j = 0; j < node_num; j++) fprintf(outfp, "%8d ", nn[ii[j]]);
    fprintf(outfp, "\n");
  }
}

static void avs_write_data_header(FILE *outfp, int tn_component,
                                  int te_component, int flag_oldUCD) {
  if (flag_oldUCD) return;
  fprintf(outfp, "%8d %8d\n", tn_component, te_component);
}

static void avs_write_node_data_header(FILE *outfp, int tn_component,
                                       struct hecmwST_result_data *data,
                                       int flag_Scalar) {
  int j, ii;

  if (flag_Scalar) {
    fprintf(outfp, "%8d", tn_component);
    for (j = 0; j < tn_component; j++) fprintf(outfp, " 1");
    fprintf(outfp, "\n");
    for (j = 0; j < data->nn_component; j++) {
      for (ii = 0; ii < data->nn_dof[j]; ii++)
        fprintf(outfp, "%s_%d, unit_unknown\n", data->node_label[j], ii + 1);
    }
  } else {
    fprintf(outfp, "%8d", data->nn_component);
    for (j = 0; j < data->nn_component; j++)
      fprintf(outfp, " %d", data->nn_dof[j]);
    fprintf(outfp, "\n");
    for (j = 0; j < data->nn_component; j++)
      fprintf(outfp, "%s, unit_unknown\n", data->node_label[j]);
  }
}

static void avs_write_node_data(FILE *outfp, int n_node, int *global_node_ID,
                                int tn_component, double *node_val_item,
                                int flag_global_ID, int nid_offset) {
  int i, k;

  HECMW_assert(!flag_global_ID || global_node_ID);

  for (i = 0; i < n_node; i++) {
    fprintf(outfp, "%8d ",
            flag_global_ID ? global_node_ID[i] : i + 1 + nid_offset);
    for (k = 0; k < tn_component; k++)
      fprintf(outfp, "%15.7e ", node_val_item[i * tn_component + k]);
    fprintf(outfp, "\n");
  }
}

int modify_element_information(const struct hecmwST_local_mesh *mesh) {
  int i, j, n, refine, max_elem;
  int *size;

  max_elem = 0;
  refine   = mesh->n_refine;
  if (refine <= 0) return 0;

  size = mesh->n_node_refine_hist;

  for (n = 0; n < refine; n++) {
    int nn = *size;
    for (i = 0; i < nn; i++) {
      if (max_elem < mesh->node_ID[2 * i]) {
        max_elem = mesh->node_ID[2 * i];
      }
    }

    for (i = 0; i < nn; i++) {
      if (mesh->node_ID[2 * i] < 0) {
        mesh->node_ID[2 * i] = max_elem - mesh->node_ID[2 * i];
      }
    }
    size++;
  }

  return 0;
}


static void
avs_write_elem_data_header(FILE *outfp, int te_component, struct hecmwST_result_data *data,
    int flag_Scalar) {
  int j, ii;

  if (flag_Scalar) {
    fprintf (outfp, "%8d", te_component);
    for (j = 0; j < te_component; j++)
      fprintf (outfp, " 1");
    fprintf (outfp, "\n");
    for (j = 0; j < data->ne_component; j++) {
      for (ii = 0; ii < data->ne_dof[j]; ii++)
        fprintf (outfp, "%s_%d, unit_unknown\n", data->elem_label[j], ii + 1);
    }
  } else {
    fprintf (outfp, "%8d", data->ne_component);
    for (j = 0; j < data->ne_component; j++)
      fprintf (outfp, " %d", data->ne_dof[j]);
    fprintf (outfp, "\n");
    for (j = 0; j < data->ne_component; j++)
      fprintf (outfp, "%s, unit_unknown\n", data->elem_label[j]);
  }
}

static void
avs_write_elem_data(FILE *outfp, int mynode,
                    int n_elem, int *elem_ID, int *global_elem_ID, int te_component, double *elem_val_item,
                    int flag_global_ID, int eid_offset) {
  int i, in, k;
  HECMW_assert(!flag_global_ID || global_elem_ID);

  in = 0;
  for (i = 0; i < n_elem; i++) {
    if (elem_ID[i * 2 + 1] != mynode)
      continue;

    fprintf (outfp, "%8d ", flag_global_ID ? global_elem_ID[i] : in+1 + eid_offset);
    for (k = 0; k < te_component; k++)
      fprintf (outfp, "%15.7e ",
          elem_val_item[i * te_component + k]);
    fprintf (outfp, "\n");
    in++;
  }
}

static void avs_output(struct hecmwST_local_mesh *mesh,
                       struct hecmwST_result_data *data, char *outfile,
                       HECMW_Comm VIS_COMM, int flag_oldUCD, int flag_global_ID,
                       int flag_Scalar) {
  int mynode, pesize;
  int total_n_node, total_n_elem;
  int tn_component, te_component;
  int *nid_offsets, *eid_offsets;
  FILE *outfp;
  HECMW_Status stat;

  int flag_skip_ext = 1;

  HECMW_Comm_rank(VIS_COMM, &mynode);
  HECMW_Comm_size(VIS_COMM, &pesize);

  /* count total number of node/elem */
  if (pesize > 1) {
    HECMW_Allreduce(&mesh->nn_internal, &total_n_node, 1, HECMW_INT, HECMW_SUM,
                    VIS_COMM);
    HECMW_Allreduce(&mesh->ne_internal, &total_n_elem, 1, HECMW_INT, HECMW_SUM,
                    VIS_COMM);
  } else {
    total_n_node = mesh->nn_internal;
    total_n_elem = mesh->n_elem;
  }

  /* count total number of data components */
  count_data_components(data, &tn_component, &te_component);

  if (mynode == 0) {
    /* allocate offset arrays */
    nid_offsets = (int *)HECMW_calloc(pesize + 1, sizeof(int));
    if (nid_offsets == NULL) HECMW_vis_memory_exit("nid_offset");
    eid_offsets = (int *)HECMW_calloc(pesize + 1, sizeof(int));
    if (eid_offsets == NULL) HECMW_vis_memory_exit("eid_offset");

    /* open file */
    outfp = fopen(outfile, "w");
    if (!outfp)
      HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0009: Cannot open output file");

    /* write header */
    avs_write_header(outfp, total_n_node, total_n_elem, tn_component, te_component,
                     flag_oldUCD);
  }

  /* write node coordinate */
  if (mynode != 0) {
    HECMW_Send(&mesh->nn_internal, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);
    if (mesh->nn_internal > 0) {
      if (flag_global_ID)
        HECMW_Send(mesh->global_node_ID, mesh->nn_internal, HECMW_INT,
                   MASTER_PE, 0, VIS_COMM);
      HECMW_Send(mesh->node, mesh->nn_internal * 3, HECMW_DOUBLE, MASTER_PE, 0,
                 VIS_COMM);
    }
  } else {
    int i;

    avs_write_node_coord(outfp, mesh->nn_internal, mesh->global_node_ID,
                         mesh->node, flag_global_ID, 0);

    nid_offsets[0] = 0;
    nid_offsets[1] = mesh->nn_internal;

    for (i = 1; i < pesize; i++) {
      int tmp_nn_internal;

      HECMW_Recv(&tmp_nn_internal, 1, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM,
                 &stat);
      nid_offsets[i + 1] = nid_offsets[i] + tmp_nn_internal;

      if (tmp_nn_internal > 0) {
        int *tmp_global_node_ID = NULL;
        double *tmp_node;

        if (flag_global_ID) {
          tmp_global_node_ID =
              (int *)HECMW_calloc(tmp_nn_internal, sizeof(int));
          if (tmp_global_node_ID == NULL)
            HECMW_vis_memory_exit("tmp recv: global_node_ID");

          HECMW_Recv(tmp_global_node_ID, tmp_nn_internal, HECMW_INT, i,
                     HECMW_ANY_TAG, VIS_COMM, &stat);
        }

        tmp_node = (double *)HECMW_calloc(tmp_nn_internal * 3, sizeof(double));
        if (tmp_node == NULL) HECMW_vis_memory_exit("tmp recv: node");

        HECMW_Recv(tmp_node, tmp_nn_internal * 3, HECMW_DOUBLE, i,
                   HECMW_ANY_TAG, VIS_COMM, &stat);

        avs_write_node_coord(outfp, tmp_nn_internal, tmp_global_node_ID,
                             tmp_node, flag_global_ID, nid_offsets[i]);

        if (flag_global_ID) HECMW_free(tmp_global_node_ID);
        HECMW_free(tmp_node);
      }
    }
  }

  /* modify element information due to refininer */
  if (mynode == 0) {
    if (modify_element_information(mesh) != 0) {
      printf("###ERROR: modify element information due to refininer \n");
    }
  }

  /* write element connectivity */
  if (mynode != 0) {
    HECMW_Send(&mesh->n_elem, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);
    HECMW_Send(&mesh->ne_internal, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);

    if (mesh->n_elem > 0) {
      HECMW_Send(&mesh->n_node, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);

      HECMW_Send(mesh->elem_ID, mesh->n_elem * 2, HECMW_INT, MASTER_PE, 0,
                 VIS_COMM);
      HECMW_Send(mesh->elem_type, mesh->n_elem, HECMW_INT, MASTER_PE, 0,
                 VIS_COMM);
      HECMW_Send(mesh->elem_node_index, mesh->n_elem + 1, HECMW_INT, MASTER_PE,
                 0, VIS_COMM);
      HECMW_Send(mesh->elem_node_item, mesh->elem_node_index[mesh->n_elem],
                 HECMW_INT, MASTER_PE, 0, VIS_COMM);

      if (flag_global_ID) {
        HECMW_Send(mesh->global_elem_ID, mesh->n_elem, HECMW_INT, MASTER_PE, 0,
                   VIS_COMM);
        HECMW_Send(mesh->global_node_ID, mesh->n_node, HECMW_INT, MASTER_PE, 0,
                   VIS_COMM);
      } else {
        HECMW_Send(mesh->node_ID, mesh->n_node * 2, HECMW_INT, MASTER_PE, 0,
                   VIS_COMM);
      }
    }
  } else {
    int j;

    avs_write_elem_conn(outfp, mynode, mesh->n_elem, mesh->elem_ID,
                        mesh->elem_type, mesh->elem_node_index,
                        mesh->elem_node_item, mesh->global_elem_ID,
                        mesh->global_node_ID, mesh->node_ID, flag_global_ID,
                        flag_skip_ext, flag_oldUCD, nid_offsets, 0);

    eid_offsets[0] = 0;
    eid_offsets[1] = mesh->ne_internal;

    for (j = 1; j < pesize; j++) {
      int tmp_n_elem;
      int tmp_ne_internal;

      HECMW_Recv(&tmp_n_elem, 1, HECMW_INT, j, HECMW_ANY_TAG, VIS_COMM, &stat);
      HECMW_Recv(&tmp_ne_internal, 1, HECMW_INT, j, HECMW_ANY_TAG, VIS_COMM,
                 &stat);
      eid_offsets[j + 1] = eid_offsets[j] + tmp_ne_internal;

      if (tmp_n_elem > 0) {
        int tmp_n_node;
        int *tmp_elem_ID;
        int *tmp_elem_type;
        int *tmp_elem_node_index;
        int *tmp_elem_node_item;
        int *tmp_global_elem_ID;
        int *tmp_global_node_ID;
        int *tmp_node_ID;

        HECMW_Recv(&tmp_n_node, 1, HECMW_INT, j, HECMW_ANY_TAG, VIS_COMM,
                   &stat);

        tmp_elem_ID         = (int *)HECMW_calloc(tmp_n_elem * 2, sizeof(int));
        tmp_elem_type       = (int *)HECMW_calloc(tmp_n_elem, sizeof(int));
        tmp_elem_node_index = (int *)HECMW_calloc(tmp_n_elem + 1, sizeof(int));
        if ((tmp_elem_ID == NULL) || (tmp_elem_type == NULL) ||
            (tmp_elem_node_index == NULL))
          HECMW_vis_memory_exit(
              "tmp recv: elem_ID, elem_type, elem_node_index");

        HECMW_Recv(tmp_elem_ID, tmp_n_elem * 2, HECMW_INT, j, HECMW_ANY_TAG,
                   VIS_COMM, &stat);
        HECMW_Recv(tmp_elem_type, tmp_n_elem, HECMW_INT, j, HECMW_ANY_TAG,
                   VIS_COMM, &stat);
        HECMW_Recv(tmp_elem_node_index, tmp_n_elem + 1, HECMW_INT, j,
                   HECMW_ANY_TAG, VIS_COMM, &stat);

        tmp_elem_node_item =
            (int *)HECMW_calloc(tmp_elem_node_index[tmp_n_elem], sizeof(int));
        if (tmp_elem_node_item == NULL)
          HECMW_vis_memory_exit("tmp recv: elem_node_item");

        HECMW_Recv(tmp_elem_node_item, tmp_elem_node_index[tmp_n_elem],
                   HECMW_INT, j, HECMW_ANY_TAG, VIS_COMM, &stat);

        if (flag_global_ID) {
          tmp_global_elem_ID = (int *)HECMW_calloc(tmp_n_elem, sizeof(int));
          tmp_global_node_ID = (int *)HECMW_calloc(tmp_n_node, sizeof(int));
          if ((tmp_global_elem_ID == NULL) || (tmp_global_node_ID == NULL))
            HECMW_vis_memory_exit("tmp recv: global_elem_ID, global_node_ID");

          HECMW_Recv(tmp_global_elem_ID, tmp_n_elem, HECMW_INT, j,
                     HECMW_ANY_TAG, VIS_COMM, &stat);
          HECMW_Recv(tmp_global_node_ID, tmp_n_node, HECMW_INT, j,
                     HECMW_ANY_TAG, VIS_COMM, &stat);

          tmp_node_ID = NULL;
        } else {
          tmp_global_elem_ID = NULL;
          tmp_global_node_ID = NULL;

          tmp_node_ID = (int *)HECMW_calloc(tmp_n_node * 2, sizeof(int));
          if (tmp_node_ID == NULL) HECMW_vis_memory_exit("tmp recv: node_ID");

          HECMW_Recv(tmp_node_ID, tmp_n_node * 2, HECMW_INT, j, HECMW_ANY_TAG,
                     VIS_COMM, &stat);
        }

        avs_write_elem_conn(outfp, j, tmp_n_elem, tmp_elem_ID, tmp_elem_type,
                            tmp_elem_node_index, tmp_elem_node_item,
                            tmp_global_elem_ID, tmp_global_node_ID, tmp_node_ID,
                            flag_global_ID, flag_skip_ext, flag_oldUCD,
                            nid_offsets, eid_offsets[j]);

        HECMW_free(tmp_elem_ID);
        HECMW_free(tmp_elem_type);
        HECMW_free(tmp_elem_node_index);
        HECMW_free(tmp_elem_node_item);

        if (flag_global_ID) {
          HECMW_free(tmp_global_elem_ID);
          HECMW_free(tmp_global_node_ID);
        } else {
          HECMW_free(tmp_node_ID);
        }
      }
    }
  }

  /* write header for data */
  if (mynode == 0) {
    avs_write_data_header(outfp, tn_component, te_component, flag_oldUCD);
    avs_write_node_data_header(outfp, tn_component, data, flag_Scalar);
  }

  /* write node data */
  if (mynode != 0) {
    HECMW_Send(&mesh->nn_internal, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);
    if (mesh->nn_internal > 0) {
      if (flag_global_ID)
        HECMW_Send(mesh->global_node_ID, mesh->nn_internal, HECMW_INT,
                   MASTER_PE, 0, VIS_COMM);
      HECMW_Send(data->node_val_item, mesh->nn_internal * tn_component,
                 HECMW_DOUBLE, MASTER_PE, 0, VIS_COMM);
    }
  } else {
    int i;

    avs_write_node_data(outfp, mesh->nn_internal, mesh->global_node_ID,
                        tn_component, data->node_val_item, flag_global_ID, 0);

    for (i = 1; i < pesize; i++) {
      int tmp_nn_internal;

      HECMW_Recv(&tmp_nn_internal, 1, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM,
                 &stat);

      if (tmp_nn_internal > 0) {
        int *tmp_global_node_ID = NULL;
        double *tmp_node_val_item;

        if (flag_global_ID) {
          tmp_global_node_ID =
              (int *)HECMW_calloc(tmp_nn_internal, sizeof(int));
          if (tmp_global_node_ID == NULL)
            HECMW_vis_memory_exit("tmp recv: global_node_ID (for data)");

          HECMW_Recv(tmp_global_node_ID, tmp_nn_internal, HECMW_INT, i,
                     HECMW_ANY_TAG, VIS_COMM, &stat);
        }

        tmp_node_val_item = (double *)HECMW_calloc(
            tmp_nn_internal * tn_component, sizeof(double));
        if (tmp_node_val_item == NULL)
          HECMW_vis_memory_exit("tmp recv: node_val_item");

        HECMW_Recv(tmp_node_val_item, tmp_nn_internal * tn_component,
                   HECMW_DOUBLE, i, HECMW_ANY_TAG, VIS_COMM, &stat);

        avs_write_node_data(outfp, tmp_nn_internal, tmp_global_node_ID,
                            tn_component, tmp_node_val_item, flag_global_ID,
                            nid_offsets[i]);

        if (flag_global_ID) HECMW_free(tmp_global_node_ID);
        HECMW_free(tmp_node_val_item);
      }
    }
  }

  /* write elem data header*/
  if (mynode == 0){
    if (te_component > 0){
      avs_write_elem_data_header(outfp, te_component, data, flag_Scalar);
    }
  }

  /* write elem data */
  if (mynode != 0) {
    if (te_component > 0){
      HECMW_Send (&mesh->n_elem, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);
      if (mesh->n_elem > 0){
        HECMW_Send (mesh->elem_ID, mesh->n_elem * 2, HECMW_INT, MASTER_PE, 0, VIS_COMM);
        if (flag_global_ID)
          HECMW_Send (mesh->global_elem_ID, mesh->n_elem, HECMW_INT, MASTER_PE, 0, VIS_COMM);
        HECMW_Send (data->elem_val_item, mesh->n_elem * te_component, HECMW_DOUBLE, MASTER_PE, 0, VIS_COMM);
      }
    }
  } else {
    int i;
    if (te_component > 0){
      avs_write_elem_data(outfp, mynode, mesh->n_elem, mesh->elem_ID, mesh->global_elem_ID, te_component, data->elem_val_item, flag_global_ID, 0);
      for (i = 1; i < pesize; i++){
        int tmp_n_elem;
        HECMW_Recv (&tmp_n_elem, 1, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM, &stat);

        if (tmp_n_elem > 0){
          int *tmp_elem_ID;
          double *tmp_elem_val_item;
          int *tmp_global_elem_ID = NULL;

          tmp_elem_ID = (int *) HECMW_calloc (tmp_n_elem * 2, sizeof (int));
          if (tmp_elem_ID == NULL)
            HECMW_vis_memory_exit ("tmp recv: elem_ID");

          HECMW_Recv (tmp_elem_ID, tmp_n_elem * 2, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM, &stat);

          if (flag_global_ID){
            tmp_global_elem_ID = (int *) HECMW_calloc (tmp_n_elem, sizeof (int));
            if (tmp_global_elem_ID == NULL)
              HECMW_vis_memory_exit ("tmp recv: global_elem_ID (for data)");

            HECMW_Recv (tmp_global_elem_ID, tmp_n_elem, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM, &stat);
          }

          tmp_elem_val_item = (double *) HECMW_calloc (tmp_n_elem * te_component, sizeof (double));
          if (tmp_elem_val_item == NULL)
            HECMW_vis_memory_exit ("tmp recv: elem_val_item");

          HECMW_Recv (tmp_elem_val_item, tmp_n_elem * te_component, HECMW_DOUBLE, i, HECMW_ANY_TAG, VIS_COMM, &stat);

          avs_write_elem_data(outfp, i, tmp_n_elem, tmp_elem_ID, tmp_global_elem_ID, te_component, tmp_elem_val_item, flag_global_ID, eid_offsets[i]);

          HECMW_free (tmp_elem_ID);
          if (flag_global_ID)
            HECMW_free (tmp_global_elem_ID);
          HECMW_free (tmp_elem_val_item);
        }
      }
    }
  }

  if (mynode == 0) {
    HECMW_free(nid_offsets);
    HECMW_free(eid_offsets);
    fclose(outfp);
  }

  return;
}

void HECMW_avs_output(struct hecmwST_local_mesh *mesh,
                      struct hecmwST_result_data *data, char *outfile,
                      HECMW_Comm VIS_COMM) {
  int flag_oldUCD    = 0;
  int flag_global_ID = 1;
  int flag_Scalar    = 0;
  avs_output(mesh, data, outfile, VIS_COMM, flag_oldUCD, flag_global_ID,
             flag_Scalar);
}

void HECMW_reorder_avs_output(struct hecmwST_local_mesh *mesh,
                              struct hecmwST_result_data *data, char *outfile,
                              HECMW_Comm VIS_COMM) {
  int flag_oldUCD    = 1;
  int flag_global_ID = 0;
  int flag_Scalar    = 0;
  avs_output(mesh, data, outfile, VIS_COMM, flag_oldUCD, flag_global_ID,
             flag_Scalar);
}

void HECMW_microavs_output(struct hecmwST_local_mesh *mesh,
                           struct hecmwST_result_data *data, char *outfile,
                           HECMW_Comm VIS_COMM) {
  int flag_oldUCD    = 0;
  int flag_global_ID = 1;
  int flag_Scalar    = 1;
  avs_output(mesh, data, outfile, VIS_COMM, flag_oldUCD, flag_global_ID,
             flag_Scalar);
}

void HECMW_bin_avs_output(struct hecmwST_local_mesh *mesh,
                          struct hecmwST_result_data *data, char *outfile,
                          HECMW_Comm VIS_COMM) {
  int i, j, k, ii, m;
  int mynode, pesize;
  HECMW_Status stat;
  int ielm, nn[20], tmp_int, tn_component, tmp_int2, te_component, tmp_nn[20];
  double *tmp_recv_d, *tmp_send_d;
  int *tmp_recv_i, *tmp_elem_ID, *tmp_elem_type, *tmp_elem_global_ID,
      *tmp_elem_node_index, *tmp_elem_node_item, *tmp_section_ID,
      *tmp_node_global_ID;
  FILE *fp, *fp2;
  int total_n_node, total_n_elem;

  int icell_type[1];
  float xyz[3];
  char keyword[7];
  char title[70];
  float version;
  int stepno;
  float timeval;
  int num_nodeveclen;
  int num_celldata;
  int null_flag;
  float null_value;
  char nodedata_label[16];
  char nodedata_unit[16];
  int ztype;
  int tmp_int_conv, *tmp_global_id, count, elem_type_bin;
  float tmp_d_conv;
  int node_num;

  HECMW_Comm_rank(VIS_COMM, &mynode);
  HECMW_Comm_size(VIS_COMM, &pesize);
  if ((mesh->elem_type[0] == 231) || (mesh->elem_type[0] == 232) ||
      (mesh->elem_type[0] == 731) || (mesh->elem_type[0] == 732)) {
    elem_type_bin = 2;
    node_num      = 3;
  } else if ((mesh->elem_type[0] == 241) || (mesh->elem_type[0] == 242) ||
             (mesh->elem_type[0] == 741) || (mesh->elem_type[0] == 742) ||
             (mesh->elem_type[0] == 743)) {
    elem_type_bin = 3;
    node_num      = 4;
  } else if ((mesh->elem_type[0] == 341) || (mesh->elem_type[0] == 342) ||
             (mesh->elem_type[0] == 3414)) {
    elem_type_bin = 4;
    node_num      = 4;
  } else if ((mesh->elem_type[0] == 351) || (mesh->elem_type[0] == 352)) {
    elem_type_bin = 6;
    node_num      = 6;
  } else if ((mesh->elem_type[0] == 361) || (mesh->elem_type[0] == 362)) {
    elem_type_bin = 7;
    node_num      = 8;
  }
  if (pesize > 1) {
    HECMW_Allreduce(&mesh->nn_internal, &total_n_node, 1, HECMW_INT, HECMW_SUM,
                    VIS_COMM);
    HECMW_Allreduce(&mesh->ne_internal, &total_n_elem, 1, HECMW_INT, HECMW_SUM,
                    VIS_COMM);
  } else {
    total_n_node = mesh->nn_internal;
    total_n_elem = mesh->n_elem;
  }
  tn_component = 0;
  for (i = 0; i < data->nn_component; i++) tn_component += data->nn_dof[i];
  te_component = 0;
  for (i = 0; i < data->ne_component; i++) te_component += data->ne_dof[i];
  if (mynode == 0) {
    fp2 = fopen(outfile, "w");
    fprintf(fp2, "#UCD Binary format\n");
    fprintf(fp2, "#\n");
    fprintf(fp2, "data\n");
    fprintf(fp2, "bin_data.dat\n");
    fclose(fp2);

    fp = fopen("bin_data.dat", "wb");

    if (fp == NULL)
      HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0009: Cannot open output file");
    fprintf(stderr, "Start writing binary output for AVS UCD format\n");
    snprintf(keyword, 7, "AVS UCD");
    version = 1.0;
#ifdef CONVERSE_ORDER
    SWAP_FLOAT(version);
#endif
    strcpy(title, "ucd binary test data");
    stepno = 1;
#ifdef CONVERSE_ORDER
    SWAP_INT(stepno);
#endif
    timeval = 1.0;
#ifdef CONVERSE_ORDER
    SWAP_FLOAT(timeval);
#endif
    tmp_int_conv = total_n_node;
#ifdef CONVERSE_ORDER
    SWAP_INT(tmp_int_conv);
#endif

    fwrite(keyword, 7, 1, fp);
    fwrite(&version, 4, 1, fp);
    fwrite(title, 70, 1, fp);
    fwrite(&stepno, 4, 1, fp);
    fwrite(&timeval, 4, 1, fp);
    fwrite(&tmp_int_conv, 4, 1, fp);
    ztype = 1;
#ifdef CONVERSE_ORDER
    SWAP_INT(ztype);
#endif

    fwrite(&ztype, 4, 1, fp);
  }

  if (mynode != 0) {
    HECMW_Send(&mesh->nn_internal, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);
    if (mesh->nn_internal > 0) {
      HECMW_Send(mesh->global_node_ID, mesh->nn_internal, HECMW_INT, MASTER_PE,
                 0, VIS_COMM);
      HECMW_Send(mesh->node, mesh->nn_internal * 3, HECMW_DOUBLE, MASTER_PE, 0,
                 VIS_COMM);
    }
  }
  if (mynode == 0) {
    for (i = 0; i < mesh->nn_internal; i++) {
      tmp_int_conv = mesh->global_node_ID[i];
#ifdef CONVERSE_ORDER
      SWAP_INT(tmp_int_conv);
#endif
      fwrite(&tmp_int_conv, 4, 1, fp);
      xyz[0] = (float)mesh->node[i * 3];
      xyz[1] = (float)mesh->node[i * 3 + 1];
      xyz[2] = (float)mesh->node[i * 3 + 2];
#ifdef CONVERSE_ORDER
      SWAP_FLOAT(xyz[0]);
      SWAP_FLOAT(xyz[1]);
      SWAP_FLOAT(xyz[2]);
#endif
      fwrite(xyz, 4, 3, fp);
    }
    for (i = 1; i < pesize; i++) {
      HECMW_Recv(&tmp_int, 1, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM, &stat);
      if (tmp_int > 0) {
        tmp_recv_i = (int *)HECMW_calloc(tmp_int, sizeof(int));
        tmp_recv_d = (double *)HECMW_calloc(tmp_int * 3, sizeof(double));
        if ((tmp_recv_i == NULL) || (tmp_recv_d == NULL))
          HECMW_vis_memory_exit("tmp_recv");
        HECMW_Recv(tmp_recv_i, tmp_int, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM,
                   &stat);
        HECMW_Recv(tmp_recv_d, tmp_int * 3, HECMW_DOUBLE, i, HECMW_ANY_TAG,
                   VIS_COMM, &stat);
        for (j = 0; j < tmp_int; j++) {
          tmp_int_conv = tmp_recv_i[j];
#ifdef CONVERSE_ORDER
          SWAP_INT(tmp_int_conv);
#endif
          fwrite(&tmp_int_conv, 4, 1, fp);
          xyz[0] = (float)tmp_recv_d[j * 3];
          xyz[1] = (float)tmp_recv_d[j * 3 + 1];
          xyz[2] = (float)tmp_recv_d[j * 3 + 2];
#ifdef CONVERSE_ORDER
          SWAP_FLOAT(xyz[0]);
          SWAP_FLOAT(xyz[1]);
          SWAP_FLOAT(xyz[2]);
#endif
          fwrite(xyz, 4, 3, fp);
        }
        HECMW_free(tmp_recv_i);
        HECMW_free(tmp_recv_d);
      }
    }
  }
  if (mynode != 0) {
    HECMW_Send(&mesh->n_elem, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);
    if (mesh->n_elem > 0) {
      HECMW_Send(mesh->global_elem_ID, mesh->n_elem, HECMW_INT, MASTER_PE, 0,
                 VIS_COMM);
      HECMW_Send(mesh->elem_ID, mesh->n_elem * 2, HECMW_INT, MASTER_PE, 0,
                 VIS_COMM);
    }
  }
  if (mynode == 0) {
    tmp_global_id = (int *)HECMW_calloc(total_n_elem, sizeof(int));
    if (tmp_global_id == NULL) HECMW_vis_memory_exit("tmp_global_id");
    count = 0;

    for (i = 0; i < mesh->n_elem; i++) {
      if (mesh->elem_ID[i * 2 + 1] == mynode) {
        tmp_global_id[count] = mesh->global_elem_ID[i];
        count++;
      }
    }
    for (j = 1; j < pesize; j++) {
      HECMW_Recv(&tmp_int, 1, HECMW_INT, j, HECMW_ANY_TAG, VIS_COMM, &stat);
      if (tmp_int > 0) {
        tmp_elem_global_ID = (int *)HECMW_calloc(tmp_int, sizeof(int));
        tmp_elem_ID        = (int *)HECMW_calloc(tmp_int * 2, sizeof(int));
        if ((tmp_elem_global_ID == NULL) || (tmp_elem_ID == NULL))
          HECMW_vis_memory_exit("tmp recv");
        HECMW_Recv(tmp_elem_global_ID, tmp_int, HECMW_INT, j, HECMW_ANY_TAG,
                   VIS_COMM, &stat);
        HECMW_Recv(tmp_elem_ID, tmp_int * 2, HECMW_INT, j, HECMW_ANY_TAG,
                   VIS_COMM, &stat);

        for (i = 0; i < tmp_int; i++) {
          if (tmp_elem_ID[i * 2 + 1] == j) {
            tmp_global_id[count] = tmp_elem_global_ID[i];
            count++;
          }
        }
        HECMW_free(tmp_elem_global_ID);
        HECMW_free(tmp_elem_ID);
      }
    }
  }

  if (mynode != 0) {
    HECMW_Send(&mesh->n_elem, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);
    if (mesh->n_elem > 0) {
      HECMW_Send(&mesh->n_node, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);

      HECMW_Send(mesh->global_elem_ID, mesh->n_elem, HECMW_INT, MASTER_PE, 0,
                 VIS_COMM);
      HECMW_Send(mesh->elem_ID, mesh->n_elem * 2, HECMW_INT, MASTER_PE, 0,
                 VIS_COMM);
      HECMW_Send(mesh->elem_node_index, mesh->n_elem + 1, HECMW_INT, MASTER_PE,
                 0, VIS_COMM);
      HECMW_Send(mesh->elem_node_item, mesh->elem_node_index[mesh->n_elem],
                 HECMW_INT, MASTER_PE, 0, VIS_COMM);
      HECMW_Send(mesh->global_node_ID, mesh->n_node, HECMW_INT, MASTER_PE, 0,
                 VIS_COMM);
    }
  }

  if (mynode == 0) {
    tmp_int_conv = total_n_elem;
#ifdef CONVERSE_ORDER
    SWAP_INT(tmp_int_conv);
#endif
    fwrite(&tmp_int_conv, 4, 1, fp);
    for (i = 0; i < total_n_elem; i++) {
      tmp_int_conv = tmp_global_id[i];
#ifdef CONVERSE_ORDER
      SWAP_INT(tmp_int_conv);
#endif
      fwrite(&tmp_int_conv, 4, 1, fp);
    }
    for (i = 0; i < total_n_elem; i++) {
      tmp_int_conv = 1;
#ifdef CONVERSE_ORDER
      SWAP_INT(tmp_int_conv);
#endif
      fwrite(&tmp_int_conv, 4, 1, fp);
    }
    icell_type[0] = (char)elem_type_bin;
    for (i = 0; i < total_n_elem; i++) {
      fwrite(&icell_type[0], 1, 1, fp);
    }

    for (i = 0; i < mesh->n_elem; i++) {
      if (mesh->elem_ID[i * 2 + 1] == mynode) {
        ielm = mesh->global_elem_ID[i];
        for (j = 0; j < node_num; j++)
          tmp_nn[j] =
              mesh->global_node_ID
                  [mesh->elem_node_item[mesh->elem_node_index[i] + j] - 1];
        /*		if(mesh->elem_type[0]==342) {
for(m=0;m<4;m++)
nn[m]=tmp_nn[m];
nn[4]=tmp_nn[6];
nn[5]=tmp_nn[5];
nn[6]=tmp_nn[7];
nn[7]=tmp_nn[4];
nn[8]=tmp_nn[9];
nn[9]=tmp_nn[8];
}
else {
         */
        for (m = 0; m < node_num; m++) nn[m] = tmp_nn[m];
        /*		}
         */

        for (m = 0; m < node_num; m++) {
          tmp_int_conv = nn[m];
#ifdef CONVERSE_ORDER
          SWAP_INT(tmp_int_conv);
#endif
          fwrite(&tmp_int_conv, 4, 1, fp);
        }
      }
    }
    for (j = 1; j < pesize; j++) {
      HECMW_Recv(&tmp_int, 1, HECMW_INT, j, HECMW_ANY_TAG, VIS_COMM, &stat);
      if (tmp_int > 0) {
        HECMW_Recv(&tmp_int2, 1, HECMW_INT, j, HECMW_ANY_TAG, VIS_COMM, &stat);
        tmp_elem_global_ID  = (int *)HECMW_calloc(tmp_int, sizeof(int));
        tmp_elem_ID         = (int *)HECMW_calloc(tmp_int * 2, sizeof(int));
        tmp_elem_type       = (int *)HECMW_calloc(tmp_int, sizeof(int));
        tmp_section_ID      = (int *)HECMW_calloc(tmp_int, sizeof(int));
        tmp_elem_node_index = (int *)HECMW_calloc(tmp_int + 1, sizeof(int));
        if ((tmp_elem_global_ID == NULL) || (tmp_elem_ID == NULL) ||
            (tmp_elem_type == NULL) || (tmp_section_ID == NULL))
          HECMW_vis_memory_exit("tmp recv");
        HECMW_Recv(tmp_elem_global_ID, tmp_int, HECMW_INT, j, HECMW_ANY_TAG,
                   VIS_COMM, &stat);
        HECMW_Recv(tmp_elem_ID, tmp_int * 2, HECMW_INT, j, HECMW_ANY_TAG,
                   VIS_COMM, &stat);
        HECMW_Recv(tmp_elem_node_index, tmp_int + 1, HECMW_INT, j,
                   HECMW_ANY_TAG, VIS_COMM, &stat);
        tmp_elem_node_item =
            (int *)HECMW_calloc(tmp_elem_node_index[tmp_int], sizeof(int));
        tmp_node_global_ID = (int *)HECMW_calloc(tmp_int2, sizeof(int));
        HECMW_Recv(tmp_elem_node_item, tmp_elem_node_index[tmp_int], HECMW_INT,
                   j, HECMW_ANY_TAG, VIS_COMM, &stat);
        HECMW_Recv(tmp_node_global_ID, tmp_int2, HECMW_INT, j, HECMW_ANY_TAG,
                   VIS_COMM, &stat);

        for (i = 0; i < tmp_int; i++) {
          if (tmp_elem_ID[i * 2 + 1] == j) {
            ielm = tmp_elem_global_ID[i];
            for (m      = 0; m < node_num; m++)
              tmp_nn[m] = tmp_node_global_ID
                  [tmp_elem_node_item[tmp_elem_node_index[i] + m] - 1];
            /*		if(mesh->elem_type[0]==342) {
for(m=0;m<4;m++)
nn[m]=tmp_nn[m];
nn[4]=tmp_nn[6];
nn[5]=tmp_nn[5];
nn[6]=tmp_nn[7];
nn[7]=tmp_nn[4];
nn[8]=tmp_nn[9];
nn[9]=tmp_nn[8];
}
else {
             */
            for (m = 0; m < node_num; m++) nn[m] = tmp_nn[m];
            /*		}
             */
            for (m = 0; m < node_num; m++) {
              tmp_int_conv = nn[m];
#ifdef CONVERSE_ORDER
              SWAP_INT(tmp_int_conv);
#endif
              fwrite(&tmp_int_conv, 4, 1, fp);
            }
          }
        }
        HECMW_free(tmp_elem_global_ID);
        HECMW_free(tmp_elem_ID);
        HECMW_free(tmp_elem_type);
        HECMW_free(tmp_section_ID);
        HECMW_free(tmp_elem_node_index);
        HECMW_free(tmp_elem_node_item);
        HECMW_free(tmp_node_global_ID);
      }
    }
    HECMW_free(tmp_global_id);
  }

  if (mynode == 0) {
    tmp_int_conv = tn_component;
#ifdef CONVERSE_ORDER
    SWAP_INT(tmp_int_conv);
#endif
    fwrite(&tmp_int_conv, 4, 1, fp);
    tmp_int_conv = 1;
#ifdef CONVERSE_ORDER
    SWAP_INT(tmp_int_conv);
#endif
    fwrite(&tmp_int_conv, 4, 1, fp);
    for (j = 0; j < data->nn_component; j++) {
      for (ii = 0; ii < data->nn_dof[j]; ii++) {
        sprintf(nodedata_label, "%s%d\n", data->node_label[j], ii + 1);
        sprintf(nodedata_unit, "%s", "m/s");
        num_nodeveclen = 1;
#ifdef CONVERSE_ORDER
        SWAP_INT(num_nodeveclen);
#endif
        null_flag = 0;
#ifdef CONVERSE_ORDER
        SWAP_INT(null_flag);
#endif
        null_value = -999.999;
#ifdef CONVERSE_ORDER
        SWAP_FLOAT(null_value);
#endif
        fwrite(nodedata_label, 16, 1, fp);
        fwrite(nodedata_unit, 16, 1, fp);
        fwrite(&num_nodeveclen, 4, 1, fp);
        fwrite(&null_flag, 4, 1, fp);
        fwrite(&null_value, 4, 1, fp);
      }
    }
  }
  if (mynode != 0) {
    HECMW_Send(&mesh->nn_internal, 1, HECMW_INT, MASTER_PE, 0, VIS_COMM);
    if (mesh->nn_internal > 0) {
      HECMW_Send(mesh->global_node_ID, mesh->nn_internal, HECMW_INT, MASTER_PE,
                 0, VIS_COMM);
      tmp_send_d = (double *)HECMW_calloc(mesh->nn_internal * tn_component,
                                          sizeof(double));
      for (i          = 0; i < mesh->nn_internal * tn_component; i++)
        tmp_send_d[i] = data->node_val_item[i];

      HECMW_Send(tmp_send_d, mesh->nn_internal * tn_component, HECMW_DOUBLE,
                 MASTER_PE, 0, VIS_COMM);

      HECMW_free(tmp_send_d);
    }
  }

  if (mynode == 0) {
    for (i = 0; i < mesh->nn_internal; i++) {
      for (k = 0; k < tn_component; k++) {
        tmp_d_conv = (float)data->node_val_item[i * tn_component + k];
#ifdef CONVERSE_ORDER
        SWAP_FLOAT(tmp_d_conv);
#endif
        fwrite(&tmp_d_conv, 4, 1, fp);
      }
    }
    for (i = 1; i < pesize; i++) {
      HECMW_Recv(&tmp_int, 1, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM, &stat);
      if (tmp_int > 0) {
        tmp_recv_i = (int *)HECMW_calloc(tmp_int, sizeof(int));
        tmp_recv_d =
            (double *)HECMW_calloc(tmp_int * tn_component, sizeof(double));
        if ((tmp_recv_i == NULL) || (tmp_recv_d == NULL))
          HECMW_vis_memory_exit("tmp_recv");
        HECMW_Recv(tmp_recv_i, tmp_int, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM,
                   &stat);
        HECMW_Recv(tmp_recv_d, tmp_int * tn_component, HECMW_DOUBLE, i,
                   HECMW_ANY_TAG, VIS_COMM, &stat);
        for (j = 0; j < tmp_int; j++) {
          for (k = 0; k < tn_component; k++) {
            tmp_d_conv = (float)tmp_recv_d[j * tn_component + k];
#ifdef CONVERSE_ORDER
            SWAP_FLOAT(tmp_d_conv);
#endif
            fwrite(&tmp_d_conv, 4, 1, fp);
          }
        }
        HECMW_free(tmp_recv_i);
        HECMW_free(tmp_recv_d);
      }
    }
    num_celldata = 0;
    fwrite(&num_celldata, 4, 1, fp);
  }
  if (mynode == 0) fclose(fp);
  return;
}

void HECMW_separate_avs_output(struct hecmwST_local_mesh *mesh,
                               struct hecmwST_result_data *data,
                               char *outfile) {
  FILE *outfp;
  int mynode;
  int total_n_node, total_n_elem;
  int tn_component, te_component;

  int flag_oldUCD    = 1;
  int flag_global_ID = 0;
  int flag_Scalar    = 0;
  int flag_skip_ext  = 0;

  mynode = mesh->my_rank;

  /* count total number of node/elem */
  total_n_node = mesh->n_node;
  total_n_elem = mesh->n_elem;

  /* count total number of data components */
  count_data_components(data, &tn_component, &te_component);

  /* open file */
  outfp = fopen(outfile, "w");
  if (!outfp)
    HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0009: Cannot open output file");

  /* write header */
  avs_write_header(outfp, total_n_node, total_n_elem, tn_component, 0,
                   flag_oldUCD);

  /* write node coordinate */
  avs_write_node_coord(outfp, mesh->n_node, mesh->global_node_ID, mesh->node,
                       flag_global_ID, 0);

  /* write element connectivity */
  avs_write_elem_conn(outfp, mynode, mesh->n_elem, mesh->elem_ID,
                      mesh->elem_type, mesh->elem_node_index,
                      mesh->elem_node_item, mesh->global_elem_ID,
                      mesh->global_node_ID, mesh->node_ID, flag_global_ID,
                      flag_skip_ext, flag_oldUCD, NULL, 0);

  /* write total number of data components */
  avs_write_data_header(outfp, tn_component, 0, flag_oldUCD);

  /* write header for node data */
  avs_write_node_data_header(outfp, tn_component, data, flag_Scalar);

  /* write node data */
  avs_write_node_data(outfp, mesh->n_node, mesh->global_node_ID, tn_component,
                      data->node_val_item, flag_global_ID, 0);

  fclose(outfp);
  return;
}
