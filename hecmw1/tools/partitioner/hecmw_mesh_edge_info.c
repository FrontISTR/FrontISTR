/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>

#include "hecmw_util.h"
#include "hecmw_common.h"

#include "hecmw_part_struct.h"
#include "hecmw_mesh_hash_sort.h"
#include "hecmw_mesh_edge_info.h"

static int edge_info_rod1(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[2], node_index;
  long long int edge[1];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 2; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[0] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_rod2(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[3], node_index;
  long long int edge[2];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 3; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[1], node[2]);
    if (edge[1] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_tri1(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[3], node_index;
  long long int edge[3];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 3; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[1], node[2]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[2], node[0]);
    if (edge[2] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_tri2(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[6], node_index;
  long long int edge[6];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 6; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[5]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[5], node[1]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[1], node[3]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[3], node[2]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[2], node[4]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[4], node[0]);
    if (edge[5] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_qua1(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[4], node_index;
  long long int edge[4];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 4; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[1], node[2]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[2], node[3]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[3], node[0]);
    if (edge[3] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_qua2(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[8], node_index;
  long long int edge[8];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 8; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[4]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[4], node[1]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[1], node[5]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[5], node[2]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[2], node[6]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[6], node[3]);
    if (edge[5] < 0) goto error;
    edge[6] = HECMW_mesh_hsort_edge(node[3], node[7]);
    if (edge[6] < 0) goto error;
    edge[7] = HECMW_mesh_hsort_edge(node[7], node[0]);
    if (edge[7] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_tet1(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[4], node_index;
  long long int edge[6];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 4; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[1], node[2]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[2], node[0]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[0], node[3]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[1], node[3]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[2], node[3]);
    if (edge[5] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_tet2(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[10], node_index;
  long long int edge[12];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 10; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[4]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[4], node[1]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[1], node[5]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[5], node[2]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[2], node[6]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[6], node[0]);
    if (edge[5] < 0) goto error;
    edge[6] = HECMW_mesh_hsort_edge(node[0], node[7]);
    if (edge[6] < 0) goto error;
    edge[7] = HECMW_mesh_hsort_edge(node[7], node[3]);
    if (edge[7] < 0) goto error;
    edge[8] = HECMW_mesh_hsort_edge(node[1], node[8]);
    if (edge[8] < 0) goto error;
    edge[9] = HECMW_mesh_hsort_edge(node[8], node[3]);
    if (edge[9] < 0) goto error;
    edge[10] = HECMW_mesh_hsort_edge(node[2], node[9]);
    if (edge[10] < 0) goto error;
    edge[11] = HECMW_mesh_hsort_edge(node[9], node[3]);
    if (edge[11] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_pyr1(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[5], node_index;
  long long int edge[8];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 5; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[1], node[2]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[2], node[3]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[3], node[0]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[0], node[4]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[1], node[4]);
    if (edge[5] < 0) goto error;
    edge[6] = HECMW_mesh_hsort_edge(node[2], node[4]);
    if (edge[6] < 0) goto error;
    edge[7] = HECMW_mesh_hsort_edge(node[3], node[4]);
    if (edge[7] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_pyr2(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[13], node_index;
  long long int edge[16];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 13; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[5]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[5], node[1]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[1], node[6]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[6], node[2]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[2], node[7]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[7], node[3]);
    if (edge[5] < 0) goto error;
    edge[6] = HECMW_mesh_hsort_edge(node[3], node[8]);
    if (edge[6] < 0) goto error;
    edge[7] = HECMW_mesh_hsort_edge(node[8], node[0]);
    if (edge[7] < 0) goto error;
    edge[8] = HECMW_mesh_hsort_edge(node[0], node[9]);
    if (edge[8] < 0) goto error;
    edge[9] = HECMW_mesh_hsort_edge(node[9], node[4]);
    if (edge[9] < 0) goto error;
    edge[10] = HECMW_mesh_hsort_edge(node[1], node[10]);
    if (edge[10] < 0) goto error;
    edge[11] = HECMW_mesh_hsort_edge(node[10], node[4]);
    if (edge[11] < 0) goto error;
    edge[12] = HECMW_mesh_hsort_edge(node[2], node[11]);
    if (edge[12] < 0) goto error;
    edge[13] = HECMW_mesh_hsort_edge(node[11], node[4]);
    if (edge[13] < 0) goto error;
    edge[14] = HECMW_mesh_hsort_edge(node[3], node[12]);
    if (edge[14] < 0) goto error;
    edge[15] = HECMW_mesh_hsort_edge(node[12], node[4]);
    if (edge[15] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_pri1(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[6], node_index;
  long long int edge[9];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 6; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[1], node[2]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[2], node[0]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[3], node[4]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[4], node[5]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[5], node[3]);
    if (edge[5] < 0) goto error;
    edge[6] = HECMW_mesh_hsort_edge(node[0], node[3]);
    if (edge[6] < 0) goto error;
    edge[7] = HECMW_mesh_hsort_edge(node[1], node[4]);
    if (edge[7] < 0) goto error;
    edge[8] = HECMW_mesh_hsort_edge(node[2], node[5]);
    if (edge[8] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_pri2(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[15], node_index;
  long long int edge[18];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 15; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[8]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[8], node[1]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[1], node[6]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[6], node[2]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[2], node[7]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[7], node[0]);
    if (edge[5] < 0) goto error;
    edge[6] = HECMW_mesh_hsort_edge(node[3], node[11]);
    if (edge[6] < 0) goto error;
    edge[7] = HECMW_mesh_hsort_edge(node[11], node[4]);
    if (edge[7] < 0) goto error;
    edge[8] = HECMW_mesh_hsort_edge(node[4], node[9]);
    if (edge[8] < 0) goto error;
    edge[9] = HECMW_mesh_hsort_edge(node[9], node[5]);
    if (edge[9] < 0) goto error;
    edge[10] = HECMW_mesh_hsort_edge(node[5], node[10]);
    if (edge[10] < 0) goto error;
    edge[11] = HECMW_mesh_hsort_edge(node[10], node[3]);
    if (edge[11] < 0) goto error;
    edge[12] = HECMW_mesh_hsort_edge(node[0], node[12]);
    if (edge[12] < 0) goto error;
    edge[13] = HECMW_mesh_hsort_edge(node[12], node[3]);
    if (edge[13] < 0) goto error;
    edge[14] = HECMW_mesh_hsort_edge(node[1], node[13]);
    if (edge[14] < 0) goto error;
    edge[15] = HECMW_mesh_hsort_edge(node[13], node[4]);
    if (edge[15] < 0) goto error;
    edge[16] = HECMW_mesh_hsort_edge(node[2], node[14]);
    if (edge[16] < 0) goto error;
    edge[17] = HECMW_mesh_hsort_edge(node[14], node[5]);
    if (edge[17] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_hex1(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[8], node_index;
  long long int edge[12];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 8; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[1], node[2]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[2], node[3]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[3], node[0]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[4], node[5]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[5], node[6]);
    if (edge[5] < 0) goto error;
    edge[6] = HECMW_mesh_hsort_edge(node[6], node[7]);
    if (edge[6] < 0) goto error;
    edge[7] = HECMW_mesh_hsort_edge(node[7], node[4]);
    if (edge[7] < 0) goto error;
    edge[8] = HECMW_mesh_hsort_edge(node[0], node[4]);
    if (edge[8] < 0) goto error;
    edge[9] = HECMW_mesh_hsort_edge(node[1], node[5]);
    if (edge[9] < 0) goto error;
    edge[10] = HECMW_mesh_hsort_edge(node[2], node[6]);
    if (edge[10] < 0) goto error;
    edge[11] = HECMW_mesh_hsort_edge(node[3], node[7]);
    if (edge[11] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_hex2(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[20], node_index;
  long long int edge[24];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 20; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[8]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[8], node[1]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[1], node[9]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[9], node[2]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[2], node[10]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[10], node[3]);
    if (edge[5] < 0) goto error;
    edge[6] = HECMW_mesh_hsort_edge(node[3], node[11]);
    if (edge[6] < 0) goto error;
    edge[7] = HECMW_mesh_hsort_edge(node[11], node[0]);
    if (edge[7] < 0) goto error;
    edge[8] = HECMW_mesh_hsort_edge(node[4], node[12]);
    if (edge[8] < 0) goto error;
    edge[9] = HECMW_mesh_hsort_edge(node[12], node[5]);
    if (edge[9] < 0) goto error;
    edge[10] = HECMW_mesh_hsort_edge(node[5], node[13]);
    if (edge[10] < 0) goto error;
    edge[11] = HECMW_mesh_hsort_edge(node[13], node[6]);
    if (edge[11] < 0) goto error;
    edge[12] = HECMW_mesh_hsort_edge(node[6], node[14]);
    if (edge[12] < 0) goto error;
    edge[13] = HECMW_mesh_hsort_edge(node[14], node[7]);
    if (edge[13] < 0) goto error;
    edge[14] = HECMW_mesh_hsort_edge(node[7], node[15]);
    if (edge[14] < 0) goto error;
    edge[15] = HECMW_mesh_hsort_edge(node[15], node[4]);
    if (edge[15] < 0) goto error;
    edge[16] = HECMW_mesh_hsort_edge(node[0], node[16]);
    if (edge[16] < 0) goto error;
    edge[17] = HECMW_mesh_hsort_edge(node[16], node[4]);
    if (edge[17] < 0) goto error;
    edge[18] = HECMW_mesh_hsort_edge(node[1], node[17]);
    if (edge[18] < 0) goto error;
    edge[19] = HECMW_mesh_hsort_edge(node[17], node[5]);
    if (edge[19] < 0) goto error;
    edge[20] = HECMW_mesh_hsort_edge(node[2], node[18]);
    if (edge[20] < 0) goto error;
    edge[21] = HECMW_mesh_hsort_edge(node[18], node[6]);
    if (edge[21] < 0) goto error;
    edge[22] = HECMW_mesh_hsort_edge(node[3], node[19]);
    if (edge[22] < 0) goto error;
    edge[23] = HECMW_mesh_hsort_edge(node[19], node[7]);
    if (edge[23] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_mst1(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[4], node_index;
  long long int edge[6];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 4; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[1], node[2]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[2], node[3]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[3], node[1]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[0], node[2]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[0], node[3]);
    if (edge[5] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_mst2(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[7], node_index;
  long long int edge[9];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 7; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[1], node[6]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[6], node[2]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[2], node[4]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[4], node[3]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[3], node[5]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[5], node[1]);
    if (edge[5] < 0) goto error;
    edge[6] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[6] < 0) goto error;
    edge[7] = HECMW_mesh_hsort_edge(node[0], node[2]);
    if (edge[7] < 0) goto error;
    edge[8] = HECMW_mesh_hsort_edge(node[0], node[3]);
    if (edge[8] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_msq1(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[5], node_index;
  long long int edge[8];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 5; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[1], node[2]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[2], node[3]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[3], node[4]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[4], node[1]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[0], node[2]);
    if (edge[5] < 0) goto error;
    edge[6] = HECMW_mesh_hsort_edge(node[0], node[3]);
    if (edge[6] < 0) goto error;
    edge[7] = HECMW_mesh_hsort_edge(node[0], node[4]);
    if (edge[7] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_msq2(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[9], node_index;
  long long int edge[12];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 9; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[1], node[5]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[5], node[2]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[2], node[6]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[6], node[3]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[3], node[7]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[7], node[4]);
    if (edge[5] < 0) goto error;
    edge[6] = HECMW_mesh_hsort_edge(node[4], node[8]);
    if (edge[6] < 0) goto error;
    edge[7] = HECMW_mesh_hsort_edge(node[8], node[1]);
    if (edge[7] < 0) goto error;
    edge[8] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[8] < 0) goto error;
    edge[9] = HECMW_mesh_hsort_edge(node[0], node[2]);
    if (edge[9] < 0) goto error;
    edge[10] = HECMW_mesh_hsort_edge(node[0], node[3]);
    if (edge[10] < 0) goto error;
    edge[11] = HECMW_mesh_hsort_edge(node[0], node[4]);
    if (edge[11] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_jtb1(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[2], node_index;
  long long int edge[1];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 2; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[0] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_jtt1(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[6], node_index;
  long long int edge[9];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 6; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[1], node[2]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[2], node[0]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[3], node[4]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[4], node[5]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[5], node[3]);
    if (edge[5] < 0) goto error;
    edge[6] = HECMW_mesh_hsort_edge(node[0], node[3]);
    if (edge[6] < 0) goto error;
    edge[7] = HECMW_mesh_hsort_edge(node[1], node[4]);
    if (edge[7] < 0) goto error;
    edge[8] = HECMW_mesh_hsort_edge(node[2], node[5]);
    if (edge[8] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_jtt2(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[12], node_index;
  long long int edge[15];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 12; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[8]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[8], node[1]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[1], node[6]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[6], node[2]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[2], node[7]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[7], node[0]);
    if (edge[5] < 0) goto error;
    edge[6] = HECMW_mesh_hsort_edge(node[3], node[11]);
    if (edge[6] < 0) goto error;
    edge[7] = HECMW_mesh_hsort_edge(node[11], node[4]);
    if (edge[7] < 0) goto error;
    edge[8] = HECMW_mesh_hsort_edge(node[4], node[9]);
    if (edge[8] < 0) goto error;
    edge[9] = HECMW_mesh_hsort_edge(node[9], node[5]);
    if (edge[9] < 0) goto error;
    edge[10] = HECMW_mesh_hsort_edge(node[5], node[10]);
    if (edge[10] < 0) goto error;
    edge[11] = HECMW_mesh_hsort_edge(node[10], node[3]);
    if (edge[11] < 0) goto error;
    edge[12] = HECMW_mesh_hsort_edge(node[0], node[3]);
    if (edge[12] < 0) goto error;
    edge[13] = HECMW_mesh_hsort_edge(node[1], node[4]);
    if (edge[13] < 0) goto error;
    edge[14] = HECMW_mesh_hsort_edge(node[2], node[5]);
    if (edge[14] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_jtq1(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[8], node_index;
  long long int edge[12];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 8; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[1], node[2]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[2], node[3]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[3], node[0]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[4], node[5]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[5], node[6]);
    if (edge[5] < 0) goto error;
    edge[6] = HECMW_mesh_hsort_edge(node[6], node[7]);
    if (edge[6] < 0) goto error;
    edge[7] = HECMW_mesh_hsort_edge(node[7], node[4]);
    if (edge[7] < 0) goto error;
    edge[8] = HECMW_mesh_hsort_edge(node[0], node[4]);
    if (edge[8] < 0) goto error;
    edge[9] = HECMW_mesh_hsort_edge(node[1], node[5]);
    if (edge[9] < 0) goto error;
    edge[10] = HECMW_mesh_hsort_edge(node[2], node[6]);
    if (edge[10] < 0) goto error;
    edge[11] = HECMW_mesh_hsort_edge(node[3], node[7]);
    if (edge[11] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_jtq2(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[16], node_index;
  long long int edge[20];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 16; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[8]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[8], node[1]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[1], node[9]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[9], node[2]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[2], node[10]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[10], node[3]);
    if (edge[5] < 0) goto error;
    edge[6] = HECMW_mesh_hsort_edge(node[3], node[11]);
    if (edge[6] < 0) goto error;
    edge[7] = HECMW_mesh_hsort_edge(node[11], node[0]);
    if (edge[7] < 0) goto error;
    edge[8] = HECMW_mesh_hsort_edge(node[4], node[12]);
    if (edge[8] < 0) goto error;
    edge[9] = HECMW_mesh_hsort_edge(node[12], node[5]);
    if (edge[9] < 0) goto error;
    edge[10] = HECMW_mesh_hsort_edge(node[5], node[13]);
    if (edge[10] < 0) goto error;
    edge[11] = HECMW_mesh_hsort_edge(node[13], node[6]);
    if (edge[11] < 0) goto error;
    edge[12] = HECMW_mesh_hsort_edge(node[6], node[14]);
    if (edge[12] < 0) goto error;
    edge[13] = HECMW_mesh_hsort_edge(node[14], node[7]);
    if (edge[13] < 0) goto error;
    edge[14] = HECMW_mesh_hsort_edge(node[7], node[15]);
    if (edge[14] < 0) goto error;
    edge[15] = HECMW_mesh_hsort_edge(node[15], node[4]);
    if (edge[15] < 0) goto error;
    edge[16] = HECMW_mesh_hsort_edge(node[0], node[4]);
    if (edge[16] < 0) goto error;
    edge[17] = HECMW_mesh_hsort_edge(node[1], node[5]);
    if (edge[17] < 0) goto error;
    edge[18] = HECMW_mesh_hsort_edge(node[2], node[6]);
    if (edge[18] < 0) goto error;
    edge[19] = HECMW_mesh_hsort_edge(node[3], node[7]);
    if (edge[19] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_bem1(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[2], node_index;
  long long int edge[1];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 2; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[0] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_bem2(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[3], node_index;
  long long int edge[2];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 3; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[1], node[2]);
    if (edge[1] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_bem3(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[4], node_index;
  long long int edge[6];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 4; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[1], node[2]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[2], node[0]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[0], node[3]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[1], node[3]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[2], node[3]);
    if (edge[5] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_sht1(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[3], node_index;
  long long int edge[3];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 3; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[1], node[2]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[2], node[0]);
    if (edge[2] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_sht2(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[6], node_index;
  long long int edge[6];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 6; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[5]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[5], node[1]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[1], node[3]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[3], node[2]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[2], node[4]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[4], node[0]);
    if (edge[5] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_shq1(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[4], node_index;
  long long int edge[4];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 4; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[1], node[2]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[2], node[3]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[3], node[0]);
    if (edge[3] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_shq2(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[8], node_index;
  long long int edge[8];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 8; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[4]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[4], node[1]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[1], node[5]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[5], node[2]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[2], node[6]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[6], node[3]);
    if (edge[5] < 0) goto error;
    edge[6] = HECMW_mesh_hsort_edge(node[3], node[7]);
    if (edge[6] < 0) goto error;
    edge[7] = HECMW_mesh_hsort_edge(node[7], node[0]);
    if (edge[7] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_sht6(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[6], node_index;
  long long int edge[9];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 6; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[1], node[2]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[2], node[0]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[3], node[4]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[4], node[5]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[5], node[3]);
    if (edge[5] < 0) goto error;
    edge[6] = HECMW_mesh_hsort_edge(node[0], node[3]);
    if (edge[6] < 0) goto error;
    edge[7] = HECMW_mesh_hsort_edge(node[1], node[4]);
    if (edge[7] < 0) goto error;
    edge[8] = HECMW_mesh_hsort_edge(node[2], node[5]);
    if (edge[8] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

static int edge_info_shq8(struct hecmwST_local_mesh *local_mesh, const int is,
                          const int ie) {
  int node[8], node_index;
  long long int edge[12];
  int i, j;

  for (i = is; i < ie; i++) {
    node_index = local_mesh->elem_node_index[i];
    for (j = 0; j < 8; j++) {
      node[j] = local_mesh->elem_node_item[node_index + j];
    }

    edge[0] = HECMW_mesh_hsort_edge(node[0], node[1]);
    if (edge[0] < 0) goto error;
    edge[1] = HECMW_mesh_hsort_edge(node[1], node[2]);
    if (edge[1] < 0) goto error;
    edge[2] = HECMW_mesh_hsort_edge(node[2], node[3]);
    if (edge[2] < 0) goto error;
    edge[3] = HECMW_mesh_hsort_edge(node[3], node[0]);
    if (edge[3] < 0) goto error;
    edge[4] = HECMW_mesh_hsort_edge(node[4], node[5]);
    if (edge[4] < 0) goto error;
    edge[5] = HECMW_mesh_hsort_edge(node[5], node[6]);
    if (edge[5] < 0) goto error;
    edge[6] = HECMW_mesh_hsort_edge(node[6], node[7]);
    if (edge[6] < 0) goto error;
    edge[7] = HECMW_mesh_hsort_edge(node[7], node[4]);
    if (edge[7] < 0) goto error;
    edge[8] = HECMW_mesh_hsort_edge(node[0], node[4]);
    if (edge[8] < 0) goto error;
    edge[9] = HECMW_mesh_hsort_edge(node[1], node[5]);
    if (edge[9] < 0) goto error;
    edge[10] = HECMW_mesh_hsort_edge(node[2], node[6]);
    if (edge[10] < 0) goto error;
    edge[11] = HECMW_mesh_hsort_edge(node[3], node[7]);
    if (edge[11] < 0) goto error;
  }

  return 0;

error:
  return -1;
}

extern int HECMW_mesh_edge_info(struct hecmwST_local_mesh *local_mesh,
                                struct hecmw_part_edge_data *edge_data) {
  int rtc;
  int i, is, ie;

  if (local_mesh == NULL) {
    HECMW_set_error(HECMW_PART_E_NULL_POINTER, "\'local_mesh\' is NULL");
    goto error;
  }
  if (edge_data == NULL) {
    HECMW_set_error(HECMW_PART_E_NULL_POINTER, "\'edge_data\' is NULL");
    goto error;
  }

  rtc = HECMW_mesh_hsort_edge_init(local_mesh->n_node, local_mesh->n_elem);
  if (rtc != 0) goto error;

  for (i = 0; i < local_mesh->n_elem_type; i++) {
    is = local_mesh->elem_type_index[i];
    ie = local_mesh->elem_type_index[i + 1];

    switch (local_mesh->elem_type_item[i]) {
      case HECMW_ETYPE_ROD1: /* line ( 1st order ) */
      case HECMW_ETYPE_ROD31:
      case HECMW_ETYPE_SPGDPT1:
        if (edge_info_rod1(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_ROD2: /* line ( 2nd order ) */
        if (edge_info_rod2(local_mesh, is, ie)) goto error;
        break;

      case HECMW_ETYPE_TRI1: /* triangle ( 1st order ) */
        if (edge_info_tri1(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_TRI2: /* triangle ( 2nd order ) */
        if (edge_info_tri2(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_QUA1: /* quadrangle ( 1st order ) */
        if (edge_info_qua1(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_QUA2: /* quadrangle ( 2nd order ) */
        if (edge_info_qua2(local_mesh, is, ie)) goto error;
        break;

      case HECMW_ETYPE_TET1: /* tetrahedron ( 1st order ) */
        if (edge_info_tet1(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_TET1_4: /* tetrahedron ( 1st order ) */
        if (edge_info_tet1(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_TET2: /* tetrahedron ( 2nd order ) */
        if (edge_info_tet2(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_PYR1: /* pyramid ( 1st order ) */
        if (edge_info_pyr1(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_PYR2: /* pyramid ( 2nd order ) */
        if (edge_info_pyr2(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_PRI1: /* prism ( 1st order ) */
        if (edge_info_pri1(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_PRI2: /* prism ( 2nd order ) */
        if (edge_info_pri2(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_HEX1: /* hexahedron ( 1st order ) */
        if (edge_info_hex1(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_HEX1_4: /* hexahedron ( 1st order ) */
        if (edge_info_hex1(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_HEX2: /* hexahedron ( 2nd order ) */
        if (edge_info_hex2(local_mesh, is, ie)) goto error;
        break;

      case HECMW_ETYPE_MST1: /* triangluar master-slave type ( 1st order ) */
        if (edge_info_mst1(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_MST2: /* triangluar master-slave type ( 2nd order ) */
        if (edge_info_mst2(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_MSQ1: /* quadrilateral master-slave type ( 1st order ) */
        if (edge_info_msq1(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_MSQ2: /* quadrilateral master-slave type ( 2nd order ) */
        if (edge_info_msq2(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_JTT1: /* triangular interface ( 1st order ) */
        if (edge_info_jtt1(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_JTT2: /* triangluar interface ( 2nd order ) */
        if (edge_info_jtt2(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_JTQ1: /* quadrilateral interface ( 1st order ) */
        if (edge_info_jtq1(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_JTQ2: /* quadrilateral interface ( 2nd order ) */
        if (edge_info_jtq2(local_mesh, is, ie)) goto error;
        break;

      case HECMW_ETYPE_BEM1: /* beam ( 1st order ) */
        if (edge_info_bem1(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_BEM2: /* beam ( 2nd order ) */
        if (edge_info_bem2(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_BEM3: /* beam ( Mixed beam 341) */
        if (edge_info_bem3(local_mesh, is, ie)) goto error;
        break;

      case HECMW_ETYPE_SHT1: /* triangluar shell ( 1st order ) */
        if (edge_info_sht1(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_SHT2: /* triangular shell ( 2nd order ) */
        if (edge_info_sht2(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_SHQ1: /* quadrilateral shell ( 1st order ) */
        if (edge_info_shq1(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_SHQ2: /* quadrilateral shell ( 2nd order ) */
        if (edge_info_shq2(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_SHT6: /* triangluar shell ( Mixed solid ) */
        if (edge_info_sht6(local_mesh, is, ie)) goto error;
        break;
      case HECMW_ETYPE_SHQ8: /* quadrilateral shell ( Mixed solid ) */
        if (edge_info_shq8(local_mesh, is, ie)) goto error;
        break;

      case HECMW_ETYPE_PTT1: /* triangular patch ( 1st order ) */
      case HECMW_ETYPE_PTT2: /* triangular patch ( 2nd order ) */
      case HECMW_ETYPE_PTQ1: /* quadrilateral patch ( 1st order ) */
      case HECMW_ETYPE_PTQ2: /* quadrilateral patch ( 2nd order ) */
        break; /* no need to add any edge */

      case HECMW_ETYPE_LN11:
      case HECMW_ETYPE_LN12:
      case HECMW_ETYPE_LN13:
      case HECMW_ETYPE_LN14:
      case HECMW_ETYPE_LN15:
      case HECMW_ETYPE_LN16:
      case HECMW_ETYPE_LN21:
      case HECMW_ETYPE_LN22:
      case HECMW_ETYPE_LN23:
      case HECMW_ETYPE_LN24:
      case HECMW_ETYPE_LN25:
      case HECMW_ETYPE_LN26:
      case HECMW_ETYPE_LN31:
      case HECMW_ETYPE_LN32:
      case HECMW_ETYPE_LN33:
      case HECMW_ETYPE_LN34:
      case HECMW_ETYPE_LN35:
      case HECMW_ETYPE_LN36:
      case HECMW_ETYPE_LN41:
      case HECMW_ETYPE_LN42:
      case HECMW_ETYPE_LN43:
      case HECMW_ETYPE_LN44:
      case HECMW_ETYPE_LN45:
      case HECMW_ETYPE_LN46:
      case HECMW_ETYPE_LN51:
      case HECMW_ETYPE_LN52:
      case HECMW_ETYPE_LN53:
      case HECMW_ETYPE_LN54:
      case HECMW_ETYPE_LN55:
      case HECMW_ETYPE_LN56:
      case HECMW_ETYPE_LN61:
      case HECMW_ETYPE_LN62:
      case HECMW_ETYPE_LN63:
      case HECMW_ETYPE_LN64:
      case HECMW_ETYPE_LN65:
      case HECMW_ETYPE_LN66:
        if (edge_info_rod1(local_mesh, is, ie)) goto error;
        break;

      default:
        HECMW_set_error(HECMW_PART_E_INVALID_ETYPE, "%d",
                        local_mesh->elem_type_item[i]);
        goto error;
    }
  }

  edge_data->n_edge = HECMW_mesh_hsort_edge_get_n();
  if (edge_data->n_edge < 0) goto error;

  edge_data->edge_node_item = HECMW_mesh_hsort_edge_get_v();
  if (edge_data->edge_node_item == NULL) goto error;

  HECMW_mesh_hsort_edge_final();

  return 0;

error:
  HECMW_mesh_hsort_edge_final();

  return -1;
}
