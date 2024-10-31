/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#define INAGAKI_PARTITIONER

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <math.h>

#include "hecmw_util.h"
#include "hecmw_common.h"
#include "hecmw_io.h"

#include "hecmw_part_define.h"
#include "hecmw_part_struct.h"
#include "hecmw_part_log.h"
#include "hecmw_mesh_hash_sort.h"
#include "hecmw_mesh_edge_info.h"
#include "hecmw_part_get_control.h"
#include "hecmw_partition.h"
#include "hecmw_ucd_print.h"
#include "hecmw_graph.h"
#include "hecmw_common_define.h"

#ifdef HECMW_PART_WITH_METIS
#include "metis.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#define INTERNAL 1

#define EXTERNAL 2

#define BOUNDARY 4

#define OVERLAP 8

#define MASK 16

#define MARK 32

#define MY_DOMAIN 1

#define NEIGHBOR_DOMAIN 2

#define MPC_BLOCK 4

#define CANDIDATE 8

#define EPS (1.0E-12)

#define F_1_2 (0.5)

#define F_6_10 (0.6)

#define QSORT_LOWER 50

#define MASK_BIT(map, bit) ((map) |= (bit))

#define EVAL_BIT(map, bit) ((map) & (bit))

#define INV_BIT(map, bit) ((map) ^= (bit))

#define CLEAR_BIT(map, bit) \
  ((map) |= (bit));         \
  ((map) ^= (bit))

#define CLEAR_IEB(map) \
  ((map) |= (7));      \
  ((map) ^= (7))

#define CLEAR_MM(map) \
  ((map) |= (48));    \
  ((map) ^= (48))

#define DSWAP(a, aa) \
  atemp = (a);       \
  (a)   = (aa);      \
  (aa)  = atemp;

#define ISWAP(b, bb) \
  btemp = (b);       \
  (b)   = (bb);      \
  (bb)  = btemp;

#define RTC_NORMAL 0

#define RTC_ERROR (-1)

#define RTC_WARN 1

#define MAX_NODE_SIZE 20

struct link_unit {
  int id;

  struct link_unit *next;
};

struct link_list {
  int n;

  struct link_unit *list;

  struct link_unit *last;
};

/*===== internal/boundary node/element list of each domain =======*/
static int *n_int_nlist = NULL;
static int *n_bnd_nlist = NULL;
static int *n_int_elist = NULL;
static int *n_bnd_elist = NULL;
static int **int_nlist  = NULL;
static int **bnd_nlist  = NULL;
static int **int_elist  = NULL;
static int **bnd_elist  = NULL;
static int **ngrp_idx   = NULL;
static int **ngrp_item  = NULL;
static int **egrp_idx   = NULL;
static int **egrp_item  = NULL;

/*===== speed up (K. Inagaki )=======*/
static int spdup_clear_MMbnd(char *node_flag, char *elem_flag,
                             int current_domain) {
  int i, node, elem;

  for (i = 0; i < n_bnd_nlist[2 * current_domain + 1]; i++) {
    node = bnd_nlist[current_domain][i];
    CLEAR_MM(node_flag[node - 1]);
  }
  for (i = 0; i < n_bnd_elist[2 * current_domain + 1]; i++) {
    elem = bnd_elist[current_domain][i];
    CLEAR_MM(elem_flag[elem - 1]);
  }
  return RTC_NORMAL;
}

static int spdup_clear_IEB(char *node_flag, char *elem_flag,
                           int current_domain) {
  int i, node, elem;

  for (i = 0; i < n_int_nlist[current_domain]; i++) {
    node = int_nlist[current_domain][i];
    CLEAR_IEB(node_flag[node - 1]);
  }
  for (i = 0; i < n_bnd_nlist[2 * current_domain + 1]; i++) {
    node = bnd_nlist[current_domain][i];
    CLEAR_IEB(node_flag[node - 1]);
  }
  for (i = 0; i < n_int_elist[current_domain]; i++) {
    elem = int_elist[current_domain][i];
    CLEAR_IEB(elem_flag[elem - 1]);
  }
  for (i = 0; i < n_bnd_elist[2 * current_domain + 1]; i++) {
    elem = bnd_elist[current_domain][i];
    CLEAR_IEB(elem_flag[elem - 1]);
  }

  return RTC_NORMAL;
}

static int spdup_init_list(const struct hecmwST_local_mesh *global_mesh) {
  int i, j, k;
  int js, je;
  int node, n_domain, domain[20], flag;

  /*init lists for count (calloc) */
  n_int_nlist = (int *)HECMW_calloc(global_mesh->n_subdomain, sizeof(int));
  if (n_int_nlist == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  n_bnd_nlist = (int *)HECMW_calloc(2 * global_mesh->n_subdomain, sizeof(int));
  if (n_bnd_nlist == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  n_int_elist = (int *)HECMW_calloc(global_mesh->n_subdomain, sizeof(int));
  if (n_int_elist == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  n_bnd_elist = (int *)HECMW_calloc(2 * global_mesh->n_subdomain, sizeof(int));
  if (n_bnd_elist == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  int_nlist = (int **)HECMW_malloc(global_mesh->n_subdomain * sizeof(int *));
  if (int_nlist == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  bnd_nlist = (int **)HECMW_malloc(global_mesh->n_subdomain * sizeof(int *));
  if (bnd_nlist == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  int_elist = (int **)HECMW_malloc(global_mesh->n_subdomain * sizeof(int *));
  if (int_elist == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  bnd_elist = (int **)HECMW_malloc(global_mesh->n_subdomain * sizeof(int *));
  if (bnd_elist == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  /* count internal node */
  for (i = 0; i < global_mesh->n_node; i++) {
    n_int_nlist[global_mesh->node_ID[2 * i + 1]]++;
  }

  /*count internal elem */
  for (i = 0; i < global_mesh->n_elem; i++) {
    n_int_elist[global_mesh->elem_ID[2 * i + 1]]++;
  }

  /*count boundary node and elem */
  for (i = 0; i < global_mesh->n_elem; i++) {
    js        = global_mesh->elem_node_index[i];
    je        = global_mesh->elem_node_index[i + 1];
    node      = global_mesh->elem_node_item[js];
    n_domain  = 1;
    domain[0] = global_mesh->node_ID[2 * node - 1];
    for (j = js + 1; j < je; j++) {
      node = global_mesh->elem_node_item[j];
      for (flag = 0, k = 0; k < n_domain; k++) {
        if (global_mesh->node_ID[2 * node - 1] == domain[k]) {
          flag++;
          break;
        }
      }
      if (flag == 0) {
        domain[n_domain] = global_mesh->node_ID[2 * node - 1];
        n_domain++;
      }
    }

    if (n_domain > 1) {
      for (j = 0; j < n_domain; j++) {
        n_bnd_elist[domain[j]]++;
        n_bnd_nlist[domain[j]] += je - js;
      }
    }
  }

  /*allocate node/element list of each domain */
  for (i = 0; i < global_mesh->n_subdomain; i++) {
    int_nlist[i] = (int *)HECMW_calloc(n_int_nlist[i], sizeof(int));
    if (int_nlist[i] == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    }
    bnd_nlist[i] = (int *)HECMW_calloc(n_bnd_nlist[i], sizeof(int));
    if (bnd_nlist[i] == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    }
    int_elist[i] = (int *)HECMW_calloc(n_int_elist[i], sizeof(int));
    if (int_elist[i] == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    }
    bnd_elist[i] = (int *)HECMW_calloc(n_bnd_elist[i], sizeof(int));
    if (bnd_elist[i] == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    }
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int int_cmp(const void *v1, const void *v2) {
  const int *i1, *i2;

  i1 = (const int *)v1;
  i2 = (const int *)v2;

  if (*i1 < *i2) return -1;
  if (*i1 > *i2) return 1;
  return 0;
}

static int get_boundary_nodelist(const struct hecmwST_local_mesh *global_mesh,
                                 int domain) {
  int i, j, k;
  int ks, ke, node, elem, counter;

  for (counter = 0, j = 0; j < n_bnd_elist[2 * domain + 1]; j++) {
    elem = bnd_elist[domain][j];
    ks   = global_mesh->elem_node_index[elem - 1];
    ke   = global_mesh->elem_node_index[elem];
    for (k = ks; k < ke; k++) {
      node                       = global_mesh->elem_node_item[k];
      bnd_nlist[domain][counter] = node;
      counter++;
    }
  }

  qsort(bnd_nlist[domain], counter, sizeof(int), int_cmp);

  if (counter > 1) {
    i = 1;
    for (j = 1; j < counter; j++) {
      if (bnd_nlist[domain][j - 1] != bnd_nlist[domain][j]) {
        bnd_nlist[domain][i] = bnd_nlist[domain][j];
        i++;
      }
    }
  } else {
    i = counter;
  }

  n_bnd_nlist[2 * domain + 1] = i;

  return RTC_NORMAL;
}

static int sort_and_resize_bndlist(const struct hecmwST_local_mesh *global_mesh,
                                   int domain) {
  int i, node, elem;
  int *work = NULL;
  int bnd_and_int, bnd_not_int;
  int n_nlist, n_elist;

  /*boundary node list */
  n_nlist = n_bnd_nlist[2 * domain + 1];
  work    = (int *)HECMW_malloc(n_nlist * sizeof(int));
  if (work == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  /*sort */
  bnd_and_int = 0;
  bnd_not_int = 0;
  for (i = 0; i < n_nlist; i++) {
    node = bnd_nlist[domain][i];
    if (global_mesh->node_ID[2 * node - 1] == domain) {
      work[bnd_and_int] = node;
      bnd_and_int++;
    }
  }
  for (i = 0; i < n_nlist; i++) {
    node = bnd_nlist[domain][i];
    if (global_mesh->node_ID[2 * node - 1] != domain) {
      work[bnd_and_int + bnd_not_int] = node;
      bnd_not_int++;
    }
  }
  n_bnd_nlist[2 * domain]     = bnd_and_int;
  n_bnd_nlist[2 * domain + 1] = bnd_and_int + bnd_not_int;
  HECMW_assert(n_nlist == n_bnd_nlist[2 * domain + 1]);

  /*resize */
  HECMW_free(bnd_nlist[domain]);
  bnd_nlist[domain] = (int *)HECMW_calloc(n_nlist, sizeof(int));
  if (bnd_nlist[domain] == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  for (i = 0; i < n_nlist; i++) {
    bnd_nlist[domain][i] = work[i];
  }
  HECMW_free(work);

  /*boundary element list */
  n_elist = n_bnd_elist[2 * domain + 1];
  work    = (int *)HECMW_malloc(n_elist * sizeof(int));
  if (work == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  /*sort */
  bnd_and_int = 0;
  bnd_not_int = 0;
  for (i = 0; i < n_elist; i++) {
    elem = bnd_elist[domain][i];
    if (global_mesh->elem_ID[2 * elem - 1] == domain) {
      work[bnd_and_int] = elem;
      bnd_and_int++;
    }
  }
  for (i = 0; i < n_elist; i++) {
    elem = bnd_elist[domain][i];
    if (global_mesh->elem_ID[2 * elem - 1] != domain) {
      work[bnd_and_int + bnd_not_int] = elem;
      bnd_not_int++;
    }
  }
  n_bnd_elist[2 * domain]     = bnd_and_int;
  n_bnd_elist[2 * domain + 1] = bnd_and_int + bnd_not_int;
  for (i = 0; i < n_elist; i++) {
    bnd_elist[domain][i] = work[i];
  }
  HECMW_free(work);
  HECMW_assert(n_elist == n_bnd_elist[2 * domain + 1]);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int spdup_make_list(const struct hecmwST_local_mesh *global_mesh) {
  int i, j, k;
  int js, je, ks, ke;
  int node, elem, n_domain, domain[20], flag;
  int current_domain;
  int rtc;

  /*clear counters */
  for (i = 0; i < global_mesh->n_subdomain; i++) {
    n_int_nlist[i]         = 0;
    n_bnd_nlist[2 * i]     = 0;
    n_bnd_nlist[2 * i + 1] = 0;
    n_int_elist[i]         = 0;
    n_bnd_elist[2 * i]     = 0;
    n_bnd_elist[2 * i + 1] = 0;
  }

  /* internal nodelist for each domain */
  for (i = 0; i < global_mesh->n_node; i++) {
    current_domain = global_mesh->node_ID[2 * i + 1];
    int_nlist[current_domain][n_int_nlist[current_domain]] = i + 1;
    n_int_nlist[current_domain]++;
  }

  /* internal elemlist for each domain */
  for (i = 0; i < global_mesh->n_elem; i++) {
    current_domain = global_mesh->elem_ID[2 * i + 1];
    int_elist[current_domain][n_int_elist[current_domain]] = i + 1;
    n_int_elist[current_domain]++;
  }

  /* boundary elemlist for each domain */
  for (i = 0; i < global_mesh->n_elem; i++) {
    js        = global_mesh->elem_node_index[i];
    je        = global_mesh->elem_node_index[i + 1];
    node      = global_mesh->elem_node_item[js];
    n_domain  = 1;
    domain[0] = global_mesh->node_ID[2 * node - 1];
    for (j = js + 1; j < je; j++) {
      node = global_mesh->elem_node_item[j];
      for (flag = 0, k = 0; k < n_domain; k++) {
        if (global_mesh->node_ID[2 * node - 1] == domain[k]) {
          flag++;
          break;
        }
      }
      if (flag == 0) {
        domain[n_domain] = global_mesh->node_ID[2 * node - 1];
        n_domain++;
      }
    }

    if (n_domain > 1) {
      for (j = 0; j < n_domain; j++) {
        bnd_elist[domain[j]][n_bnd_elist[2 * domain[j] + 1]] = i + 1;
        n_bnd_elist[2 * domain[j] + 1]++;
      }
    }
  }

  /* boundary nodelist for each domain */
  for (i = 0; i < global_mesh->n_subdomain; i++) {
    rtc = get_boundary_nodelist(global_mesh, i);
    if (rtc != RTC_NORMAL) goto error;
  }

  for (i = 0; i < global_mesh->n_subdomain; i++) {
    rtc = sort_and_resize_bndlist(global_mesh, i);
    if (rtc != RTC_NORMAL) goto error;
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int spdup_make_node_grouplist(
    const struct hecmwST_local_mesh *global_mesh) {
  struct hecmwST_node_grp *node_group_global = global_mesh->node_group;
  int i, j, k, node, n_bnd, n_out;
  int *n_domain = NULL;
  int **domain  = NULL;
  int current_domain;
  int counter[global_mesh->n_subdomain];

  /*make list of node to domain(both internal and boundary) */
  n_domain = (int *)HECMW_calloc(global_mesh->n_node, sizeof(int));
  if (n_domain == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  /*count outer node(boundary and not internal) */
  for (i = 0; i < global_mesh->n_subdomain; i++) {
    n_bnd = n_bnd_nlist[2 * i];
    n_out = n_bnd_nlist[2 * i + 1] - n_bnd_nlist[2 * i];
    if (n_out == 0) continue;
    for (j = 0; j < n_out; j++) {
      node = bnd_nlist[i][n_bnd + j];
      n_domain[node - 1]++;
    }
  }
  /*make list */
  domain = (int **)HECMW_malloc(global_mesh->n_node * sizeof(int *));
  if (domain == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  for (i = 0; i < global_mesh->n_node; i++) {
    domain[i] = (int *)HECMW_malloc((n_domain[i] + 1) *
                                    sizeof(int)); /*+1 means internal node */
    if (domain[i] == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    }
    domain[i][0] = global_mesh->node_ID[2 * i + 1];
    n_domain[i]  = 1;
  }
  for (i = 0; i < global_mesh->n_subdomain; i++) {
    n_bnd = n_bnd_nlist[2 * i];
    n_out = n_bnd_nlist[2 * i + 1] - n_bnd_nlist[2 * i];
    if (n_out == 0) continue;
    for (j = 0; j < n_out; j++) {
      node                                 = bnd_nlist[i][n_bnd + j];
      domain[node - 1][n_domain[node - 1]] = i;
      n_domain[node - 1]++;
    }
  }

  /*make ngroup index list */
  ngrp_idx = (int **)HECMW_malloc(global_mesh->n_subdomain * sizeof(int *));
  if (ngrp_idx == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  for (i = 0; i < global_mesh->n_subdomain; i++) {
    ngrp_idx[i] =
        (int *)HECMW_calloc((node_group_global->n_grp + 1), sizeof(int));
    if (ngrp_idx[i] == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    }
  }
  for (i = 0; i < node_group_global->n_grp; i++) { /*skip group "ALL" */
    for (j = 0; j < global_mesh->n_subdomain; j++) {
      ngrp_idx[j][i + 1] = ngrp_idx[j][i];
    }
    if (node_group_global->grp_index[i + 1] - node_group_global->grp_index[i] ==
        global_mesh->n_node) {
      continue;
    }
    for (j = node_group_global->grp_index[i];
         j < node_group_global->grp_index[i + 1]; j++) {
      node = node_group_global->grp_item[j];
      for (k = 0; k < n_domain[node - 1]; k++) {
        current_domain = domain[node - 1][k];
        ngrp_idx[current_domain][i + 1]++;
      }
    }
  }

  /*make ngroup item list */
  ngrp_item = (int **)HECMW_malloc(global_mesh->n_subdomain * sizeof(int *));
  if (ngrp_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  for (i = 0; i < global_mesh->n_subdomain; i++) {
    ngrp_item[i] = (int *)HECMW_malloc(ngrp_idx[i][node_group_global->n_grp] *
                                       sizeof(int));
    if (ngrp_item[i] == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    }
    counter[i] = 0;
  }
  for (i = 0; i < node_group_global->n_grp; i++) { /*skip group "ALL" */
    if (node_group_global->grp_index[i + 1] - node_group_global->grp_index[i] ==
        global_mesh->n_node) {
      continue;
    }
    for (j = node_group_global->grp_index[i];
         j < node_group_global->grp_index[i + 1]; j++) {
      node = node_group_global->grp_item[j];
      for (k = 0; k < n_domain[node - 1]; k++) {
        current_domain = domain[node - 1][k];
        ngrp_item[current_domain][counter[current_domain]] = node;
        counter[current_domain]++;
      }
    }
  }

  for (i = 0; i < global_mesh->n_node; i++) {
    HECMW_free(domain[i]);
  }
  HECMW_free(n_domain);
  HECMW_free(domain);
  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int spdup_make_element_grouplist(
    const struct hecmwST_local_mesh *global_mesh) {
  struct hecmwST_elem_grp *elem_group_global = global_mesh->elem_group;
  int i, j, k, elem, n_bnd, n_out;
  int *n_domain = NULL;
  int **domain  = NULL;
  int current_domain;
  int counter[global_mesh->n_subdomain];

  /*make list of elem to domain(both internal and boundary) */
  n_domain = (int *)HECMW_calloc(global_mesh->n_elem, sizeof(int));
  if (n_domain == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  /*count outer elem(boundary and not internal) */
  for (i = 0; i < global_mesh->n_subdomain; i++) {
    n_bnd = n_bnd_elist[2 * i];
    n_out = n_bnd_elist[2 * i + 1] - n_bnd_elist[2 * i];
    if (n_out == 0) continue;
    for (j = 0; j < n_out; j++) {
      elem = bnd_elist[i][n_bnd + j];
      n_domain[elem - 1]++;
    }
  }
  /*make list */
  domain = (int **)HECMW_malloc(global_mesh->n_elem * sizeof(int *));
  if (domain == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  for (i = 0; i < global_mesh->n_elem; i++) {
    domain[i] = (int *)HECMW_malloc((n_domain[i] + 1) *
                                    sizeof(int)); /*+1 means internal elem */
    if (domain[i] == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    }
    domain[i][0] = global_mesh->elem_ID[2 * i + 1];
    n_domain[i]  = 1;
  }
  for (i = 0; i < global_mesh->n_subdomain; i++) {
    n_bnd = n_bnd_elist[2 * i];
    n_out = n_bnd_elist[2 * i + 1] - n_bnd_elist[2 * i];
    if (n_out == 0) continue;
    for (j = 0; j < n_out; j++) {
      elem                                 = bnd_elist[i][n_bnd + j];
      domain[elem - 1][n_domain[elem - 1]] = i;
      n_domain[elem - 1]++;
    }
  }

  /*make egroup index list */
  egrp_idx = (int **)HECMW_malloc(global_mesh->n_subdomain * sizeof(int *));
  if (egrp_idx == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  for (i = 0; i < global_mesh->n_subdomain; i++) {
    egrp_idx[i] =
        (int *)HECMW_calloc((elem_group_global->n_grp + 1), sizeof(int));
    if (egrp_idx[i] == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    }
  }
  for (i = 0; i < elem_group_global->n_grp; i++) { /*skip group "ALL" */
    for (j = 0; j < global_mesh->n_subdomain; j++) {
      egrp_idx[j][i + 1] = egrp_idx[j][i];
    }
    if (elem_group_global->grp_index[i + 1] - elem_group_global->grp_index[i] ==
        global_mesh->n_elem) {
      continue;
    }
    for (j = elem_group_global->grp_index[i];
         j < elem_group_global->grp_index[i + 1]; j++) {
      elem = elem_group_global->grp_item[j];
      for (k = 0; k < n_domain[elem - 1]; k++) {
        current_domain = domain[elem - 1][k];
        egrp_idx[current_domain][i + 1]++;
      }
    }
  }

  /*make egroup item list */
  egrp_item = (int **)HECMW_malloc(global_mesh->n_subdomain * sizeof(int *));
  if (egrp_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  for (i = 0; i < global_mesh->n_subdomain; i++) {
    egrp_item[i] = (int *)HECMW_malloc(egrp_idx[i][elem_group_global->n_grp] *
                                       sizeof(int));
    if (egrp_item[i] == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    }
    counter[i] = 0;
  }
  for (i = 0; i < elem_group_global->n_grp; i++) { /*skip group "ALL" */
    if (elem_group_global->grp_index[i + 1] - elem_group_global->grp_index[i] ==
        global_mesh->n_elem) {
      continue;
    }
    for (j = elem_group_global->grp_index[i];
         j < elem_group_global->grp_index[i + 1]; j++) {
      elem = elem_group_global->grp_item[j];
      for (k = 0; k < n_domain[elem - 1]; k++) {
        current_domain = domain[elem - 1][k];
        egrp_item[current_domain][counter[current_domain]] = elem;
        counter[current_domain]++;
      }
    }
  }

  for (i = 0; i < global_mesh->n_elem; i++) {
    HECMW_free(domain[i]);
  }
  HECMW_free(n_domain);
  HECMW_free(domain);
  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int spdup_makelist_main(const struct hecmwST_local_mesh *global_mesh) {
  int rtc;

  rtc = spdup_init_list(global_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = spdup_make_list(global_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = spdup_make_node_grouplist(global_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = spdup_make_element_grouplist(global_mesh);
  if (rtc != RTC_NORMAL) goto error;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static void spdup_freelist(const struct hecmwST_local_mesh *global_mesh) {
  int i;

  HECMW_free(n_int_nlist);
  HECMW_free(n_bnd_nlist);
  HECMW_free(n_int_elist);
  HECMW_free(n_bnd_elist);

  for (i = 0; i < global_mesh->n_subdomain; i++) {
    HECMW_free(int_nlist[i]);
    HECMW_free(bnd_nlist[i]);
    HECMW_free(int_elist[i]);
    HECMW_free(bnd_elist[i]);
    HECMW_free(ngrp_idx[i]);
    HECMW_free(ngrp_item[i]);
    HECMW_free(egrp_idx[i]);
    HECMW_free(egrp_item[i]);
  }

  HECMW_free(int_nlist);
  HECMW_free(bnd_nlist);
  HECMW_free(int_elist);
  HECMW_free(bnd_elist);
  HECMW_free(ngrp_idx);
  HECMW_free(ngrp_item);
  HECMW_free(egrp_idx);
  HECMW_free(egrp_item);
}

static int is_spdup_available(const struct hecmwST_local_mesh *global_mesh) {
  return global_mesh->hecmw_flag_parttype == HECMW_FLAG_PARTTYPE_NODEBASED &&
         global_mesh->hecmw_flag_partdepth == 1 &&
         global_mesh->mpc->n_mpc == 0 && global_mesh->contact_pair->n_pair == 0;
}

/*================================================================================================*/

static char *get_dist_file_name(char *header, int domain, char *fname) {
  char s_domain[HECMW_NAME_LEN + 1];

  sprintf(s_domain, "%d", domain);

  strcpy(fname, header);
  strcat(fname, ".");
  strcat(fname, s_domain);

  return fname;
}

static void free_link_list(struct link_unit *llist) {
  struct link_unit *p, *q;

  for (p = llist; p; p = q) {
    q = p->next;
    HECMW_free(p);
  }
  llist = NULL;
}

/*================================================================================================*/

static int init_struct_global(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'local_mesh\' is NULL");
    goto error;
  }

  memset(local_mesh->gridfile, 0, HECMW_NAME_LEN + 1);
  local_mesh->hecmw_n_file = 0;
  local_mesh->files        = NULL;
  memset(local_mesh->header, 0, HECMW_HEADER_LEN + 1);

  local_mesh->hecmw_flag_adapt       = 0;
  local_mesh->hecmw_flag_initcon     = 0;
  local_mesh->hecmw_flag_parttype    = 0;
  local_mesh->hecmw_flag_partdepth   = 0;
  local_mesh->hecmw_flag_version     = 0;
  local_mesh->hecmw_flag_partcontact = 0;

  local_mesh->zero_temp = 0.0;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int init_struct_node(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'local_mesh\' is NULL");
    goto error;
  }

  local_mesh->n_node             = 0;
  local_mesh->n_node_gross       = 0;
  local_mesh->nn_internal        = 0;
  local_mesh->node_internal_list = NULL;

  local_mesh->node           = NULL;
  local_mesh->node_ID        = NULL;
  local_mesh->global_node_ID = NULL;

  local_mesh->n_dof          = 0;
  local_mesh->n_dof_grp      = 0;
  local_mesh->node_dof_index = NULL;
  local_mesh->node_dof_item  = NULL;

  local_mesh->node_val_index = NULL;
  local_mesh->node_val_item  = NULL;

  local_mesh->node_init_val_index = NULL;
  local_mesh->node_init_val_item  = NULL;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int init_struct_elem(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'local_mesh\' is NULL");
    goto error;
  }

  local_mesh->n_elem             = 0;
  local_mesh->n_elem_gross       = 0;
  local_mesh->ne_internal        = 0;
  local_mesh->elem_internal_list = NULL;

  local_mesh->elem_ID        = NULL;
  local_mesh->global_elem_ID = NULL;

  local_mesh->n_elem_type     = 0;
  local_mesh->elem_type       = NULL;
  local_mesh->elem_type_index = NULL;
  local_mesh->elem_type_item  = NULL;

  local_mesh->elem_node_index = NULL;
  local_mesh->elem_node_item  = NULL;

  local_mesh->section_ID = NULL;

  local_mesh->n_elem_mat_ID     = 0;
  local_mesh->elem_mat_ID_index = NULL;
  local_mesh->elem_mat_ID_item  = NULL;

  local_mesh->elem_mat_int_index = NULL;
  local_mesh->elem_mat_int_val   = NULL;

  local_mesh->elem_val_index = NULL;
  local_mesh->elem_val_item  = NULL;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int init_struct_comm(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'local_mesh\' is NULL");
    goto error;
  }

  local_mesh->zero        = 0;
  local_mesh->PETOT       = 0;
  local_mesh->PEsmpTOT    = 0;
  local_mesh->my_rank     = 0;
  local_mesh->errnof      = 0;
  local_mesh->n_subdomain = 0;

  local_mesh->n_neighbor_pe = 0;
  local_mesh->neighbor_pe   = NULL;

  local_mesh->import_index = NULL;
  local_mesh->import_item  = NULL;
  local_mesh->export_index = NULL;
  local_mesh->export_item  = NULL;
  local_mesh->shared_index = NULL;
  local_mesh->shared_item  = NULL;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int init_struct_adapt(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'local_mesh\' is NULL");
    goto error;
  }

  local_mesh->coarse_grid_level       = 0;
  local_mesh->n_adapt                 = 0;
  local_mesh->when_i_was_refined_node = NULL;
  local_mesh->when_i_was_refined_elem = NULL;
  local_mesh->adapt_parent_type       = NULL;
  local_mesh->adapt_type              = NULL;
  local_mesh->adapt_level             = NULL;
  local_mesh->adapt_parent            = NULL;
  local_mesh->adapt_children_index    = NULL;
  local_mesh->adapt_children_item     = NULL;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int init_struct_sect(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'local_mesh\' is NULL");
    goto error;
  }
  if (local_mesh->section == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'local_mesh->section\' is NULL");
    goto error;
  }

  local_mesh->section->n_sect            = 0;
  local_mesh->section->sect_type         = NULL;
  local_mesh->section->sect_opt          = NULL;
  local_mesh->section->sect_mat_ID_index = NULL;
  local_mesh->section->sect_mat_ID_item  = NULL;
  local_mesh->section->sect_I_index      = NULL;
  local_mesh->section->sect_I_item       = NULL;
  local_mesh->section->sect_R_index      = NULL;
  local_mesh->section->sect_R_item       = NULL;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int init_struct_mat(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'local_mesh\' is NULL");
    goto error;
  }
  if (local_mesh->material == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'local_mesh->material\' is NULL");
    goto error;
  }

  local_mesh->material->n_mat             = 0;
  local_mesh->material->n_mat_item        = 0;
  local_mesh->material->n_mat_subitem     = 0;
  local_mesh->material->n_mat_table       = 0;
  local_mesh->material->mat_name          = NULL;
  local_mesh->material->mat_item_index    = NULL;
  local_mesh->material->mat_subitem_index = NULL;
  local_mesh->material->mat_table_index   = NULL;
  local_mesh->material->mat_val           = NULL;
  local_mesh->material->mat_temp          = NULL;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int init_struct_mpc(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'local_mesh\' is NULL");
    return -1;
  }
  if (local_mesh->mpc == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'local_mesh->mpc\' is NULL");
    goto error;
  }

  local_mesh->mpc->n_mpc     = 0;
  local_mesh->mpc->mpc_index = NULL;
  local_mesh->mpc->mpc_item  = NULL;
  local_mesh->mpc->mpc_dof   = NULL;
  local_mesh->mpc->mpc_val   = NULL;
  local_mesh->mpc->mpc_const = NULL;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int init_struct_amp(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'local_mesh\' is NULL");
    goto error;
  }

  if (local_mesh->amp == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'local_mesh->amp\' is NULL");
    goto error;
  }

  local_mesh->amp->n_amp               = 0;
  local_mesh->amp->amp_name            = NULL;
  local_mesh->amp->amp_type_definition = NULL;
  local_mesh->amp->amp_type_time       = NULL;
  local_mesh->amp->amp_type_value      = NULL;
  local_mesh->amp->amp_index           = NULL;
  local_mesh->amp->amp_val             = NULL;
  local_mesh->amp->amp_table           = NULL;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int init_struct_node_grp(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'local_mesh\' is NULL");
    goto error;
  }

  if (local_mesh->node_group == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'local_mesh->node_group\' is NULL");
    goto error;
  }

  local_mesh->node_group->n_grp     = 0;
  local_mesh->node_group->grp_name  = NULL;
  local_mesh->node_group->grp_index = NULL;
  local_mesh->node_group->grp_item  = NULL;

  local_mesh->node_group->n_bc         = 0;
  local_mesh->node_group->bc_grp_ID    = 0;
  local_mesh->node_group->bc_grp_type  = 0;
  local_mesh->node_group->bc_grp_index = 0;
  local_mesh->node_group->bc_grp_dof   = 0;
  local_mesh->node_group->bc_grp_val   = 0;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int init_struct_elem_grp(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'local_mesh\' is NULL");
    goto error;
  }

  if (local_mesh->elem_group == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'local_mesh->elem_group\' is NULL");
    goto error;
  }

  local_mesh->elem_group->n_grp     = 0;
  local_mesh->elem_group->grp_name  = NULL;
  local_mesh->elem_group->grp_index = NULL;
  local_mesh->elem_group->grp_item  = NULL;

  local_mesh->elem_group->n_bc         = 0;
  local_mesh->elem_group->bc_grp_ID    = NULL;
  local_mesh->elem_group->bc_grp_type  = NULL;
  local_mesh->elem_group->bc_grp_index = NULL;
  local_mesh->elem_group->bc_grp_val   = NULL;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int init_struct_surf_grp(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'local_mesh\' is NULL");
    goto error;
  }

  if (local_mesh->surf_group == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'local_mesh->surf_group\' is NULL");
    goto error;
  }

  local_mesh->surf_group->n_grp     = 0;
  local_mesh->surf_group->grp_name  = NULL;
  local_mesh->surf_group->grp_index = NULL;
  local_mesh->surf_group->grp_item  = NULL;

  local_mesh->surf_group->n_bc         = 0;
  local_mesh->surf_group->bc_grp_ID    = NULL;
  local_mesh->surf_group->bc_grp_type  = NULL;
  local_mesh->surf_group->bc_grp_index = NULL;
  local_mesh->surf_group->bc_grp_val   = NULL;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int init_struct_contact_pair(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'local_mesh\' is NULL");
    goto error;
  }

  if (local_mesh->contact_pair == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG,
                    "\'local_mesh->contact_pair\' is NULL");
    goto error;
  }

  local_mesh->contact_pair->n_pair        = 0;
  local_mesh->contact_pair->name          = NULL;
  local_mesh->contact_pair->type          = NULL;
  local_mesh->contact_pair->slave_grp_id  = NULL;
  local_mesh->contact_pair->slave_orisgrp_id  = NULL;
  local_mesh->contact_pair->master_grp_id = NULL;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*================================================================================================*/

static void clean_struct_global(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) return;

  init_struct_global(local_mesh);
}

static void clean_struct_node(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) return;

  if (local_mesh->node_internal_list) {
    HECMW_free(local_mesh->node_internal_list);
  }
  if (local_mesh->node) {
    HECMW_free(local_mesh->node);
  }
  if (local_mesh->node_ID) {
    HECMW_free(local_mesh->node_ID);
  }
  if (local_mesh->global_node_ID) {
    HECMW_free(local_mesh->global_node_ID);
  }
  if (local_mesh->node_dof_index) {
    HECMW_free(local_mesh->node_dof_index);
  }
  if (local_mesh->node_init_val_index) {
    HECMW_free(local_mesh->node_init_val_index);
  }
  if (local_mesh->node_init_val_item) {
    HECMW_free(local_mesh->node_init_val_item);
  }

  init_struct_node(local_mesh);
}

static void clean_struct_elem(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) return;

  if (local_mesh->elem_internal_list) {
    HECMW_free(local_mesh->elem_internal_list);
  }
  if (local_mesh->elem_ID) {
    HECMW_free(local_mesh->elem_ID);
  }
  if (local_mesh->global_elem_ID) {
    HECMW_free(local_mesh->global_elem_ID);
  }
  if (local_mesh->elem_type) {
    HECMW_free(local_mesh->elem_type);
  }
  if (local_mesh->elem_type_index) {
    HECMW_free(local_mesh->elem_type_index);
  }
  if (local_mesh->elem_node_index) {
    HECMW_free(local_mesh->elem_node_index);
  }
  if (local_mesh->elem_node_item) {
    HECMW_free(local_mesh->elem_node_item);
  }
  if (local_mesh->section_ID) {
    HECMW_free(local_mesh->section_ID);
  }
  if (local_mesh->elem_mat_ID_index) {
    HECMW_free(local_mesh->elem_mat_ID_index);
  }
  if (local_mesh->elem_mat_ID_item) {
    HECMW_free(local_mesh->elem_mat_ID_item);
  }

  init_struct_elem(local_mesh);
}

static void clean_struct_comm(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) return;

  if (local_mesh->neighbor_pe) {
    HECMW_free(local_mesh->neighbor_pe);
  }
  if (local_mesh->import_index) {
    HECMW_free(local_mesh->import_index);
  }
  if (local_mesh->import_item) {
    HECMW_free(local_mesh->import_item);
  }
  if (local_mesh->export_index) {
    HECMW_free(local_mesh->export_index);
  }
  if (local_mesh->export_item) {
    HECMW_free(local_mesh->export_item);
  }
  if (local_mesh->shared_index) {
    HECMW_free(local_mesh->shared_index);
  }
  if (local_mesh->shared_item) {
    HECMW_free(local_mesh->shared_item);
  }

  init_struct_comm(local_mesh);
}

static void clean_struct_adapt(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) return;

  init_struct_adapt(local_mesh);
}

static void clean_struct_sect(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) return;
  if (local_mesh->section == NULL) return;

  init_struct_sect(local_mesh);
}

static void clean_struct_mat(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) return;
  if (local_mesh->material == NULL) return;

  init_struct_mat(local_mesh);
}

static void clean_struct_mpc(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) return;
  if (local_mesh->mpc == NULL) return;

  HECMW_free(local_mesh->mpc->mpc_index);
  HECMW_free(local_mesh->mpc->mpc_item);
  HECMW_free(local_mesh->mpc->mpc_dof);
  HECMW_free(local_mesh->mpc->mpc_val);
  HECMW_free(local_mesh->mpc->mpc_const);

  init_struct_mpc(local_mesh);
}

static void clean_struct_amp(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) return;
  if (local_mesh->amp == NULL) return;

  init_struct_amp(local_mesh);
}

static void clean_struct_node_grp(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) return;
  if (local_mesh->node_group == NULL) return;

  if (local_mesh->node_group->grp_index) {
    HECMW_free(local_mesh->node_group->grp_index);
  }
  if (local_mesh->node_group->grp_item) {
    HECMW_free(local_mesh->node_group->grp_item);
  }

  init_struct_node_grp(local_mesh);
}

static void clean_struct_elem_grp(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) return;
  if (local_mesh->elem_group == NULL) return;

  if (local_mesh->elem_group->grp_index) {
    HECMW_free(local_mesh->elem_group->grp_index);
  }
  if (local_mesh->elem_group->grp_item) {
    HECMW_free(local_mesh->elem_group->grp_item);
  }

  init_struct_elem_grp(local_mesh);
}

static void clean_struct_surf_grp(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) return;
  if (local_mesh->surf_group == NULL) return;

  if (local_mesh->surf_group->grp_index) {
    HECMW_free(local_mesh->surf_group->grp_index);
  }
  if (local_mesh->surf_group->grp_item) {
    HECMW_free(local_mesh->surf_group->grp_item);
  }

  init_struct_surf_grp(local_mesh);
}

static void clean_struct_contact_pair(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) return;
  if (local_mesh->contact_pair == NULL) return;

  if (local_mesh->contact_pair->type) {
    HECMW_free(local_mesh->contact_pair->type);
  }
  if (local_mesh->contact_pair->slave_grp_id) {
    HECMW_free(local_mesh->contact_pair->slave_grp_id);
  }
  if (local_mesh->contact_pair->slave_orisgrp_id) {
    HECMW_free(local_mesh->contact_pair->slave_orisgrp_id);
  }
  if (local_mesh->contact_pair->master_grp_id) {
    HECMW_free(local_mesh->contact_pair->master_grp_id);
  }

  init_struct_contact_pair(local_mesh);
}

static void clean_struct_local_mesh(struct hecmwST_local_mesh *local_mesh) {
  if (local_mesh == NULL) return;

  clean_struct_global(local_mesh);
  clean_struct_node(local_mesh);
  clean_struct_elem(local_mesh);
  clean_struct_comm(local_mesh);
  clean_struct_adapt(local_mesh);
  clean_struct_sect(local_mesh);
  clean_struct_mat(local_mesh);
  clean_struct_mpc(local_mesh);
  clean_struct_amp(local_mesh);
  clean_struct_node_grp(local_mesh);
  clean_struct_elem_grp(local_mesh);
  clean_struct_surf_grp(local_mesh);
  clean_struct_contact_pair(local_mesh);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * - - - - - - - - - */

static int init_struct_result_data(struct hecmwST_result_data *result_data) {
  if (result_data == NULL) {
    HECMW_set_error(errno, "\'result_data\' is NULL");
    goto error;
  }

  result_data->nn_dof        = NULL;
  result_data->node_label    = NULL;
  result_data->node_val_item = NULL;

  result_data->ne_dof        = NULL;
  result_data->elem_label    = NULL;
  result_data->elem_val_item = NULL;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static void free_struct_result_data(struct hecmwST_result_data *result_data) {
  int i;

  if (result_data == NULL) return;

  HECMW_free(result_data->nn_dof);
  HECMW_free(result_data->ne_dof);

  if (result_data->node_label) {
    for (i = 0; i < result_data->nn_component; i++) {
      HECMW_free(result_data->node_label[i]);
    }
    HECMW_free(result_data->node_label);
  }
  if (result_data->elem_label) {
    for (i = 0; i < result_data->ne_component; i++) {
      HECMW_free(result_data->elem_label[i]);
    }
    HECMW_free(result_data->elem_label);
  }

  HECMW_free(result_data->node_val_item);
  HECMW_free(result_data->elem_val_item);

  HECMW_free(result_data);
  result_data = NULL;
}

/*================================================================================================*/

static int search_eqn_block_idx(const struct hecmwST_local_mesh *mesh) {
  int i;

  for (i = 0; i < mesh->node_group->n_grp; i++) {
    if (!strcmp(mesh->node_group->grp_name[i], HECMW_PART_EQUATION_BLOCK_NAME))
      return i;
  }

  return -1;
}

/*================================================================================================*/

static int quick_sort(int no, int n, double *arr, int *brr, int *istack) {
  double a, atemp;
  int b, btemp;
  int i, ir, j, k, l;
  int jstack = 0;
  int nstack;

  nstack = no;
  l      = 0;
  ir     = n - 1;

  for (;;) {
    if (ir - l < QSORT_LOWER) {
      for (j = l + 1; j <= ir; j++) {
        a = arr[j];
        b = brr[j];
        for (i = j - 1; i >= l; i--) {
          if (arr[i] <= a) break;
          arr[i + 1] = arr[i];
          brr[i + 1] = brr[i];
        }
        arr[i + 1] = a;
        brr[i + 1] = b;
      }

      if (!jstack) return 0;

      ir = istack[jstack];
      l  = istack[jstack - 1];
      jstack -= 2;

    } else {
      k = (l + ir) >> 1;

      DSWAP(arr[k], arr[l + 1])
      ISWAP(brr[k], brr[l + 1])

      if (arr[l] > arr[ir]) {
        DSWAP(arr[l], arr[ir])
        ISWAP(brr[l], brr[ir])
      }

      if (arr[l + 1] > arr[ir]) {
        DSWAP(arr[l + 1], arr[ir])
        ISWAP(brr[l + 1], brr[ir])
      }

      if (arr[l] > arr[l + 1]) {
        DSWAP(arr[l], arr[l + 1])
        ISWAP(brr[l], brr[l + 1])
      }

      i = l + 1;
      j = ir;
      a = arr[l + 1];
      b = brr[l + 1];

      for (;;) {
        do
          i++;
        while (arr[i] < a);
        do
          j--;
        while (arr[j] > a);

        if (j < i) break;

        DSWAP(arr[i], arr[j])
        ISWAP(brr[i], brr[j])
      }

      arr[l + 1] = arr[j];
      arr[j]     = a;
      brr[l + 1] = brr[j];
      brr[j]     = b;

      jstack += 2;

      if (jstack > nstack) {
        HECMW_set_error(HECMW_PART_E_STACK_OVERFLOW, "");
        return -1;
      }

      if (ir - i + 1 >= j - l) {
        istack[jstack]     = ir;
        istack[jstack - 1] = i;
        ir                 = j - 1;
      } else {
        istack[jstack]     = j - 1;
        istack[jstack - 1] = l;
        l                  = i;
      }
    }
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * - - - - - - - - - */

static int rcb_partition(int n, const double *coord, int *wnum,
                         const struct hecmw_part_cont_data *cont_data) {
  double *value;
  int *id, *stack;
  int rtc;
  int counter;
  int i, j, k;

  id = (int *)HECMW_malloc(sizeof(int) * n);
  if (id == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  stack = (int *)HECMW_malloc(sizeof(int) * n);
  if (stack == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  value = (double *)HECMW_malloc(sizeof(double) * n);
  if (value == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (i = 0; i < cont_data->n_rcb_div; i++) {
    for (j = 0; j < pow(2, i); j++) {
      counter = 0;

      switch (cont_data->rcb_axis[i]) {
        case HECMW_PART_RCB_X_AXIS: /* X-axis */
          for (k = 0; k < n; k++) {
            if (wnum[2 * k + 1] == j) {
              id[counter]    = k;
              value[counter] = coord[3 * k];
              counter++;
            }
          }
          break;

        case HECMW_PART_RCB_Y_AXIS: /* Y-axis */
          for (k = 0; k < n; k++) {
            if (wnum[2 * k + 1] == j) {
              id[counter]    = k;
              value[counter] = coord[3 * k + 1];
              counter++;
            }
          }
          break;

        case HECMW_PART_RCB_Z_AXIS: /* Z-axis */
          for (k = 0; k < n; k++) {
            if (wnum[2 * k + 1] == j) {
              id[counter]    = k;
              value[counter] = coord[3 * k + 2];
              counter++;
            }
          }
          break;

        default:
          HECMW_set_error(HECMW_PART_E_INVALID_RCB_DIR, "");
          goto error;
      }

      /* quick sort */
      rtc = quick_sort(n, counter, value, id, stack);
      if (rtc != RTC_NORMAL) goto error;

      /* belonging domain of node */
      for (k = 0; k < counter * F_1_2; k++) {
        wnum[2 * id[k] + 1] = j + (int)pow(2, i);
      }
    }
  }

  HECMW_free(id);
  HECMW_free(stack);
  HECMW_free(value);

  return RTC_NORMAL;

error:
  HECMW_free(id);
  HECMW_free(stack);
  HECMW_free(value);

  return RTC_ERROR;
}

/*------------------------------------------------------------------------------------------------*/

static int calc_gravity(const struct hecmwST_local_mesh *global_mesh,
                        double *coord) {
  double coord_x, coord_y, coord_z;
  int node;
  int js, je;
  int i, j;

  for (i = 0; i < global_mesh->n_elem; i++) {
    js = global_mesh->elem_node_index[i];
    je = global_mesh->elem_node_index[i + 1];

    for (coord_x = 0.0, coord_y = 0.0, coord_z = 0.0, j = js; j < je; j++) {
      node = global_mesh->elem_node_item[j];

      coord_x += global_mesh->node[3 * (node - 1)];
      coord_y += global_mesh->node[3 * (node - 1) + 1];
      coord_z += global_mesh->node[3 * (node - 1) + 2];
    }

    coord[3 * i]     = coord_x / (je - js);
    coord[3 * i + 1] = coord_y / (je - js);
    coord[3 * i + 2] = coord_z / (je - js);
  }

  return RTC_NORMAL;
}

static int rcb_partition_eb(struct hecmwST_local_mesh *global_mesh,
                            const struct hecmw_part_cont_data *cont_data) {
  double *coord = NULL;
  int rtc;

  coord = (double *)HECMW_malloc(sizeof(double) * global_mesh->n_elem * 3);
  if (coord == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  rtc = calc_gravity(global_mesh, coord);
  if (rtc != RTC_NORMAL) goto error;

  rtc = rcb_partition(global_mesh->n_elem, coord, global_mesh->elem_ID,
                      cont_data);
  if (rtc != RTC_NORMAL) goto error;

  HECMW_free(coord);

  return RTC_NORMAL;

error:
  HECMW_free(coord);

  return RTC_ERROR;
}

/*================================================================================================*/

static int create_node_graph_link_list(
    const struct hecmwST_local_mesh *global_mesh,
    const struct hecmw_part_edge_data *edge_data, struct link_list **graph) {
  int node1, node2;
  long long int i;

  for (i = 0; i < edge_data->n_edge; i++) {
    node1 = edge_data->edge_node_item[2 * i];
    node2 = edge_data->edge_node_item[2 * i + 1];

    /* node 1 */
    graph[node1 - 1]->last->next =
        (struct link_unit *)HECMW_malloc(sizeof(struct link_unit));
    if (graph[node1 - 1]->last->next == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    }

    graph[node1 - 1]->n += 1;
    graph[node1 - 1]->last->next->id   = node2;
    graph[node1 - 1]->last->next->next = NULL;
    graph[node1 - 1]->last             = graph[node1 - 1]->last->next;

    /* node 2 */
    graph[node2 - 1]->last->next =
        (struct link_unit *)HECMW_malloc(sizeof(struct link_unit));
    if (graph[node2 - 1]->last->next == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    }

    graph[node2 - 1]->n += 1;
    graph[node2 - 1]->last->next->id   = node1;
    graph[node2 - 1]->last->next->next = NULL;
    graph[node2 - 1]->last             = graph[node2 - 1]->last->next;
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int create_node_graph_compress(
    const struct hecmwST_local_mesh *global_mesh, struct link_list **graph,
    int *node_graph_index, int *node_graph_item) {
  int counter;
  int i, j;
  struct link_unit *p;

  for (counter = 0, i = 0; i < global_mesh->n_node; i++) {
    node_graph_index[i + 1] = node_graph_index[i] + graph[i]->n;

    for (p = graph[i]->list, j = 0; j < graph[i]->n; j++) {
      p                          = p->next;
      node_graph_item[counter++] = p->id - 1;
    }
  }

  return RTC_NORMAL;
}

static int create_node_graph(const struct hecmwST_local_mesh *global_mesh,
                             const struct hecmw_part_edge_data *edge_data,
                             int *node_graph_index, int *node_graph_item) {
  struct link_list **graph = NULL;
  int rtc;
  int i;

  graph = (struct link_list **)HECMW_malloc(sizeof(struct link_list *) *
                                            global_mesh->n_node);
  if (graph == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  } else {
    for (i = 0; i < global_mesh->n_node; i++) {
      graph[i] = NULL;
    }
  }
  for (i = 0; i < global_mesh->n_node; i++) {
    graph[i] = (struct link_list *)HECMW_malloc(sizeof(struct link_list));
    if (graph[i] == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    } else {
      graph[i]->list = NULL;
    }
  }
  for (i = 0; i < global_mesh->n_node; i++) {
    graph[i]->list = (struct link_unit *)HECMW_malloc(sizeof(struct link_unit));
    if (graph[i]->list == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    } else {
      graph[i]->n          = 0;
      graph[i]->list->next = NULL;
      graph[i]->last       = graph[i]->list;
    }
  }

  rtc = create_node_graph_link_list(global_mesh, edge_data, graph);
  if (rtc != RTC_NORMAL) goto error;

  rtc = create_node_graph_compress(global_mesh, graph, node_graph_index,
                                   node_graph_item);
  if (rtc != RTC_NORMAL) goto error;

  for (i = 0; i < global_mesh->n_node; i++) {
    free_link_list(graph[i]->list);
    HECMW_free(graph[i]);
  }
  HECMW_free(graph);

  return RTC_NORMAL;

error:
  if (graph) {
    for (i = 0; i < global_mesh->n_node; i++) {
      if (graph[i]) {
        free_link_list(graph[i]->list);
        HECMW_free(graph[i]);
      }
    }
    HECMW_free(graph);
  }

  return RTC_ERROR;
}

/*------------------------------------------------------------------------------------------------*/

static int set_node_belong_elem(const struct hecmwST_local_mesh *global_mesh,
                                struct hecmw_part_node_data *node_data) {
  int node, counter;
  struct link_list **node_list = NULL;
  struct link_unit *p;
  int size;
  int i, j;

  node_data->node_elem_index = NULL;
  node_data->node_elem_item  = NULL;

  node_list = (struct link_list **)HECMW_malloc(sizeof(struct link_list *) *
                                                global_mesh->n_node);
  if (node_list == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  } else {
    for (i = 0; i < global_mesh->n_node; i++) {
      node_list[i] = NULL;
    }
  }
  for (i = 0; i < global_mesh->n_node; i++) {
    node_list[i] = (struct link_list *)HECMW_malloc(sizeof(struct link_list));
    if (node_list[i] == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    } else {
      node_list[i]->list = NULL;
    }
  }
  for (i = 0; i < global_mesh->n_node; i++) {
    node_list[i]->list =
        (struct link_unit *)HECMW_malloc(sizeof(struct link_unit));
    if (node_list[i]->list == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    } else {
      node_list[i]->n          = 0;
      node_list[i]->list->next = NULL;
      node_list[i]->last       = node_list[i]->list;
    }
  }

  for (i = 0; i < global_mesh->n_elem; i++) {
    for (j = global_mesh->elem_node_index[i];
         j < global_mesh->elem_node_index[i + 1]; j++) {
      node = global_mesh->elem_node_item[j];

      size                            = sizeof(struct link_list);
      node_list[node - 1]->last->next = (struct link_unit *)HECMW_malloc(size);
      if (node_list[node - 1]->last->next == NULL) {
        HECMW_set_error(errno, "");
        goto error;
      }

      node_list[node - 1]->last       = node_list[node - 1]->last->next;
      node_list[node - 1]->last->id   = i + 1;
      node_list[node - 1]->last->next = NULL;
      node_list[node - 1]->n += 1;
    }
  }

  node_data->node_elem_index =
      (int *)HECMW_calloc(global_mesh->n_node + 1, sizeof(int));
  if (node_data->node_elem_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  for (i = 0; i < global_mesh->n_node; i++) {
    node_data->node_elem_index[i + 1] =
        node_data->node_elem_index[i] + node_list[i]->n;
  }

  size = sizeof(int) * node_data->node_elem_index[global_mesh->n_node];
  node_data->node_elem_item = (int *)HECMW_malloc(size);
  if (node_data->node_elem_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  for (counter = 0, i = 0; i < global_mesh->n_node; i++) {
    for (p = node_list[i]->list, j = 0; j < node_list[i]->n; j++) {
      p                                    = p->next;
      node_data->node_elem_item[counter++] = p->id;
    }
    HECMW_assert(counter == node_data->node_elem_index[i + 1]);
  }

  for (i = 0; i < global_mesh->n_node; i++) {
    free_link_list(node_list[i]->list);
    HECMW_free(node_list[i]);
  }
  HECMW_free(node_list);

  return RTC_NORMAL;

error:
  if (node_list) {
    for (i = 0; i < global_mesh->n_node; i++) {
      if (node_list[i]) {
        free_link_list(node_list[i]->list);
        HECMW_free(node_list[i]);
      }
    }
    HECMW_free(node_list);
  }

  HECMW_free(node_data->node_elem_index);
  HECMW_free(node_data->node_elem_item);
  node_data->node_elem_index = NULL;
  node_data->node_elem_item  = NULL;

  return RTC_ERROR;
}

static int create_elem_graph_link_list(
    const struct hecmwST_local_mesh *global_mesh,
    const struct hecmw_part_node_data *node_data, struct link_list **graph) {
  char *elem_flag = NULL;
  int elem, node;
  int size;
  int counter;
  int i, j, k;

  elem_flag = (char *)HECMW_malloc(sizeof(char) * global_mesh->n_elem);
  if (elem_flag == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (counter = 0, i = 0; i < global_mesh->n_elem; i++) {
    memset(elem_flag, 0, sizeof(char) * global_mesh->n_elem);
    MASK_BIT(elem_flag[i], MASK);

    for (j = global_mesh->elem_node_index[i];
         j < global_mesh->elem_node_index[i + 1]; j++) {
      node = global_mesh->elem_node_item[j];

      for (k = node_data->node_elem_index[node - 1];
           k < node_data->node_elem_index[node]; k++) {
        elem = node_data->node_elem_item[k];

        if (!EVAL_BIT(elem_flag[elem - 1], MASK)) {
          MASK_BIT(elem_flag[elem - 1], MASK);

          size                 = sizeof(struct link_unit);
          graph[i]->last->next = (struct link_unit *)HECMW_malloc(size);
          if (graph[i]->last->next == NULL) {
            HECMW_set_error(errno, "");
            goto error;
          }

          graph[i]->n += 1;
          graph[i]->last->next->id   = elem;
          graph[i]->last->next->next = NULL;
          graph[i]->last             = graph[i]->last->next;
          counter++;
        }
      }
    }
  }

  HECMW_free(elem_flag);

  return counter;

error:
  HECMW_free(elem_flag);

  return -1;
}

static int create_elem_graph_compress(
    const struct hecmwST_local_mesh *global_mesh, struct link_list **graph,
    int *elem_graph_index, int *elem_graph_item) {
  struct link_unit *p;
  int counter;
  int i, j;

  for (counter = 0, i = 0; i < global_mesh->n_elem; i++) {
    elem_graph_index[i + 1] = elem_graph_index[i] + graph[i]->n;

    for (p = graph[i]->list, j = 0; j < graph[i]->n; j++) {
      p                          = p->next;
      elem_graph_item[counter++] = p->id - 1;
    }
  }
  HECMW_assert(elem_graph_index[global_mesh->n_elem] == counter);

  return RTC_NORMAL;
}

static int *create_elem_graph(const struct hecmwST_local_mesh *global_mesh,
                              int *elem_graph_index) {
  struct hecmw_part_node_data *node_data = NULL;
  struct link_list **graph               = NULL;
  int *elem_graph_item                   = NULL;
  int n_graph;
  int rtc;
  int i;

  node_data = (struct hecmw_part_node_data *)HECMW_malloc(
      sizeof(struct hecmw_part_node_data));
  if (node_data == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  } else {
    node_data->node_elem_index = NULL;
    node_data->node_elem_item  = NULL;
  }

  rtc = set_node_belong_elem(global_mesh, node_data);
  if (rtc != RTC_NORMAL) goto error;

  graph = (struct link_list **)HECMW_malloc(sizeof(struct link_list *) *
                                            global_mesh->n_elem);
  if (graph == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  } else {
    for (i = 0; i < global_mesh->n_elem; i++) {
      graph[i] = NULL;
    }
  }
  for (i = 0; i < global_mesh->n_elem; i++) {
    graph[i] = (struct link_list *)HECMW_malloc(sizeof(struct link_list));
    if (graph[i] == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    } else {
      graph[i]->list = NULL;
    }
  }
  for (i = 0; i < global_mesh->n_elem; i++) {
    graph[i]->list = (struct link_unit *)HECMW_malloc(sizeof(struct link_unit));
    if (graph[i]->list == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    } else {
      graph[i]->n          = 0;
      graph[i]->list->next = NULL;
      graph[i]->last       = graph[i]->list;
    }
  }

  n_graph = create_elem_graph_link_list(global_mesh, node_data, graph);
  if (n_graph < 0) goto error;

  elem_graph_item = (int *)HECMW_malloc(sizeof(int) * n_graph);
  if (elem_graph_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  rtc = create_elem_graph_compress(global_mesh, graph, elem_graph_index,
                                   elem_graph_item);
  if (rtc != RTC_NORMAL) goto error;

  HECMW_free(node_data->node_elem_index);
  HECMW_free(node_data->node_elem_item);
  HECMW_free(node_data);
  for (i = 0; i < global_mesh->n_elem; i++) {
    free_link_list(graph[i]->list);
    HECMW_free(graph[i]);
  }
  HECMW_free(graph);

  return elem_graph_item;

error:
  if (node_data) {
    HECMW_free(node_data->node_elem_index);
    HECMW_free(node_data->node_elem_item);
    HECMW_free(node_data);
  }
  if (graph) {
    for (i = 0; i < global_mesh->n_elem; i++) {
      if (graph[i]) {
        free_link_list(graph[i]->list);
        HECMW_free(graph[i]);
      }
    }
    HECMW_free(graph);
  }
  HECMW_free(elem_graph_item);

  return NULL;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * - - - - - - - - - */

static int pmetis_interface(const int n_vertex, const int n_domain, int *xadj,
                            int *adjncy, int *part) {
  int edgecut = 0; /* number of edge-cut */
#ifdef HECMW_PART_WITH_METIS
  int n       = n_vertex; /* number of vertices */
  int *vwgt   = NULL;     /* weight for vertices */
  int *adjwgt = NULL;     /* weight for edges */
  int nparts  = n_domain; /* number of sub-domains */

#if defined(METIS_VER_MAJOR) && (METIS_VER_MAJOR == 5)
  int ncon       = 1; /* number of balancing constraints */
  int *vsize     = NULL;
  real_t *tpwgts = NULL;
  real_t *ubvec  = NULL;
  int *options   = NULL;

  HECMW_log(HECMW_LOG_DEBUG, "Entering pmetis(v5)...\n");
  METIS_PartGraphRecursive(&n, &ncon, xadj, adjncy, vwgt, vsize, adjwgt,
                           &nparts, tpwgts, ubvec, options, &edgecut, part);
  HECMW_log(HECMW_LOG_DEBUG, "Returned from pmetis(v5)\n");
#else
  int wgtflag    = 0;               /* flag of weight for edges */
  int numflag    = 0;               /* flag of stating number of index */
  int options[5] = {0, 0, 0, 0, 0}; /* options for pMETIS */

  HECMW_log(HECMW_LOG_DEBUG, "Entering pmetis(v4)...\n");
  METIS_PartGraphRecursive(&n, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                           &nparts, options, &edgecut, part);
  HECMW_log(HECMW_LOG_DEBUG, "Returned from pmetis(v4)\n");
#endif
#endif

  return edgecut;
}

static int kmetis_interface(const int n_vertex, const int n_domain, int *xadj,
                            int *adjncy, int *part) {
  int edgecut = 0; /* number of edge-cut */
#ifdef HECMW_PART_WITH_METIS
  int n       = n_vertex; /* number of vertices */
  int *vwgt   = NULL;     /* weight for vertices */
  int *adjwgt = NULL;     /* weight for edges */
  int nparts  = n_domain; /* number of sub-domains */

#if defined(METIS_VER_MAJOR) && (METIS_VER_MAJOR == 5)
  int ncon       = 1; /* number of balancing constraints */
  int *vsize     = NULL;
  real_t *tpwgts = NULL;
  real_t *ubvec  = NULL;
  int *options   = NULL;

  HECMW_log(HECMW_LOG_DEBUG, "Entering kmetis(v5)...\n");
  METIS_PartGraphKway(&n, &ncon, xadj, adjncy, vwgt, vsize, adjwgt, &nparts,
                      tpwgts, ubvec, options, &edgecut, part);
  HECMW_log(HECMW_LOG_DEBUG, "Returned from kmetis(v5)\n");
#else
  int wgtflag    = 0;               /* flag of weight for edges */
  int numflag    = 0;               /* flag of stating number of index */
  int options[5] = {0, 0, 0, 0, 0}; /* options for kMETIS */

  HECMW_log(HECMW_LOG_DEBUG, "Entering kmetis(v4)...\n");
  METIS_PartGraphKway(&n, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                      &nparts, options, &edgecut, part);
  HECMW_log(HECMW_LOG_DEBUG, "Returned from kmetis(v4)\n");
#endif
#endif

  return edgecut;
}

static int pmetis_interface_with_weight(int n_vertex, int ncon, int n_domain,
                                        const int *xadj, const int *adjncy,
                                        const int *vwgt, int *part) {
  int edgecut = 0; /* number of edge-cut */
#ifdef HECMW_PART_WITH_METIS
  int n       = n_vertex; /* number of vertices */
  int *adjwgt = NULL;     /* weight for edges */
  int nparts  = n_domain; /* number of sub-domains */

#if defined(METIS_VER_MAJOR) && (METIS_VER_MAJOR == 5)
  int *vsize     = NULL;
  real_t *tpwgts = NULL;
  real_t *ubvec  = NULL;
  int *options   = NULL;

  HECMW_log(HECMW_LOG_DEBUG, "Entering pmetis(v5)...\n");
  METIS_PartGraphRecursive(&n, &ncon, (int *)xadj, (int *)adjncy, (int *)vwgt,
                           vsize, adjwgt, &nparts, tpwgts, ubvec, options,
                           &edgecut, part);
  HECMW_log(HECMW_LOG_DEBUG, "Returned from pmetis(v5)\n");
#else
  int wgtflag    = 0;               /* flag of weight for edges */
  int numflag    = 0;               /* flag of stating number of index */
  int options[5] = {0, 0, 0, 0, 0}; /* options for pMETIS */

  if (vwgt != NULL) wgtflag = 2;

  HECMW_log(HECMW_LOG_DEBUG, "Entering pmetis(v4)...\n");
  if (ncon == 1) {
    METIS_PartGraphRecursive(&n, (int *)xadj, (int *)adjncy, (int *)vwgt,
                             adjwgt, &wgtflag, &numflag, &nparts, options,
                             &edgecut, part);
  } else {
    METIS_mCPartGraphRecursive(&n, &ncon, (int *)xadj, (int *)adjncy,
                               (int *)vwgt, adjwgt, &wgtflag, &numflag, &nparts,
                               options, &edgecut, part);
  }
  HECMW_log(HECMW_LOG_DEBUG, "Returned from pmetis(v4)\n");
#endif
#endif

  return edgecut;
}

static int kmetis_interface_with_weight(int n_vertex, int ncon, int n_domain,
                                        const int *xadj, const int *adjncy,
                                        const int *vwgt, int *part) {
  int edgecut = 0; /* number of edge-cut */
#ifdef HECMW_PART_WITH_METIS
  int n       = n_vertex; /* number of vertices */
  int *adjwgt = NULL;     /* weight for edges */
  int nparts  = n_domain; /* number of sub-domains */

#if defined(METIS_VER_MAJOR) && (METIS_VER_MAJOR == 5)
  int *vsize     = NULL;
  real_t *tpwgts = NULL;
  real_t *ubvec  = NULL;
  int *options   = NULL;

  HECMW_log(HECMW_LOG_DEBUG, "Entering kmetis(v5)...\n");
  METIS_PartGraphKway(&n, &ncon, (int *)xadj, (int *)adjncy, (int *)vwgt, vsize,
                      adjwgt, &nparts, tpwgts, ubvec, options, &edgecut, part);
  HECMW_log(HECMW_LOG_DEBUG, "Returned from kmetis(v5)\n");
#else
  int wgtflag    = 0; /* flag of weight for edges */
  int numflag    = 0; /* flag of stating number of index */
  float *ubvec   = NULL;
  int options[5] = {0, 0, 0, 0, 0}; /* options for kMETIS */

  if (vwgt != NULL) wgtflag = 2;

  if (ncon > 1) {
    ubvec = (float *)HECMW_malloc(ncon * sizeof(float));
    if (ubvec == NULL) {
      HECMW_set_error(errno, "");
      return -1;
    }
  }

  HECMW_log(HECMW_LOG_DEBUG, "Entering kmetis(v4)...\n");
  if (ncon == 1) {
    METIS_PartGraphKway(&n, (int *)xadj, (int *)adjncy, (int *)vwgt, adjwgt,
                        &wgtflag, &numflag, &nparts, options, &edgecut, part);
  } else {
    METIS_mCPartGraphKway(&n, &ncon, (int *)xadj, (int *)adjncy, (int *)vwgt,
                          adjwgt, &wgtflag, &numflag, &nparts, ubvec, options,
                          &edgecut, part);
  }
  HECMW_log(HECMW_LOG_DEBUG, "Returned from kmetis(v4)\n");

  HECMW_free(ubvec);
#endif
#endif

  return edgecut;
}

static int contact_agg_mark_node_group(int *mark,
                                       struct hecmwST_local_mesh *global_mesh,
                                       int gid, int agg_id, int *agg_dup) {
  struct hecmwST_node_grp *ngrp = global_mesh->node_group;
  int istart, iend, i;

  HECMW_assert(0 < gid && gid <= ngrp->n_grp);

  istart = ngrp->grp_index[gid - 1];
  iend   = ngrp->grp_index[gid];
  for (i = istart; i < iend; i++) {
    int nid = ngrp->grp_item[i] - 1;
    HECMW_assert(0 <= nid && nid < global_mesh->n_node);
    if (0 <= mark[nid] && mark[nid] < agg_id) {
      /* the node is included in some other contact pair */
      if (*agg_dup == -1) {
        *agg_dup = mark[nid];
      } else if (mark[nid] != *agg_dup) {
        fprintf(stderr,
                "ERROR: node included in multiple node groups in different "
                "contact pairs,\n"
                "       which is not supported by CONTACT=AGGREGATE\n");
        HECMW_abort(HECMW_comm_get_comm());
      }
    }
    mark[nid] = agg_id;
  }
  return RTC_NORMAL;
}

static int HECMW_get_num_surf_node(int etype, int sid) {
  switch (etype) {
    case HECMW_ETYPE_TET1:
    case HECMW_ETYPE_PTT1:
      return 3;
    case HECMW_ETYPE_TET2:
    case HECMW_ETYPE_PTT2:
      return 6;
    case HECMW_ETYPE_HEX1:
    case HECMW_ETYPE_PTQ1:
      return 4;
    case HECMW_ETYPE_HEX2:
    case HECMW_ETYPE_PTQ2:
      return 8;
    case HECMW_ETYPE_PRI1:
      if (1 <= sid && sid <= 3) return 4;
      if (4 <= sid && sid <= 5) return 3;
    case HECMW_ETYPE_PRI2:
      if (1 <= sid && sid <= 3) return 8;
      if (4 <= sid && sid <= 5) return 6;
    default:
      fprintf(
          stderr,
          "ERROR: parallel contact analysis of elem type %d not supported\n",
          etype);
      return -1;
  }
  return -1;
}

static const int *HECMW_get_surf_node(int etype, int sid) {
  HECMW_assert(0 < sid);

  static const int elem_surf_tet1[4][3] = {
      {1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
  static const int elem_surf_tet2[4][6] = {{1, 4, 2, 9, 3, 8},
                                           {0, 7, 3, 9, 2, 5},
                                           {0, 6, 1, 8, 3, 7},
                                           {0, 5, 2, 4, 1, 6}};
  static const int elem_surf_hex1[6][4] = {{3, 0, 4, 7}, {1, 2, 6, 5},
                                           {0, 1, 5, 4}, {2, 3, 7, 6},
                                           {3, 2, 1, 0}, {4, 5, 6, 7}};
  static const int elem_surf_hex2[6][8] = {
      {3, 11, 0, 16, 4, 15, 7, 19}, {1, 9, 2, 18, 6, 13, 5, 17},
      {0, 8, 1, 17, 5, 12, 4, 16},  {2, 10, 3, 19, 7, 14, 6, 18},
      {3, 10, 2, 9, 1, 8, 0, 11},   {4, 12, 5, 13, 6, 14, 7, 15}};
  static const int elem_surf_pri1[5][4] = {
      {1, 2, 5, 4}, {2, 0, 3, 5}, {0, 1, 4, 3}, {2, 1, 0, -1}, {3, 4, 5, -1}};
  static const int elem_surf_pri2[5][8] = {{1, 6, 2, 14, 5, 9, 4, 13},
                                           {2, 7, 0, 12, 3, 10, 5, 14},
                                           {0, 8, 1, 13, 4, 11, 3, 12},
                                           {2, 6, 1, 8, 0, 7, -1, -1},
                                           {3, 11, 4, 9, 5, 10, -1, -1}};
  static const int elem_surf_ptt1[3] = {0, 1, 2};
  static const int elem_surf_ptt2[6] = {0, 1, 2, 3, 4, 5};
  static const int elem_surf_ptq1[4] = {0, 1, 2, 3};
  static const int elem_surf_ptq2[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  switch (etype) {
    case HECMW_ETYPE_TET1:
      return elem_surf_tet1[sid - 1];
    case HECMW_ETYPE_TET2:
      return elem_surf_tet2[sid - 1];
    case HECMW_ETYPE_HEX1:
      return elem_surf_hex1[sid - 1];
    case HECMW_ETYPE_HEX2:
      return elem_surf_hex2[sid - 1];
    case HECMW_ETYPE_PRI1:
      return elem_surf_pri1[sid - 1];
    case HECMW_ETYPE_PRI2:
      return elem_surf_pri2[sid - 1];
    case HECMW_ETYPE_PTT1:
      return elem_surf_ptt1;
    case HECMW_ETYPE_PTT2:
      return elem_surf_ptt2;
    case HECMW_ETYPE_PTQ1:
      return elem_surf_ptq1;
    case HECMW_ETYPE_PTQ2:
      return elem_surf_ptq2;
  }
  fprintf(stderr,
          "ERROR: parallel contact analysis of element type %d not supported\n",
          etype);
  return NULL;
}

static int HECMW_fistr_get_num_surf_node(int etype, int sid) {
  switch (etype) {
    case HECMW_ETYPE_TET1:
    case HECMW_ETYPE_PTT1:
      return 3;
    case HECMW_ETYPE_TET2:
    case HECMW_ETYPE_PTT2:
      return 6;
    case HECMW_ETYPE_HEX1:
    case HECMW_ETYPE_PTQ1:
      return 4;
    case HECMW_ETYPE_HEX2:
    case HECMW_ETYPE_PTQ2:
      return 8;
    case HECMW_ETYPE_PRI1:
      if (1 <= sid && sid <= 2) return 3;
      if (3 <= sid && sid <= 5) return 4;
    case HECMW_ETYPE_PRI2:
      if (1 <= sid && sid <= 2) return 6;
      if (3 <= sid && sid <= 5) return 8;
    default:
      fprintf(
          stderr,
          "ERROR: parallel contact analysis of elem type %d not supported\n",
          etype);
      return -1;
  }
  return -1;
}

static const int *HECMW_fistr_get_surf_node(int etype, int sid) {
  HECMW_assert(0 < sid);

  static const int elem_surf_tet1[4][3] = {
      {0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
  static const int elem_surf_tet2[4][6] = {{0, 6, 1, 4, 2, 5},
                                           {0, 6, 1, 8, 3, 7},
                                           {1, 4, 2, 9, 3, 8},
                                           {2, 5, 0, 9, 3, 7}};
  static const int elem_surf_hex1[6][4] = {{0, 1, 2, 3}, {4, 5, 6, 7},
                                           {0, 1, 5, 4}, {1, 2, 6, 5},
                                           {2, 3, 7, 6}, {3, 0, 4, 7}};
  static const int elem_surf_hex2[6][8] = {
      {0, 8, 1, 9, 2, 10, 3, 11},   {4, 12, 5, 13, 6, 14, 7, 15},
      {0, 8, 1, 17, 5, 12, 4, 16},  {1, 9, 2, 18, 6, 13, 5, 17},
      {2, 10, 3, 19, 7, 14, 6, 18}, {3, 11, 0, 16, 4, 15, 7, 19}};
  static const int elem_surf_pri1[5][4] = {
      {0, 1, 2, -1}, {3, 4, 5, -1}, {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}};
  static const int elem_surf_pri2[5][8] = {{0, 8, 1, 6, 2, 7, -1, -1},
                                           {3, 11, 4, 9, 5, 10, -1, -1},
                                           {0, 8, 1, 13, 4, 11, 3, 12},
                                           {1, 6, 2, 14, 5, 9, 4, 13},
                                           {2, 7, 0, 12, 3, 10, 5, 14}};
  static const int elem_surf_ptt1[3] = {0, 1, 2};
  static const int elem_surf_ptt2[6] = {0, 1, 2, 3, 4, 5};
  static const int elem_surf_ptq1[4] = {0, 1, 2, 3};
  static const int elem_surf_ptq2[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  switch (etype) {
    case HECMW_ETYPE_TET1:
      return elem_surf_tet1[sid - 1];
    case HECMW_ETYPE_TET2:
      return elem_surf_tet2[sid - 1];
    case HECMW_ETYPE_HEX1:
      return elem_surf_hex1[sid - 1];
    case HECMW_ETYPE_HEX2:
      return elem_surf_hex2[sid - 1];
    case HECMW_ETYPE_PRI1:
      return elem_surf_pri1[sid - 1];
    case HECMW_ETYPE_PRI2:
      return elem_surf_pri2[sid - 1];
    case HECMW_ETYPE_PTT1:
      return elem_surf_ptt1;
    case HECMW_ETYPE_PTT2:
      return elem_surf_ptt2;
    case HECMW_ETYPE_PTQ1:
      return elem_surf_ptq1;
    case HECMW_ETYPE_PTQ2:
      return elem_surf_ptq2;
  }
  fprintf(stderr,
          "ERROR: parallel contact analysis of element type %d not supported\n",
          etype);
  return NULL;
}

static int mark_contact_master_nodes(struct hecmwST_local_mesh *global_mesh,
                                     int *mark) {
  int i, j, k;
  struct hecmwST_contact_pair *cp = global_mesh->contact_pair;
  struct hecmwST_surf_grp *sgrp   = global_mesh->surf_group;

  for (i = 0; i < global_mesh->n_node; i++) {
    mark[i] = 0;
  }

  for (i = 0; i < cp->n_pair; i++) {
    int gid    = cp->master_grp_id[i];
    int jstart = sgrp->grp_index[gid - 1];
    int jend   = sgrp->grp_index[gid];
    for (j = jstart; j < jend; j++) {
      int eid = sgrp->grp_item[j * 2] - 1;
      int sid = sgrp->grp_item[j * 2 + 1];
      int *nop =
          global_mesh->elem_node_item + global_mesh->elem_node_index[eid];
      int etype = global_mesh->elem_type[eid];

      /** IF HEC-MW NUMBERING **/
      /* int num_snode = HECMW_get_num_surf_node(etype, sid); */
      /* const int *snode = HECMW_get_surf_node(etype, sid); */
      /** ELSE IF FrontISTR NUMBERING **/
      int num_snode    = HECMW_fistr_get_num_surf_node(etype, sid);
      const int *snode = HECMW_fistr_get_surf_node(etype, sid);
      /** END IF **/

      if (num_snode < 0 || snode == NULL) return RTC_ERROR;
      for (k = 0; k < num_snode; k++) {
        int nid = nop[snode[k]] - 1;
        HECMW_assert(0 <= nid && nid < global_mesh->n_node);
        mark[nid] = 1;
      }
    }
  }
  return RTC_NORMAL;
}

static int contact_agg_mark_surf_group(int *mark,
                                       struct hecmwST_local_mesh *global_mesh,
                                       int gid, int agg_id, int *agg_dup) {
  struct hecmwST_surf_grp *sgrp = global_mesh->surf_group;
  int istart, iend, i, j;

  HECMW_assert(0 < gid && gid <= sgrp->n_grp);

  /* get all nodes in the surface and mark them!!! */
  istart = sgrp->grp_index[gid - 1];
  iend   = sgrp->grp_index[gid];
  for (i = istart; i < iend; i++) {
    int eid   = sgrp->grp_item[i * 2] - 1;
    int sid   = sgrp->grp_item[i * 2 + 1];
    int *nop  = global_mesh->elem_node_item + global_mesh->elem_node_index[eid];
    int etype = global_mesh->elem_type[eid];
    /** IF HEC-WM NUMBERING **/
    /* int num_snode = HECMW_get_num_surf_node(etype, sid); */
    /* const int *snode = HECMW_get_surf_node(etype, sid); */
    /** ELSE IF FrontISTR NUMBERING **/
    int num_snode    = HECMW_fistr_get_num_surf_node(etype, sid);
    const int *snode = HECMW_fistr_get_surf_node(etype, sid);
    /** END IF **/
    if (num_snode < 0 || snode == NULL) return RTC_ERROR;
    for (j = 0; j < num_snode; j++) {
      int nid = nop[snode[j]] - 1;
      HECMW_assert(0 <= nid && nid < global_mesh->n_node);
      if (0 <= mark[nid] && mark[nid] < agg_id) {
        /* the node is included in some other contact pair */
        if (*agg_dup == -1) {
          *agg_dup = mark[nid];
        } else if (mark[nid] != *agg_dup) {
          fprintf(stderr,
                  "ERROR: node included in multiple surface groups in "
                  "different contact pairs,\n"
                  "       which is not supported by CONTACT=AGGREGATE\n");
          HECMW_abort(HECMW_comm_get_comm());
        }
      }
      mark[nid] = agg_id;
    }
  }
  return RTC_NORMAL;
}

static int metis_partition_nb_contact_agg(
    struct hecmwST_local_mesh *global_mesh,
    const struct hecmw_part_cont_data *cont_data,
    const struct hecmw_part_edge_data *edge_data) {
  int n_edgecut;
  int *node_graph_index = NULL; /* index for nodal graph */
  int *node_graph_item  = NULL; /* member of nodal graph */
  int *belong_domain    = NULL;
  int rtc;
  int i;
  struct hecmwST_contact_pair *cp;
  int *mark;
  int agg_id, agg_dup, gid;
  int n_node2;
  const int *node_graph_index2;
  const int *node_graph_item2;
  int *node_weight2;
  struct hecmw_graph graph1, graph2;
  const int ncon = 1;

  HECMW_assert(global_mesh->hecmw_flag_partcontact ==
               HECMW_FLAG_PARTCONTACT_AGGREGATE);

  node_graph_index = (int *)HECMW_calloc(global_mesh->n_node + 1, sizeof(int));
  if (node_graph_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  node_graph_item = (int *)HECMW_malloc(sizeof(int) * edge_data->n_edge * 2);
  if (node_graph_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  HECMW_log(HECMW_LOG_DEBUG, "Starting creation of node graph...\n");

  rtc = create_node_graph(global_mesh, edge_data, node_graph_index,
                          node_graph_item);
  if (rtc != RTC_NORMAL) goto error;

  HECMW_log(HECMW_LOG_DEBUG, "Creation of node graph done\n");

  HECMW_log(HECMW_LOG_DEBUG, "Partitioning mode: contact-aggregate\n");

  HECMW_log(HECMW_LOG_DEBUG, "Starting aggregation of contact pairs...\n");

  /* aggregate contact pair if requested */
  cp   = global_mesh->contact_pair;
  mark = (int *)HECMW_malloc(global_mesh->n_node * sizeof(int));
  if (mark == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  for (i = 0; i < global_mesh->n_node; i++) {
    mark[i] = -1;
  }
  agg_id = 0;
  /* mark contact pairs */
  for (i = 0; i < cp->n_pair; i++) {
    agg_dup = -1;
    /* slave */
    if (cp->type[i] == HECMW_CONTACT_TYPE_NODE_SURF) {
      gid = cp->slave_grp_id[i];
      rtc =
          contact_agg_mark_node_group(mark, global_mesh, gid, agg_id, &agg_dup);
      if (rtc != RTC_NORMAL) goto error;
    } else if(cp->type[i] == HECMW_CONTACT_TYPE_SURF_SURF) {
      gid = cp->slave_grp_id[i];
      rtc =
          contact_agg_mark_surf_group(mark, global_mesh, gid, agg_id, &agg_dup);
      if (rtc != RTC_NORMAL) goto error;
    } else if(cp->type[i] == HECMW_CONTACT_TYPE_NODE_ELEM) {
      gid = cp->slave_grp_id[i];
      rtc =
          contact_agg_mark_surf_group(mark, global_mesh, gid, agg_id, &agg_dup);
      if (rtc != RTC_NORMAL) goto error;
    }
    /* master */
    gid = cp->master_grp_id[i];
    rtc = contact_agg_mark_surf_group(mark, global_mesh, gid, agg_id, &agg_dup);
    if (rtc != RTC_NORMAL) goto error;

    if (agg_dup >= 0) {
      for (i = 0; i < global_mesh->n_node; i++) {
        if (mark[i] == agg_id) {
          mark[i] = agg_dup;
        }
      }
    } else {
      agg_id++;
    }
  }
  /* mark other nodes */
  for (i = 0; i < global_mesh->n_node; i++) {
    if (mark[i] < 0) {
      mark[i] = agg_id++;
    }
  }
  n_node2 = agg_id;

  /* degenerate node graph */
  rtc = HECMW_graph_init_with_arrays(&graph1, global_mesh->n_node,
                                     node_graph_index, node_graph_item);
  if (rtc != RTC_NORMAL) goto error;
  rtc = HECMW_graph_init(&graph2);
  if (rtc != RTC_NORMAL) goto error;
  rtc = HECMW_graph_degeneGraph(&graph2, &graph1, n_node2, mark);
  if (rtc != RTC_NORMAL) goto error;
  HECMW_graph_finalize(&graph1);
  node_graph_index2 = HECMW_graph_getEdgeIndex(&graph2);
  node_graph_item2  = HECMW_graph_getEdgeItem(&graph2);

  node_weight2 = (int *)HECMW_calloc(n_node2, sizeof(int));
  if (node_weight2 == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  for (i = 0; i < global_mesh->n_node; i++) {
    node_weight2[mark[i]] += 1;
  }

  HECMW_log(HECMW_LOG_DEBUG, "Aggregation of contact pairs done\n");

  belong_domain = (int *)HECMW_calloc(n_node2, sizeof(int));
  if (belong_domain == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  switch (cont_data->method) {
    case HECMW_PART_METHOD_PMETIS: /* pMETIS */
      n_edgecut = pmetis_interface_with_weight(
          n_node2, ncon, global_mesh->n_subdomain, node_graph_index2,
          node_graph_item2, node_weight2, belong_domain);
      if (n_edgecut < 0) goto error;
      break;

    case HECMW_PART_METHOD_KMETIS: /* kMETIS */
      n_edgecut = kmetis_interface_with_weight(
          n_node2, ncon, global_mesh->n_subdomain, node_graph_index2,
          node_graph_item2, node_weight2, belong_domain);
      if (n_edgecut < 0) goto error;
      break;

    default:
      HECMW_set_error(HECMW_PART_E_INVALID_PMETHOD, "");
      goto error;
  }

  for (i = 0; i < global_mesh->n_node; i++) {
    global_mesh->node_ID[2 * i + 1] = belong_domain[mark[i]];
  }

  HECMW_graph_finalize(&graph2);
  HECMW_free(node_graph_index);
  HECMW_free(node_graph_item);
  HECMW_free(mark);
  HECMW_free(node_weight2);
  HECMW_free(belong_domain);

  return n_edgecut;

error:
  HECMW_free(node_graph_index);
  HECMW_free(node_graph_item);
  HECMW_free(mark);
  HECMW_free(node_weight2);
  HECMW_free(belong_domain);

  return -1;
}

static int metis_partition_nb_contact_dist(
    struct hecmwST_local_mesh *global_mesh,
    const struct hecmw_part_cont_data *cont_data,
    const struct hecmw_part_edge_data *edge_data) {
  int n_edgecut;
  int *node_graph_index = NULL; /* index for nodal graph */
  int *node_graph_item  = NULL; /* member of nodal graph */
  int *belong_domain    = NULL;
  int rtc;
  int i;
  int ncon;
  int *node_weight = NULL;
  int *mark        = NULL;

  HECMW_assert(
      global_mesh->hecmw_flag_partcontact == HECMW_FLAG_PARTCONTACT_SIMPLE ||
      global_mesh->hecmw_flag_partcontact == HECMW_FLAG_PARTCONTACT_DISTRIBUTE);

  node_graph_index = (int *)HECMW_calloc(global_mesh->n_node + 1, sizeof(int));
  if (node_graph_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  node_graph_item = (int *)HECMW_malloc(sizeof(int) * edge_data->n_edge * 2);
  if (node_graph_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  HECMW_log(HECMW_LOG_DEBUG, "Starting creation of node graph...\n");

  rtc = create_node_graph(global_mesh, edge_data, node_graph_index,
                          node_graph_item);
  if (rtc != RTC_NORMAL) goto error;

  HECMW_log(HECMW_LOG_DEBUG, "Creation of node graph done\n");

  if (global_mesh->hecmw_flag_partcontact == HECMW_FLAG_PARTCONTACT_SIMPLE) {
    HECMW_log(HECMW_LOG_DEBUG, "Partitioning mode: contact-simple\n");

    ncon        = 1;
    node_weight = NULL;
  } else /* HECMW_FLAG_PARTCONTACT_DISTRIBUTE */
  {
    HECMW_log(HECMW_LOG_DEBUG, "Partitioning mode: contact-distribute\n");

    ncon = 2;

    mark = (int *)HECMW_calloc(global_mesh->n_node, sizeof(int));
    if (mark == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    }

    rtc = mark_contact_master_nodes(global_mesh, mark);
    if (rtc != RTC_NORMAL) goto error;

    node_weight = (int *)HECMW_calloc(global_mesh->n_node * ncon, sizeof(int));
    if (node_weight == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    }

    for (i = 0; i < global_mesh->n_node; i++) {
      /* 1st condition: distribute nodes equally */
      node_weight[i * ncon] = 1;
      /* 2nd condition: distribute master nodes equally */
      node_weight[i * ncon + 1] = mark[i];
    }

    HECMW_free(mark);
  }

  belong_domain = (int *)HECMW_calloc(global_mesh->n_node, sizeof(int));
  if (belong_domain == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  switch (cont_data->method) {
    case HECMW_PART_METHOD_PMETIS: /* pMETIS */
      n_edgecut = pmetis_interface_with_weight(
          global_mesh->n_node, ncon, global_mesh->n_subdomain, node_graph_index,
          node_graph_item, node_weight, belong_domain);
      if (n_edgecut < 0) goto error;
      break;

    case HECMW_PART_METHOD_KMETIS: /* kMETIS */
      n_edgecut = kmetis_interface_with_weight(
          global_mesh->n_node, ncon, global_mesh->n_subdomain, node_graph_index,
          node_graph_item, node_weight, belong_domain);
      if (n_edgecut < 0) goto error;
      break;

    default:
      HECMW_set_error(HECMW_PART_E_INVALID_PMETHOD, "");
      goto error;
  }

  for (i = 0; i < global_mesh->n_node; i++) {
    global_mesh->node_ID[2 * i + 1] = belong_domain[i];
  }

  HECMW_free(node_graph_index);
  HECMW_free(node_graph_item);
  HECMW_free(belong_domain);
  if (node_weight) HECMW_free(node_weight);

  return n_edgecut;

error:
  HECMW_free(node_graph_index);
  HECMW_free(node_graph_item);
  HECMW_free(belong_domain);
  if (node_weight) HECMW_free(node_weight);
  if (mark) HECMW_free(mark);

  return -1;
}

static int metis_partition_nb_default(
    struct hecmwST_local_mesh *global_mesh,
    const struct hecmw_part_cont_data *cont_data,
    const struct hecmw_part_edge_data *edge_data) {
  int n_edgecut;
  int *node_graph_index = NULL; /* index for nodal graph */
  int *node_graph_item  = NULL; /* member of nodal graph */
  int *belong_domain    = NULL;
  int rtc;
  int i;

  node_graph_index = (int *)HECMW_calloc(global_mesh->n_node + 1, sizeof(int));
  if (node_graph_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  node_graph_item = (int *)HECMW_malloc(sizeof(int) * edge_data->n_edge * 2);
  if (node_graph_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  HECMW_log(HECMW_LOG_DEBUG, "Starting creation of node graph...\n");

  rtc = create_node_graph(global_mesh, edge_data, node_graph_index,
                          node_graph_item);
  if (rtc != RTC_NORMAL) goto error;

  HECMW_log(HECMW_LOG_DEBUG, "Creation of node graph done\n");

  belong_domain = (int *)HECMW_calloc(global_mesh->n_node, sizeof(int));
  if (belong_domain == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  HECMW_log(HECMW_LOG_DEBUG, "Partitioning mode: default\n");

  switch (cont_data->method) {
    case HECMW_PART_METHOD_PMETIS: /* pMETIS */
      n_edgecut =
          pmetis_interface(global_mesh->n_node, global_mesh->n_subdomain,
                           node_graph_index, node_graph_item, belong_domain);
      if (n_edgecut < 0) goto error;
      break;

    case HECMW_PART_METHOD_KMETIS: /* kMETIS */
      n_edgecut =
          kmetis_interface(global_mesh->n_node, global_mesh->n_subdomain,
                           node_graph_index, node_graph_item, belong_domain);
      if (n_edgecut < 0) goto error;
      break;

    default:
      HECMW_set_error(HECMW_PART_E_INVALID_PMETHOD, "");
      goto error;
  }

  for (i = 0; i < global_mesh->n_node; i++) {
    global_mesh->node_ID[2 * i + 1] = belong_domain[i];
  }

  HECMW_free(node_graph_index);
  HECMW_free(node_graph_item);
  HECMW_free(belong_domain);

  return n_edgecut;

error:
  HECMW_free(node_graph_index);
  HECMW_free(node_graph_item);
  HECMW_free(belong_domain);

  return -1;
}

static int metis_partition_nb(struct hecmwST_local_mesh *global_mesh,
                              const struct hecmw_part_cont_data *cont_data,
                              const struct hecmw_part_edge_data *edge_data) {
  if (global_mesh->contact_pair->n_pair > 0) {
    switch (global_mesh->hecmw_flag_partcontact) {
      case HECMW_FLAG_PARTCONTACT_AGGREGATE:
        return metis_partition_nb_contact_agg(global_mesh, cont_data,
                                              edge_data);

      case HECMW_FLAG_PARTCONTACT_DISTRIBUTE:
      case HECMW_FLAG_PARTCONTACT_SIMPLE:
        return metis_partition_nb_contact_dist(global_mesh, cont_data,
                                               edge_data);

      default:
        return -1;
    }
  } else {
    return metis_partition_nb_default(global_mesh, cont_data, edge_data);
  }
}

static int metis_partition_eb(struct hecmwST_local_mesh *global_mesh,
                              const struct hecmw_part_cont_data *cont_data,
                              int *elem_graph_index, int *elem_graph_item) {
  int n_edgecut;
  int *belong_domain = NULL;
  int i;

  belong_domain = (int *)HECMW_calloc(global_mesh->n_elem, sizeof(int));
  if (belong_domain == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  switch (cont_data->method) {
    case HECMW_PART_METHOD_PMETIS: /* pMETIS */
      n_edgecut =
          pmetis_interface(global_mesh->n_elem, global_mesh->n_subdomain,
                           elem_graph_index, elem_graph_item, belong_domain);
      if (n_edgecut < 0) goto error;
      break;

    case HECMW_PART_METHOD_KMETIS: /* kMETIS */
      n_edgecut =
          kmetis_interface(global_mesh->n_elem, global_mesh->n_subdomain,
                           elem_graph_index, elem_graph_item, belong_domain);
      if (n_edgecut < 0) goto error;
      break;

    default:
      HECMW_set_error(HECMW_PART_E_INVALID_PMETHOD, "");
      goto error;
  }

  for (i = 0; i < global_mesh->n_elem; i++) {
    global_mesh->elem_ID[2 * i + 1] = belong_domain[i];
  }

  HECMW_free(belong_domain);

  return n_edgecut;

error:
  HECMW_free(belong_domain);

  return -1;
}

/*------------------------------------------------------------------------------------------------*/

#define LINEBUF_SIZE 1023

static int read_part_file(
    const char *part_file_name,
    int n,
    int n_domain,
    int *wnum) {
  FILE *fpart;
  char linebuf[LINEBUF_SIZE + 1];
  int rtc, n_in, n_domain_in;
  int i, part;
  int *count_dom;

  fpart = fopen(part_file_name, "r");
  if (fpart == NULL) {
    HECMW_set_error(HECMW_PART_E_NO_SUCH_FILE, "%s", part_file_name);
    goto error;
  }

  /* read n and n_domain */
  if (fgets(linebuf, LINEBUF_SIZE, fpart) == NULL) {
    HECMW_set_error(HECMW_PART_E_PART_EOF, "read_part_file");
    goto error;
  }
  rtc = sscanf(linebuf, "%d %d", &n_in, &n_domain_in);
  if (rtc != 2) {
    HECMW_set_error(HECMW_PART_E_PART_INVALID_FORMAT, "");
    goto error;
  }

  if (n_in != n) {
    HECMW_set_error(HECMW_PART_E_PART_N, "");
    goto error;
  }
  if (n_domain_in != n_domain) {
    HECMW_set_error(HECMW_PART_E_PART_NDOMAIN, "");
    goto error;
  }

  count_dom = (int *) HECMW_calloc(n_domain, sizeof(int));
  if (count_dom == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  /* read part array and count members in each domain */
  for (i = 0; i < n; i++) {
    if (fgets(linebuf, LINEBUF_SIZE, fpart) == NULL) {
      HECMW_set_error(HECMW_PART_E_PART_EOF, "");
      goto error;
    }
    rtc = sscanf(linebuf, "%d", &part);
    if (rtc != 1) {
      HECMW_set_error(HECMW_PART_E_PART_INVALID_FORMAT, "");
      goto error;
    }

    if (part < 0 || n_domain <= part) {
      HECMW_set_error(HECMW_PART_E_PART_INVALID_PART, "%d", part);
      goto error;
    }

    count_dom[part]++;

    wnum[2*i+1] = part;
  }

  /* check for empty domain */
  for (i = 0; i < n_domain; i++) {
    if (count_dom[i] == 0) {
      HECMW_set_error(HECMW_PART_E_PART_EMPTY_DOMAIN, "%d", i);
      goto error;
    }
  }

  fclose(fpart);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int write_part_file(
    const char *part_file_name,
    int n,
    int n_domain,
    const int *wnum) {
  FILE *fpart;
  int i;

  fpart = fopen(part_file_name, "w");
  if (fpart == NULL) {
    HECMW_set_error(HECMW_PART_E_NO_SUCH_FILE, "%s", part_file_name);
    goto error;
  }

  fprintf(fpart, "%d %d\n", n, n_domain);

  for (i = 0; i < n; i++) {
    fprintf(fpart, "%d\n", wnum[2*i+1]);
  }

  fclose(fpart);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*------------------------------------------------------------------------------------------------*/

static int user_partition(
    int n,
    int n_domain,
    int *wnum,
    const char *part_file_name) {
  int rtc;

  rtc = read_part_file(part_file_name, n, n_domain, wnum);
  if (rtc != RTC_NORMAL) goto error;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int user_partition_nb(
    struct hecmwST_local_mesh *global_mesh,
    const struct hecmw_part_cont_data *cont_data) {
  return user_partition(global_mesh->n_node, global_mesh->n_subdomain,
                        global_mesh->node_ID, cont_data->part_file_name);
}

static int user_partition_eb(
    struct hecmwST_local_mesh *global_mesh,
    const struct hecmw_part_cont_data *cont_data) {
  return user_partition(global_mesh->n_elem, global_mesh->n_subdomain,
                        global_mesh->elem_ID, cont_data->part_file_name);
}

/*------------------------------------------------------------------------------------------------*/

static int print_part(
    struct hecmwST_local_mesh *global_mesh,
    const char *part_file_name) {
  int rtc;

  switch (global_mesh->hecmw_flag_parttype) {
  case HECMW_FLAG_PARTTYPE_NODEBASED:
    rtc = write_part_file(part_file_name, global_mesh->n_node,
			  global_mesh->n_subdomain, global_mesh->node_ID);
    if (rtc != RTC_NORMAL) goto error;

    break;

  case HECMW_FLAG_PARTTYPE_ELEMBASED:
    rtc = write_part_file(part_file_name, global_mesh->n_elem,
			  global_mesh->n_subdomain, global_mesh->elem_ID);
    if (rtc != RTC_NORMAL) goto error;

    break;

  default:
    HECMW_set_error(HECMW_PART_E_INVALID_PTYPE, "");
    goto error;
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*------------------------------------------------------------------------------------------------*/

static int count_edgecut(
    const struct hecmw_part_edge_data *edge_data,
    const int *wnum) {
  int i;
  int n_edgecut = 0;

  for (i = 0; i < edge_data->n_edge; i++) {
    if (wnum[2 * (edge_data->edge_node_item[2 * i] - 1) + 1] !=
        wnum[2 * (edge_data->edge_node_item[2 * i + 1] - 1) + 1]) {
      n_edgecut++;
    }
  }

  return n_edgecut;
}

/*------------------------------------------------------------------------------------------------*/

static int set_node_belong_domain_nb(
    struct hecmwST_local_mesh *global_mesh,
    const struct hecmw_part_cont_data *cont_data) {
  struct hecmw_part_edge_data *edge_data = NULL;
  int n_edgecut;
  int rtc;
  long long int i;

  edge_data = (struct hecmw_part_edge_data *)HECMW_malloc(
      sizeof(struct hecmw_part_edge_data));
  if (edge_data == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  } else {
    edge_data->n_edge         = 0;
    edge_data->edge_node_item = NULL;
  }

  HECMW_log(HECMW_LOG_DEBUG, "Starting creation of mesh edge info...\n");

  rtc = HECMW_mesh_edge_info(global_mesh, edge_data);
  if (rtc != 0) goto error;

  HECMW_log(HECMW_LOG_DEBUG, "Creation of mesh edge info done\n");

  switch (cont_data->method) {
    case HECMW_PART_METHOD_RCB: /* RCB */
      rtc = rcb_partition(global_mesh->n_node, global_mesh->node,
                          global_mesh->node_ID, cont_data);
      if (rtc != RTC_NORMAL) goto error;

      n_edgecut = count_edgecut(edge_data, global_mesh->node_ID);

      break;

    case HECMW_PART_METHOD_KMETIS: /* kMETIS */
    case HECMW_PART_METHOD_PMETIS: /* pMETIS */
      n_edgecut = metis_partition_nb(global_mesh, cont_data, edge_data);
      if (n_edgecut < 0) goto error;

      break;

    case HECMW_PART_METHOD_USER: /* USER */
      rtc = user_partition_nb(global_mesh, cont_data);
      if (rtc != RTC_NORMAL) goto error;

      n_edgecut = count_edgecut(edge_data, global_mesh->node_ID);

      break;

    default:
      HECMW_set_error(HECMW_PART_E_INVALID_PMETHOD, "");
      goto error;
  }

  rtc = HECMW_part_set_log_n_edgecut(edge_data->n_edge, n_edgecut);
  if (rtc != RTC_NORMAL) goto error;

  /* commented out by K.Goto; begin */
  /* rtc = eqn_block( global_mesh ); */
  /* if( rtc != RTC_NORMAL )  goto error; */
  /* commented out by K.Goto; end */

  HECMW_free(edge_data->edge_node_item);
  HECMW_free(edge_data);

  return RTC_NORMAL;

error:
  if (edge_data) {
    HECMW_free(edge_data->edge_node_item);
  }
  HECMW_free(edge_data);

  return RTC_ERROR;
}

static int set_node_belong_domain_eb(struct hecmwST_local_mesh *global_mesh) {
  int node;
  int i, j;

  for (i = 0; i < global_mesh->n_node; i++) {
    global_mesh->node_ID[2 * i + 1] = global_mesh->n_subdomain;
  }

  for (i = 0; i < global_mesh->n_elem; i++) {
    for (j = global_mesh->elem_node_index[i];
         j < global_mesh->elem_node_index[i + 1]; j++) {
      node = global_mesh->elem_node_item[j];
      if (global_mesh->elem_ID[2 * i + 1] <
          global_mesh->node_ID[2 * (node - 1) + 1]) {
        global_mesh->node_ID[2 * (node - 1) + 1] =
            global_mesh->elem_ID[2 * i + 1];
      }
    }
  }

  return RTC_NORMAL;
}

static int set_local_node_id(struct hecmwST_local_mesh *global_mesh) {
  int *counter;
  int j, domain;

  counter = (int *)HECMW_calloc(global_mesh->n_subdomain, sizeof(int));
  if (counter == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (j = 0; j < global_mesh->n_node; j++) {
    domain                      = global_mesh->node_ID[2 * j + 1];
    global_mesh->node_ID[2 * j] = ++counter[domain];
  }

  HECMW_free(counter);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int wnumbering_node(struct hecmwST_local_mesh *global_mesh,
                           const struct hecmw_part_cont_data *cont_data) {
  int rtc;
  int i;

  HECMW_free(global_mesh->node_ID);
  global_mesh->node_ID =
      (int *)HECMW_malloc(sizeof(int) * global_mesh->n_node * 2);
  if (global_mesh->node_ID == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  } else {
    for (i = 0; i < global_mesh->n_node; i++) {
      global_mesh->node_ID[2 * i]     = i + 1;
      global_mesh->node_ID[2 * i + 1] = 0;
    }
  }

  if (global_mesh->n_subdomain == 1) return RTC_NORMAL;

  switch (global_mesh->hecmw_flag_parttype) {
    case HECMW_FLAG_PARTTYPE_NODEBASED: /* for node-based partitioning */
      rtc = set_node_belong_domain_nb(global_mesh, cont_data);
      if (rtc != RTC_NORMAL) goto error;
      break;

    case HECMW_FLAG_PARTTYPE_ELEMBASED: /* for element-based partitioning */
      rtc = set_node_belong_domain_eb(global_mesh);
      if (rtc != RTC_NORMAL) goto error;
      break;

    default:
      HECMW_set_error(HECMW_PART_E_INVALID_PTYPE, "");
      goto error;
  }

  rtc = set_local_node_id(global_mesh);
  if (rtc != RTC_NORMAL) goto error;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*------------------------------------------------------------------------------------------------*/

static int set_elem_belong_domain_nb(struct hecmwST_local_mesh *global_mesh) {
  int node, node_domain, min_domain;
  int i, j;

  for (i = 0; i < global_mesh->n_elem; i++) {
    min_domain = global_mesh->n_subdomain;
    for (j = global_mesh->elem_node_index[i];
         j < global_mesh->elem_node_index[i + 1]; j++) {
      node        = global_mesh->elem_node_item[j];
      node_domain = global_mesh->node_ID[2 * (node - 1) + 1];
      if (node_domain < min_domain) {
        min_domain = node_domain;
      }
    }
    global_mesh->elem_ID[2 * i + 1] = min_domain;
  }

  return RTC_NORMAL;
}

static int count_edge_for_eb(const struct hecmwST_local_mesh *global_mesh,
                             struct hecmw_part_edge_data *elem_data,
                             int *elem_graph_index, int *elem_graph_item) {
  int rtc;
  long long int eid;
  int i, j;

  rtc = HECMW_mesh_hsort_edge_init(global_mesh->n_node, global_mesh->n_elem);
  if (rtc != RTC_NORMAL) goto error;

  for (i = 0; i < global_mesh->n_elem; i++) {
    for (j = elem_graph_index[i]; j < elem_graph_index[i + 1]; j++) {
      eid = HECMW_mesh_hsort_edge(i + 1, elem_graph_item[j] + 1);
      if (eid < 0) goto error;
    }
  }

  elem_data->n_edge = HECMW_mesh_hsort_edge_get_n();
  if (elem_data->n_edge < 0) goto error;

  elem_data->edge_node_item = HECMW_mesh_hsort_edge_get_v();
  if (elem_data->edge_node_item == NULL) goto error;

  HECMW_mesh_hsort_edge_final();

  return RTC_NORMAL;

error:
  HECMW_mesh_hsort_edge_final();

  return RTC_ERROR;
}

static int set_elem_belong_domain_eb(
    struct hecmwST_local_mesh *global_mesh,
    const struct hecmw_part_cont_data *cont_data) {
  int n_edgecut                          = 0;
  int *elem_graph_index                  = NULL;
  int *elem_graph_item                   = NULL;
  struct hecmw_part_edge_data *elem_data = NULL;
  int rtc;
  long long int i;

  elem_graph_index = (int *)HECMW_calloc(global_mesh->n_elem + 1, sizeof(int));
  if (elem_graph_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  elem_data = (struct hecmw_part_edge_data *)HECMW_malloc(
      sizeof(struct hecmw_part_edge_data));
  if (elem_data == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  } else {
    elem_data->n_edge         = 0;
    elem_data->edge_node_item = NULL;
  }

  HECMW_log(HECMW_LOG_DEBUG, "Starting creation of elem graph...\n");

  elem_graph_item = create_elem_graph(global_mesh, elem_graph_index);
  if (elem_graph_item == NULL) goto error;

  HECMW_log(HECMW_LOG_DEBUG, "Creation of elem graph done\n");

  rtc = count_edge_for_eb(global_mesh, elem_data, elem_graph_index,
                          elem_graph_item);
  if (rtc != RTC_NORMAL) goto error;

  switch (cont_data->method) {
    case HECMW_PART_METHOD_RCB: /* RCB */
      rtc = rcb_partition_eb(global_mesh, cont_data);
      if (rtc != RTC_NORMAL) goto error;

      n_edgecut = count_edgecut(elem_data, global_mesh->elem_ID);

      break;

    case HECMW_PART_METHOD_PMETIS: /* pMETIS */
    case HECMW_PART_METHOD_KMETIS: /* kMETIS */
      n_edgecut = metis_partition_eb(global_mesh, cont_data, elem_graph_index,
                                     elem_graph_item);
      if (n_edgecut < 0) goto error;

      break;

    case HECMW_PART_METHOD_USER: /* USER */
      rtc = user_partition_eb(global_mesh, cont_data);
      if (rtc != RTC_NORMAL) goto error;

      n_edgecut = count_edgecut(elem_data, global_mesh->elem_ID);

      break;

    default:
      HECMW_set_error(HECMW_PART_E_INVALID_PMETHOD, "");
      goto error;
  }

  rtc = HECMW_part_set_log_n_edgecut(elem_data->n_edge, n_edgecut);
  if (rtc != RTC_NORMAL) goto error;

  HECMW_free(elem_graph_index);
  HECMW_free(elem_graph_item);
  HECMW_free(elem_data->edge_node_item);
  HECMW_free(elem_data);

  return RTC_NORMAL;

error:
  HECMW_free(elem_graph_index);
  HECMW_free(elem_graph_item);
  if (elem_data) {
    HECMW_free(elem_data->edge_node_item);
  }
  HECMW_free(elem_data);

  return RTC_ERROR;
}

static int set_local_elem_id(struct hecmwST_local_mesh *global_mesh) {
  int *counter;
  int j, domain;

  counter = (int *)HECMW_calloc(global_mesh->n_subdomain, sizeof(int));
  if (counter == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (j = 0; j < global_mesh->n_elem; j++) {
    domain                      = global_mesh->elem_ID[2 * j + 1];
    global_mesh->elem_ID[2 * j] = ++counter[domain];
  }

  HECMW_free(counter);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int wnumbering_elem(struct hecmwST_local_mesh *global_mesh,
                           const struct hecmw_part_cont_data *cont_data) {
  int rtc;
  int i;

  HECMW_free(global_mesh->elem_ID);
  global_mesh->elem_ID =
      (int *)HECMW_malloc(sizeof(int) * global_mesh->n_elem * 2);
  if (global_mesh->elem_ID == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  } else {
    for (i = 0; i < global_mesh->n_elem; i++) {
      global_mesh->elem_ID[2 * i]     = i + 1;
      global_mesh->elem_ID[2 * i + 1] = 0;
    }
  }

  if (global_mesh->n_subdomain == 1) return RTC_NORMAL;

  switch (global_mesh->hecmw_flag_parttype) {
    case HECMW_FLAG_PARTTYPE_NODEBASED: /* for node-based partitioning */
      rtc = set_elem_belong_domain_nb(global_mesh);
      if (rtc != RTC_NORMAL) goto error;

      break;

    case HECMW_FLAG_PARTTYPE_ELEMBASED: /* for element-based partitioning */
      rtc = set_elem_belong_domain_eb(global_mesh, cont_data);
      if (rtc != RTC_NORMAL) goto error;

      break;

    default:
      HECMW_set_error(HECMW_PART_E_INVALID_PTYPE, "");
      goto error;
  }

  rtc = set_local_elem_id(global_mesh);
  if (rtc != RTC_NORMAL) goto error;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int wnumbering(struct hecmwST_local_mesh *global_mesh,
                      const struct hecmw_part_cont_data *cont_data) {
  int rtc;

  HECMW_assert(global_mesh);
  HECMW_assert(cont_data);

  HECMW_log(HECMW_LOG_DEBUG, "Starting double numbering...");

  switch (global_mesh->hecmw_flag_parttype) {
    case HECMW_FLAG_PARTTYPE_NODEBASED: /* for node-based partitioning */
      rtc = wnumbering_node(global_mesh, cont_data);
      if (rtc != RTC_NORMAL) goto error;

      rtc = wnumbering_elem(global_mesh, cont_data);
      if (rtc != RTC_NORMAL) goto error;

      break;

    case HECMW_FLAG_PARTTYPE_ELEMBASED: /* for element-based partitioning */

      rtc = wnumbering_elem(global_mesh, cont_data);
      if (rtc != RTC_NORMAL) goto error;

      rtc = wnumbering_node(global_mesh, cont_data);
      if (rtc != RTC_NORMAL) goto error;

      break;

    default:
      HECMW_set_error(HECMW_PART_E_INVALID_PTYPE, "");
      goto error;
  }

  HECMW_log(HECMW_LOG_DEBUG, "Double numbering done");

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*==================================================================================================


  create neighboring domain & communication information


==================================================================================================*/

/*K. Inagaki */
static int mask_node_by_domain(const struct hecmwST_local_mesh *global_mesh,
                               char *node_flag, int current_domain) {
  int i, node;

  for (i = 0; i < n_int_nlist[current_domain]; i++) {
    node = int_nlist[current_domain][i];
    MASK_BIT(node_flag[node - 1], INTERNAL);
  }

  return RTC_NORMAL;
}

static int mask_elem_by_domain(const struct hecmwST_local_mesh *global_mesh,
                               char *elem_flag, int current_domain) {
  int i;

  for (i = 0; i < global_mesh->n_elem; i++) {
    (global_mesh->elem_ID[2 * i + 1] == current_domain)
        ? MASK_BIT(elem_flag[i], INTERNAL)
        : MASK_BIT(elem_flag[i], EXTERNAL);
  }

  return RTC_NORMAL;
}

/*K. Inagaki */
static int mask_elem_by_domain_mod(char *elem_flag, int current_domain) {
  int i, elem;

  for (i = 0; i < n_int_elist[current_domain]; i++) {
    elem = int_elist[current_domain][i];
    MASK_BIT(elem_flag[elem - 1], INTERNAL);
  }

  return RTC_NORMAL;
}

#if 0
/* For Additional overlap for explicit DOF elimination for MPC */
/* NO LONGER NEEDED because node-migration implemented */
static int mask_slave_node(const struct hecmwST_local_mesh *global_mesh,
                           char *node_flag, int current_domain) {
  int i;

  for (i = 0; i < global_mesh->mpc->n_mpc; i++) {
    int j0, je, slave, master, j, evalsum;
    j0    = global_mesh->mpc->mpc_index[i];
    je    = global_mesh->mpc->mpc_index[i + 1];
    slave = global_mesh->mpc->mpc_item[j0];

    /* mask all slave nodes */
    MASK_BIT(node_flag[slave - 1], MASK);

    /* mark slave nodes that have mpc-link across the boundary */
    evalsum = 0;
    for (j = j0 + 1; j < je; j++) {
      master = global_mesh->mpc->mpc_item[j];
      if (EVAL_BIT(node_flag[slave - 1], INTERNAL) ^ /* exclusive or */
          EVAL_BIT(node_flag[master - 1], INTERNAL)) {
        evalsum++;
      }
    }
    if (evalsum) {
      MASK_BIT(node_flag[slave - 1], MARK);
    }
  }
  return RTC_NORMAL;
}
#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * - - - - - - - - - */

/*K. Inagaki */
static int mask_overlap_elem(char *elem_flag, int domain) {
  int i, elem;

  for (i = 0; i < n_bnd_elist[2 * domain + 1]; i++) {
    elem = bnd_elist[domain][i];
    MASK_BIT(elem_flag[elem - 1], OVERLAP);
    MASK_BIT(elem_flag[elem - 1], BOUNDARY);
  }

  return RTC_NORMAL;
}

static int mask_boundary_node(const struct hecmwST_local_mesh *global_mesh,
                              char *node_flag, const char *elem_flag) {
  int node;
  int i, j;

  for (i = 0; i < global_mesh->n_elem; i++) {
    if (EVAL_BIT(elem_flag[i], BOUNDARY)) {
      for (j = global_mesh->elem_node_index[i];
           j < global_mesh->elem_node_index[i + 1]; j++) {
        node = global_mesh->elem_node_item[j];
        MASK_BIT(node_flag[node - 1], OVERLAP);
        MASK_BIT(node_flag[node - 1], BOUNDARY);
      }
    }
  }

  return RTC_NORMAL;
}

/*K. Inagaki */
static int mask_boundary_node_mod(const struct hecmwST_local_mesh *global_mesh,
                                  char *node_flag, char *elem_flag,
                                  int domain) {
  int i, node;

  for (i = 0; i < n_bnd_nlist[2 * domain + 1]; i++) {
    node = bnd_nlist[domain][i];
    MASK_BIT(node_flag[node - 1], OVERLAP);
    MASK_BIT(node_flag[node - 1], BOUNDARY);
  }

  return RTC_NORMAL;
}

#if 0
/* For Additional overlap for explicit DOF elimination for MPC */
/* NO LONGER NEEDED because node-migration implemented */
static int mask_boundary_elem_with_slave(
    const struct hecmwST_local_mesh *global_mesh, const char *node_flag,
    char *elem_flag, int *added) {
  int node, evalsum;
  int i, j;

  *added = 0;

  for (i = 0; i < global_mesh->n_elem; i++) {
    if (EVAL_BIT(elem_flag[i], BOUNDARY)) continue;
    if (HECMW_is_etype_link(global_mesh->elem_type[i]))
      continue; /* skip link elements */

    evalsum = 0;
    for (j = global_mesh->elem_node_index[i];
         j < global_mesh->elem_node_index[i + 1]; j++) {
      node = global_mesh->elem_node_item[j];
      /* check if the node is on boundary and a slave having mpc-link across the
       * boundary */
      if (EVAL_BIT(node_flag[node - 1], BOUNDARY) &&
          EVAL_BIT(node_flag[node - 1], MASK) &&
          EVAL_BIT(node_flag[node - 1], MARK)) {
        evalsum++;
      }
    }

    if (evalsum) {
      MASK_BIT(elem_flag[i], OVERLAP);
      MASK_BIT(elem_flag[i], BOUNDARY);
      (*added)++;
    }
  }

  return RTC_NORMAL;
}

static int mask_boundary_link_elem_with_slave(
    const struct hecmwST_local_mesh *global_mesh, const char *node_flag,
    char *elem_flag, int *added) {
  int node, evalsum;
  int i, j;

  *added = 0;

  for (i = 0; i < global_mesh->n_elem; i++) {
    if (EVAL_BIT(elem_flag[i], BOUNDARY)) continue;
    if (!HECMW_is_etype_link(global_mesh->elem_type[i]))
      continue; /* check only link elements */

    evalsum = 0;
    for (j = global_mesh->elem_node_index[i];
         j < global_mesh->elem_node_index[i + 1]; j++) {
      node = global_mesh->elem_node_item[j];
      /* check if the node is on boundary and a slave */
      if (EVAL_BIT(node_flag[node - 1], BOUNDARY) &&
          EVAL_BIT(node_flag[node - 1], MASK)) {
        evalsum++;
      }
    }

    if (evalsum) {
      MASK_BIT(elem_flag[i], OVERLAP);
      MASK_BIT(elem_flag[i], BOUNDARY);
      (*added)++;
    }
  }

  return RTC_NORMAL;
}
#endif

static int mask_additional_overlap_elem(
    const struct hecmwST_local_mesh *global_mesh, const char *node_flag,
    char *elem_flag) {
  int node, evalsum;
  int i, j;

  for (i = 0; i < global_mesh->n_elem; i++) {
    evalsum = 0;
    for (j = global_mesh->elem_node_index[i];
         j < global_mesh->elem_node_index[i + 1]; j++) {
      node = global_mesh->elem_node_item[j];
      evalsum += (EVAL_BIT(node_flag[node - 1], BOUNDARY));
    }

    if (evalsum) {
      MASK_BIT(elem_flag[i], OVERLAP);
      MASK_BIT(elem_flag[i], BOUNDARY);
    }
  }

  return RTC_NORMAL;
}

static int mask_contact_slave_surf(const struct hecmwST_local_mesh *global_mesh,
                                   char *elem_flag, char *node_flag) {
  int i, j, k;
  int elem, node, selem;
  int evalsum, evalsum2;
  int master_gid, slave_gid;
  int jstart, jend;
  struct hecmwST_contact_pair *cp;
  struct hecmwST_surf_grp *sgrp;
  struct hecmwST_node_grp *ngrp;
  struct hecmwST_elem_grp *egrp;

  cp   = global_mesh->contact_pair;
  sgrp = global_mesh->surf_group;
  ngrp = global_mesh->node_group;
  egrp = global_mesh->elem_group;

  for (i = 0; i < cp->n_pair; i++) {
    switch (cp->type[i]) {
    case HECMW_CONTACT_TYPE_NODE_SURF:
      /* if any elem of master surf is internal */
      evalsum    = 0;
      master_gid = cp->master_grp_id[i];
      jstart     = sgrp->grp_index[master_gid - 1];
      jend       = sgrp->grp_index[master_gid];
      for (j = jstart; j < jend; j++) {
        elem = sgrp->grp_item[j * 2];
        if (EVAL_BIT(elem_flag[elem - 1], INTERNAL)) {
          evalsum++;
          break;
        }
      }
      if (evalsum) {
        /* mask all external slave nodes as BOUNDARY (but not OVERLAP) */
        slave_gid = cp->slave_grp_id[i];
        jstart    = ngrp->grp_index[slave_gid - 1];
        jend      = ngrp->grp_index[slave_gid];
        for (j = jstart; j < jend; j++) {
          node = ngrp->grp_item[j];
          if (!EVAL_BIT(node_flag[node - 1], INTERNAL)) {
            MASK_BIT(node_flag[node - 1], BOUNDARY);
          }
        }
      }
      /* if any elem of master surf is external */
      evalsum    = 0;
      master_gid = cp->master_grp_id[i];
      jstart     = sgrp->grp_index[master_gid - 1];
      jend       = sgrp->grp_index[master_gid];
      for (j = jstart; j < jend; j++) {
        elem = sgrp->grp_item[j * 2];
        if (!EVAL_BIT(elem_flag[elem - 1], INTERNAL)) {
          evalsum++;
          break;
        }
      }
      if (evalsum) {
        /* mask all internal slave nodes as BOUNDARY (but not OVERLAP) */
        slave_gid = cp->slave_grp_id[i];
        jstart    = ngrp->grp_index[slave_gid - 1];
        jend      = ngrp->grp_index[slave_gid];
        for (j = jstart; j < jend; j++) {
          node = ngrp->grp_item[j];
          if (EVAL_BIT(node_flag[node - 1], INTERNAL)) {
            MASK_BIT(node_flag[node - 1], BOUNDARY);
          }
        }
      }
      break;

    case HECMW_CONTACT_TYPE_SURF_SURF:
      /* if any elem of master surf is internal or boundary */
      evalsum    = 0;
      master_gid = cp->master_grp_id[i];
      jstart     = sgrp->grp_index[master_gid - 1];
      jend       = sgrp->grp_index[master_gid];
      for (j = jstart; j < jend; j++) {
        elem = sgrp->grp_item[j * 2];
        if (EVAL_BIT(elem_flag[elem - 1], INTERNAL)
            || EVAL_BIT(elem_flag[elem - 1], BOUNDARY)) {
          evalsum++;
          break;
        }
      }
      if (evalsum) {
        /* mask all external slave elems/nodes as BOUNDARY (but not OVERLAP) */
        slave_gid = cp->slave_grp_id[i];
        jstart    = sgrp->grp_index[slave_gid - 1];
        jend      = sgrp->grp_index[slave_gid];
        for (j = jstart; j < jend; j++) {
          selem = sgrp->grp_item[j * 2];
          if (!EVAL_BIT(elem_flag[selem - 1], INTERNAL)) {
            MASK_BIT(elem_flag[selem - 1], BOUNDARY);
            for (k = global_mesh->elem_node_index[selem - 1];
                 k < global_mesh->elem_node_index[selem]; k++) {
              node = global_mesh->elem_node_item[k];
              MASK_BIT(node_flag[node - 1], BOUNDARY);
            }
          }
        }
      }
      /* if any elem of master surf is external or boundary */
      evalsum    = 0;
      master_gid = cp->master_grp_id[i];
      jstart     = sgrp->grp_index[master_gid - 1];
      jend       = sgrp->grp_index[master_gid];
      for (j = jstart; j < jend; j++) {
        elem = sgrp->grp_item[j * 2];
        if (!EVAL_BIT(elem_flag[elem - 1], INTERNAL)
            || EVAL_BIT(elem_flag[elem - 1], BOUNDARY)) {
          evalsum++;
          break;
        }
      }
      if (evalsum) {
        /* mask all internal slave nodes as BOUNDARY (but not OVERLAP) */
        slave_gid = cp->slave_grp_id[i];
        jstart    = sgrp->grp_index[slave_gid - 1];
        jend      = sgrp->grp_index[slave_gid];
        for (j = jstart; j < jend; j++) {
          evalsum2 = 0;
          selem = sgrp->grp_item[j * 2];
          for (k = global_mesh->elem_node_index[selem - 1];
               k < global_mesh->elem_node_index[selem]; k++) {
            node = global_mesh->elem_node_item[k];
            if (EVAL_BIT(node_flag[node - 1], INTERNAL)) {
              evalsum2++;
              break;
            }
          }
          if (evalsum2) {
            MASK_BIT(elem_flag[selem - 1], BOUNDARY);
            for (k = global_mesh->elem_node_index[selem - 1];
                 k < global_mesh->elem_node_index[selem]; k++) {
              node = global_mesh->elem_node_item[k];
              MASK_BIT(node_flag[node - 1], BOUNDARY);
            }
          }
        }
      }
      break;

    case HECMW_CONTACT_TYPE_NODE_ELEM:
      /* if any elem of master surf is internal */
      evalsum    = 0;
      master_gid = cp->master_grp_id[i];
      jstart     = egrp->grp_index[master_gid - 1];
      jend       = egrp->grp_index[master_gid];
      for (j = jstart; j < jend; j++) {
        elem = egrp->grp_item[j];
        if (EVAL_BIT(elem_flag[elem - 1], INTERNAL)) {
          evalsum++;
          break;
        }
      }
      if (evalsum) {
        /* mask all external slave nodes as BOUNDARY (but not OVERLAP) */
        slave_gid = cp->slave_grp_id[i];
        jstart    = ngrp->grp_index[slave_gid - 1];
        jend      = ngrp->grp_index[slave_gid];
        for (j = jstart; j < jend; j++) {
          node = ngrp->grp_item[j];
          if (!EVAL_BIT(node_flag[node - 1], INTERNAL)) {
            MASK_BIT(node_flag[node - 1], BOUNDARY);
          }
        }
      }
      /* if any elem of master surf is external */
      evalsum    = 0;
      master_gid = cp->master_grp_id[i];
      jstart     = egrp->grp_index[master_gid - 1];
      jend       = egrp->grp_index[master_gid];
      for (j = jstart; j < jend; j++) {
        elem = egrp->grp_item[j];
        if (!EVAL_BIT(elem_flag[elem - 1], INTERNAL)) {
          evalsum++;
          break;
        }
      }
      if (evalsum) {
        /* mask all internal slave nodes as BOUNDARY (but not OVERLAP) */
        slave_gid = cp->slave_grp_id[i];
        jstart    = ngrp->grp_index[slave_gid - 1];
        jend      = ngrp->grp_index[slave_gid];
        for (j = jstart; j < jend; j++) {
          node = ngrp->grp_item[j];
          if (EVAL_BIT(node_flag[node - 1], INTERNAL)) {
            MASK_BIT(node_flag[node - 1], BOUNDARY);
          }
        }
      }
      break;

    default:
      return RTC_ERROR;
    }
  }

  return RTC_NORMAL;
}

static int mask_mesh_status_nb(const struct hecmwST_local_mesh *global_mesh,
                               char *node_flag, char *elem_flag,
                               int current_domain) {
  int rtc;
  int i;

  rtc = mask_node_by_domain(global_mesh, node_flag, current_domain);
  if (rtc != RTC_NORMAL) goto error;

  rtc = mask_elem_by_domain_mod(elem_flag, current_domain);
  if (rtc != RTC_NORMAL) goto error;

  rtc = mask_overlap_elem(elem_flag, current_domain);
  if (rtc != RTC_NORMAL) goto error;

  rtc =
      mask_boundary_node_mod(global_mesh, node_flag, elem_flag, current_domain);
  if (rtc != RTC_NORMAL) goto error;

#if 0
  /* Additional overlap for explicit DOF elimination for MPC */
  /* NO LONGER NEEDED because node-migration implemented */
  if (global_mesh->mpc->n_mpc > 0) {
    int added = 0;

    rtc = mask_slave_node(global_mesh, node_flag, current_domain);
    if (rtc != RTC_NORMAL) goto error;

    rtc = mask_boundary_elem_with_slave(global_mesh, node_flag, elem_flag,
                                        &added);
    if (rtc != RTC_NORMAL) goto error;

    if (added > 0) {
      rtc = mask_boundary_node(global_mesh, node_flag, elem_flag);
      if (rtc != RTC_NORMAL) goto error;
    }

    added = 0;
    rtc = mask_boundary_link_elem_with_slave(global_mesh, node_flag, elem_flag,
                                             &added);
    if (rtc != RTC_NORMAL) goto error;

    if (added > 0) {
      rtc = mask_boundary_node(global_mesh, node_flag, elem_flag);
      if (rtc != RTC_NORMAL) goto error;
    }

    for (i = 0; i < global_mesh->n_node; i++) {
      CLEAR_BIT(node_flag[i], MASK);
      CLEAR_BIT(node_flag[i], MARK);
    }
  }
#endif

  for (i = 1; i < global_mesh->hecmw_flag_partdepth; i++) {
    rtc = mask_additional_overlap_elem(global_mesh, node_flag, elem_flag);
    if (rtc != RTC_NORMAL) goto error;

    rtc = mask_boundary_node(global_mesh, node_flag, elem_flag);
    if (rtc != RTC_NORMAL) goto error;
  }

  if (global_mesh->contact_pair->n_pair > 0) {
    rtc = mask_contact_slave_surf(global_mesh, elem_flag, node_flag);
    if (rtc != RTC_NORMAL) goto error;
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * - - - - - - - - - */

static int mask_overlap_node_mark(const struct hecmwST_local_mesh *global_mesh,
                                  char *node_flag, const char *elem_flag) {
  int node;
  int i, j;

  for (i = 0; i < global_mesh->n_elem; i++) {
    if (EVAL_BIT(elem_flag[i], INTERNAL)) {
      for (j = global_mesh->elem_node_index[i];
           j < global_mesh->elem_node_index[i + 1]; j++) {
        node = global_mesh->elem_node_item[j];
        MASK_BIT(node_flag[node - 1], MARK);
      }

    } else {
      for (j = global_mesh->elem_node_index[i];
           j < global_mesh->elem_node_index[i + 1]; j++) {
        node = global_mesh->elem_node_item[j];
        MASK_BIT(node_flag[node - 1], MASK);
      }
    }
  }

  return RTC_NORMAL;
}

static int mask_overlap_node_inner(const struct hecmwST_local_mesh *global_mesh,
                                   char *node_flag) {
  int i;

  for (i = 0; i < global_mesh->n_node; i++) {
    if (EVAL_BIT(node_flag[i], MARK) && EVAL_BIT(node_flag[i], MASK)) {
      MASK_BIT(node_flag[i], OVERLAP);
      MASK_BIT(node_flag[i], BOUNDARY);
    }
  }

  return RTC_NORMAL;
}

static int mask_overlap_node(const struct hecmwST_local_mesh *global_mesh,
                             char *node_flag, const char *elem_flag) {
  int rtc;
  int i;

  rtc = mask_overlap_node_mark(global_mesh, node_flag, elem_flag);
  if (rtc != RTC_NORMAL) goto error;

  rtc = mask_overlap_node_inner(global_mesh, node_flag);
  if (rtc != RTC_NORMAL) goto error;

  for (i = 0; i < global_mesh->n_node; i++) {
    CLEAR_BIT(node_flag[i], MASK);
    CLEAR_BIT(node_flag[i], MARK);
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int mask_boundary_elem(const struct hecmwST_local_mesh *global_mesh,
                              const char *node_flag, char *elem_flag) {
  int node, evalsum;
  int i, j;

  for (i = 0; i < global_mesh->n_elem; i++) {
    evalsum = 0;
    for (j = global_mesh->elem_node_index[i];
         j < global_mesh->elem_node_index[i + 1]; j++) {
      node = global_mesh->elem_node_item[j];
      if (EVAL_BIT(node_flag[node - 1], BOUNDARY)) evalsum++;
    }

    if (evalsum) {
      MASK_BIT(elem_flag[i], OVERLAP);
      MASK_BIT(elem_flag[i], BOUNDARY);
    }
  }

  return RTC_NORMAL;
}

static int mask_mesh_status_eb(const struct hecmwST_local_mesh *global_mesh,
                               char *node_flag, char *elem_flag,
                               int current_domain) {
  int rtc;
  int i;

  for (i = 0; i < global_mesh->n_node; i++) {
    CLEAR_BIT(node_flag[i], INTERNAL);
    CLEAR_BIT(node_flag[i], EXTERNAL);
    CLEAR_BIT(node_flag[i], BOUNDARY);
  }
  for (i = 0; i < global_mesh->n_elem; i++) {
    CLEAR_BIT(elem_flag[i], INTERNAL);
    CLEAR_BIT(elem_flag[i], EXTERNAL);
    CLEAR_BIT(elem_flag[i], BOUNDARY);
  }

  rtc = mask_node_by_domain(global_mesh, node_flag, current_domain);
  if (rtc != RTC_NORMAL) goto error;

  rtc = mask_elem_by_domain(global_mesh, elem_flag, current_domain);
  if (rtc != RTC_NORMAL) goto error;

  rtc = mask_overlap_node(global_mesh, node_flag, elem_flag);
  if (rtc != RTC_NORMAL) goto error;

  rtc = mask_boundary_elem(global_mesh, node_flag, elem_flag);
  if (rtc != RTC_NORMAL) goto error;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*------------------------------------------------------------------------------------------------*/

static int mask_neighbor_domain_nb(const struct hecmwST_local_mesh *global_mesh,
                                   const char *node_flag, char *domain_flag) {
  int i;

  for (i = 0; i < global_mesh->n_node; i++) {
    if (!EVAL_BIT(node_flag[i], INTERNAL) && EVAL_BIT(node_flag[i], BOUNDARY)) {
      MASK_BIT(domain_flag[global_mesh->node_ID[2 * i + 1]], MASK);
    }
  }

  return RTC_NORMAL;
}

/*K. Inagaki */
static int mask_neighbor_domain_nb_mod(
    const struct hecmwST_local_mesh *global_mesh, const char *node_flag,
    char *domain_flag, int domain) {
  int i, node;

  for (i = n_bnd_nlist[2 * domain]; i < n_bnd_nlist[2 * domain + 1]; i++) {
    node = bnd_nlist[domain][i];
    MASK_BIT(domain_flag[global_mesh->node_ID[2 * node - 1]], MASK);
  }

  return RTC_NORMAL;
}

static int mask_neighbor_domain_nb_contact(
    const struct hecmwST_local_mesh *global_mesh, const char *node_flag,
    const char *elem_flag, char *domain_flag) {
  int i, j, k;
  int elem, node, selem;
  int evalsum;
  int master_gid, slave_gid;
  int jstart, jend;
  struct hecmwST_contact_pair *cp;
  struct hecmwST_surf_grp *sgrp;
  struct hecmwST_node_grp *ngrp;
  struct hecmwST_elem_grp *egrp;

  cp   = global_mesh->contact_pair;
  sgrp = global_mesh->surf_group;
  ngrp = global_mesh->node_group;
  egrp = global_mesh->elem_group;

  for (i = 0; i < cp->n_pair; i++) {
    /* if any slave node is internal */
    evalsum = 0;
    switch (cp->type[i]) {
    case HECMW_CONTACT_TYPE_NODE_SURF:
      slave_gid = cp->slave_grp_id[i];
      jstart    = ngrp->grp_index[slave_gid - 1];
      jend      = ngrp->grp_index[slave_gid];
      for (j = jstart; j < jend; j++) {
        node = ngrp->grp_item[j];
        if (EVAL_BIT(node_flag[node - 1], INTERNAL)) {
          evalsum++;
          break;
        }
      }
      break;
    case HECMW_CONTACT_TYPE_SURF_SURF:
      slave_gid = cp->slave_grp_id[i];
      jstart    = sgrp->grp_index[slave_gid - 1];
      jend      = sgrp->grp_index[slave_gid];
      for (j = jstart; j < jend; j++) {
        selem = sgrp->grp_item[j * 2];
        for (k = global_mesh->elem_node_index[selem - 1];
             k < global_mesh->elem_node_index[selem]; k++) {
          node = global_mesh->elem_node_item[k];
          if (EVAL_BIT(node_flag[node - 1], INTERNAL)) {
            evalsum++;
            break;
          }
        }
        if (evalsum) break;
      }
      break;
    case HECMW_CONTACT_TYPE_NODE_ELEM:
      slave_gid = cp->slave_grp_id[i];
      jstart    = ngrp->grp_index[slave_gid - 1];
      jend      = ngrp->grp_index[slave_gid];
      for (j = jstart; j < jend; j++) {
        node = ngrp->grp_item[j];
        if (EVAL_BIT(node_flag[node - 1], INTERNAL)) {
          evalsum++;
          break;
        }
      }
      break;
    default:
      return RTC_ERROR;
    }
    /* the domain to which elems of the master surf belong is neighbor */
    if (evalsum) {
      master_gid = cp->master_grp_id[i];
      if( cp->type[i] == HECMW_CONTACT_TYPE_NODE_ELEM ) {
        jstart     = egrp->grp_index[master_gid - 1];
        jend       = egrp->grp_index[master_gid];
        for (j = jstart; j < jend; j++) {
          elem = egrp->grp_item[j];
          if (!EVAL_BIT(elem_flag[elem - 1], INTERNAL)) {
            MASK_BIT(domain_flag[global_mesh->elem_ID[2 * (elem - 1) + 1]], MASK);
          }
        }
      } else {
        jstart     = sgrp->grp_index[master_gid - 1];
        jend       = sgrp->grp_index[master_gid];
        for (j = jstart; j < jend; j++) {
          elem = sgrp->grp_item[j * 2];
          if (!EVAL_BIT(elem_flag[elem - 1], INTERNAL)) {
            MASK_BIT(domain_flag[global_mesh->elem_ID[2 * (elem - 1) + 1]], MASK);
          }
        }
      }
    }
  }

  return RTC_NORMAL;
}

static int mask_neighbor_domain_eb(const struct hecmwST_local_mesh *global_mesh,
                                   const char *elem_flag, char *domain_flag) {
  int i;

  for (i = 0; i < global_mesh->n_elem; i++) {
    if (EVAL_BIT(elem_flag[i], EXTERNAL) && EVAL_BIT(elem_flag[i], BOUNDARY)) {
      MASK_BIT(domain_flag[global_mesh->elem_ID[2 * i + 1]], MASK);
    }
  }

  return RTC_NORMAL;
}

static int count_neighbor_domain(const struct hecmwST_local_mesh *global_mesh,
                                 const char *domain_flag) {
  int counter;
  int i;

  for (counter = 0, i = 0; i < global_mesh->n_subdomain; i++) {
    if (EVAL_BIT(domain_flag[i], MASK)) counter++;
  }

  return counter;
}

static int set_neighbor_domain(const struct hecmwST_local_mesh *global_mesh,
                               struct hecmwST_local_mesh *local_mesh,
                               const char *domain_flag) {
  int counter;
  int i;

  for (counter = 0, i = 0; i < global_mesh->n_subdomain; i++) {
    if (EVAL_BIT(domain_flag[i], MASK)) {
      local_mesh->neighbor_pe[counter++] = i;
    }
  }

  return counter;
}

static int create_neighbor_info(const struct hecmwST_local_mesh *global_mesh,
                                struct hecmwST_local_mesh *local_mesh,
                                char *node_flag, char *elem_flag,
                                int current_domain) {
  int rtc;
  char *domain_flag = NULL;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(node_flag);
  HECMW_assert(elem_flag);

  HECMW_log(HECMW_LOG_DEBUG,
            "Starting creation of neighboring domain information...");

  local_mesh->n_neighbor_pe = 0;
  local_mesh->neighbor_pe   = NULL;

  domain_flag = (char *)HECMW_calloc(global_mesh->n_subdomain, sizeof(char));
  if (domain_flag == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  switch (global_mesh->hecmw_flag_parttype) {
    case HECMW_FLAG_PARTTYPE_NODEBASED: /* for node-based partitioning */
      rtc = mask_mesh_status_nb(global_mesh, node_flag, elem_flag,
                                current_domain);
      if (rtc != RTC_NORMAL) goto error;

      if (is_spdup_available(global_mesh)) {
        rtc = mask_neighbor_domain_nb_mod(global_mesh, node_flag, domain_flag,
                                          current_domain);
      } else {
        rtc = mask_neighbor_domain_nb(global_mesh, node_flag, domain_flag);
      }
      if (rtc != RTC_NORMAL) goto error;

      rtc = mask_neighbor_domain_nb_contact(global_mesh, node_flag, elem_flag,
                                            domain_flag);
      if (rtc != RTC_NORMAL) goto error;

      break;

    case HECMW_FLAG_PARTTYPE_ELEMBASED: /* for element-based partitioning */
      rtc = mask_mesh_status_eb(global_mesh, node_flag, elem_flag,
                                current_domain);
      if (rtc != RTC_NORMAL) goto error;

      rtc = mask_neighbor_domain_eb(global_mesh, elem_flag, domain_flag);
      if (rtc != RTC_NORMAL) goto error;

      break;

    default:
      HECMW_set_error(HECMW_PART_E_INVALID_PTYPE, "");
      goto error;
  }

  local_mesh->n_neighbor_pe = count_neighbor_domain(global_mesh, domain_flag);
  if (local_mesh->n_neighbor_pe < 0) {
    HECMW_set_error(HECMW_PART_E_NNEIGHBORPE_LOWER, "");
    goto error;
  }

  if (local_mesh->n_neighbor_pe == 0) {
    local_mesh->neighbor_pe = NULL;
    HECMW_free(domain_flag);
    return RTC_NORMAL;
  }

  local_mesh->neighbor_pe =
      (int *)HECMW_malloc(sizeof(int) * local_mesh->n_neighbor_pe);
  if (local_mesh->neighbor_pe == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  rtc = set_neighbor_domain(global_mesh, local_mesh, domain_flag);
  HECMW_assert(rtc == local_mesh->n_neighbor_pe);

  HECMW_free(domain_flag);

  HECMW_log(HECMW_LOG_DEBUG, "Creation of neighboring domain information done");

  return RTC_NORMAL;

error:
  HECMW_free(domain_flag);
  HECMW_free(local_mesh->neighbor_pe);
  local_mesh->n_neighbor_pe = 0;
  local_mesh->neighbor_pe   = NULL;

  return RTC_ERROR;
}

/*================================================================================================*/

static int mask_comm_node(const struct hecmwST_local_mesh *global_mesh,
                          char *node_flag_current, char *node_flag_neighbor) {
  int i;

  for (i = 0; i < global_mesh->n_node; i++) {
    if (EVAL_BIT(node_flag_current[i], BOUNDARY) &&
        EVAL_BIT(node_flag_neighbor[i], BOUNDARY)) {
      MASK_BIT(node_flag_current[i], MASK);
    }
  }

  return RTC_NORMAL;
}

/*K. Inagaki */
static int mask_comm_node_mod(const struct hecmwST_local_mesh *global_mesh,
                              char *node_flag_current, char *node_flag_neighbor,
                              int current_domain) {
  int i, node;

  for (i = 0; i < n_bnd_nlist[2 * current_domain + 1]; i++) {
    node = bnd_nlist[current_domain][i];
    if (EVAL_BIT(node_flag_neighbor[node - 1], BOUNDARY)) {
      MASK_BIT(node_flag_current[node - 1], MASK);
    }
  }

  return RTC_NORMAL;
}

static int mask_comm_elem(const struct hecmwST_local_mesh *global_mesh,
                          char *elem_flag_current, char *elem_flag_neighbor) {
  int i;

  for (i = 0; i < global_mesh->n_elem; i++) {
    if (EVAL_BIT(elem_flag_current[i], BOUNDARY) &&
        EVAL_BIT(elem_flag_neighbor[i], BOUNDARY)) {
      MASK_BIT(elem_flag_current[i], MASK);
    }
  }

  return RTC_NORMAL;
}

/*K. Inagaki */
static int mask_comm_elem_mod(const struct hecmwST_local_mesh *global_mesh,
                              char *elem_flag_current, char *elem_flag_neighbor,
                              int current_domain) {
  int i, elem;

  for (i = 0; i < n_bnd_elist[2 * current_domain + 1]; i++) {
    elem = bnd_elist[current_domain][i];
    if (EVAL_BIT(elem_flag_neighbor[elem - 1], BOUNDARY)) {
      MASK_BIT(elem_flag_current[elem - 1], MASK);
    }
  }

  return RTC_NORMAL;
}

/*K. Inagaki */
static int count_masked_comm_node(const struct hecmwST_local_mesh *global_mesh,
                                  const char *node_flag, int domain) {
  int counter;
  int i, node;

  for (counter = 0, i = 0; i < n_int_nlist[domain]; i++) {
    node = int_nlist[domain][i];
    if (EVAL_BIT(node_flag[node - 1], MASK)) counter++;
  }

  return counter;
}

static int count_masked_comm_elem(const struct hecmwST_local_mesh *global_mesh,
                                  const char *elem_flag, int domain) {
  int counter;
  int i;

  for (counter = 0, i = 0; i < global_mesh->n_elem; i++) {
    if (EVAL_BIT(elem_flag[i], MASK) &&
        global_mesh->elem_ID[2 * i + 1] == domain)
      counter++;
  }

  return counter;
}

static int count_masked_shared_node(
    const struct hecmwST_local_mesh *global_mesh, const char *node_flag) {
  int counter;
  int i;

  for (counter = 0, i = 0; i < global_mesh->n_node; i++) {
    if (EVAL_BIT(node_flag[i], MASK)) counter++;
  }

  return counter;
}

static int count_masked_shared_elem(
    const struct hecmwST_local_mesh *global_mesh, const char *elem_flag) {
  int counter;
  int i;

  for (counter = 0, i = 0; i < global_mesh->n_elem; i++) {
    if (EVAL_BIT(elem_flag[i], MASK)) counter++;
  }

  return counter;
}

/*K. Inagaki */
static int count_masked_shared_elem_mod(
    const struct hecmwST_local_mesh *global_mesh, const char *elem_flag,
    int domain) {
  int counter;
  int i, elem;

  for (counter = 0, i = 0; i < n_bnd_elist[2 * domain + 1]; i++) {
    elem = bnd_elist[domain][i];
    if (EVAL_BIT(elem_flag[elem - 1], MASK)) counter++;
  }

  return counter;
}

/*K. Inagaki */
static int create_comm_node_pre(const struct hecmwST_local_mesh *global_mesh,
                                const char *node_flag, int **comm_node,
                                int neighbor_idx, int domain) {
  int counter;
  int i, node;

  for (counter = 0, i = 0; i < n_int_nlist[domain]; i++) {
    node = int_nlist[domain][i];
    if (EVAL_BIT(node_flag[node - 1], MASK)) {
      comm_node[neighbor_idx][counter++] = node;
    }
  }

  return counter;
}

static int create_comm_elem_pre(const struct hecmwST_local_mesh *global_mesh,
                                const char *elem_flag, int **comm_elem,
                                int neighbor_idx, int domain) {
  int counter;
  int i;

  for (counter = 0, i = 0; i < global_mesh->n_elem; i++) {
    if (EVAL_BIT(elem_flag[i], MASK) &&
        global_mesh->elem_ID[2 * i + 1] == domain) {
      comm_elem[neighbor_idx][counter++] = i + 1;
    }
  }

  return counter;
}

static int create_shared_node_pre(const struct hecmwST_local_mesh *global_mesh,
                                  const char *node_flag, int **shared_node,
                                  int neighbor_idx) {
  int counter;
  int i;

  for (counter = 0, i = 0; i < global_mesh->n_node; i++) {
    if (EVAL_BIT(node_flag[i], MASK)) {
      shared_node[neighbor_idx][counter++] = i + 1;
    }
  }

  return counter;
}

static int create_shared_elem_pre(const struct hecmwST_local_mesh *global_mesh,
                                  const char *elem_flag, int **shared_elem,
                                  int neighbor_idx) {
  int counter;
  int i;

  for (counter = 0, i = 0; i < global_mesh->n_elem; i++) {
    if (EVAL_BIT(elem_flag[i], MASK)) {
      shared_elem[neighbor_idx][counter++] = i + 1;
    }
  }

  return counter;
}

/*K. Inagaki */
static int create_shared_elem_pre_mod(
    const struct hecmwST_local_mesh *global_mesh, const char *elem_flag,
    int **shared_elem, int neighbor_idx, int neighbor_domain) {
  int counter;
  int i, idx1, idx2, elem1, elem2, n_bnd, n_out, maxe;

  n_bnd = n_bnd_elist[2 * neighbor_domain];
  n_out =
      n_bnd_elist[2 * neighbor_domain + 1] - n_bnd_elist[2 * neighbor_domain];
  maxe = global_mesh->n_elem + 1;

  elem1 = (n_bnd == 0) ? maxe : bnd_elist[neighbor_domain][0];
  elem2 = (n_out == 0) ? maxe : bnd_elist[neighbor_domain][n_bnd];
  for (counter = 0, idx1 = 0, idx2 = 0, i = 0; i < n_bnd + n_out; i++) {
    if (elem1 < elem2) {
      if (EVAL_BIT(elem_flag[elem1 - 1], MASK)) {
        shared_elem[neighbor_idx][counter++] = elem1;
      }
      idx1++;
      elem1 = (idx1 == n_bnd) ? maxe : bnd_elist[neighbor_domain][idx1];
    } else {
      if (EVAL_BIT(elem_flag[elem2 - 1], MASK)) {
        shared_elem[neighbor_idx][counter++] = elem2;
      }
      idx2++;
      elem2 = (idx2 == n_out) ? maxe : bnd_elist[neighbor_domain][idx2 + n_bnd];
    }
  }

  return counter;
}

static int create_comm_item(int n_neighbor_pe, int **comm_item_pre,
                            int *comm_index, int *comm_item) {
  int i, j, js, je;

  for (i = 0; i < n_neighbor_pe; i++) {
    js = comm_index[i];
    je = comm_index[i + 1];

    for (j = 0; j < je - js; j++) {
      comm_item[js + j] = comm_item_pre[i][j];
    }
  }

  return RTC_NORMAL;
}

/*------------------------------------------------------------------------------------------------*/

static int create_import_info_nb(const struct hecmwST_local_mesh *global_mesh,
                                 struct hecmwST_local_mesh *local_mesh,
                                 const char *node_flag, int **import_node,
                                 int neighbor_idx, int neighbor_domain) {
  int n_import_node, rtc;

  n_import_node =
      count_masked_comm_node(global_mesh, node_flag, neighbor_domain);
  HECMW_assert(n_import_node >= 0);

  local_mesh->import_index[neighbor_idx + 1] =
      local_mesh->import_index[neighbor_idx] + n_import_node;

  import_node[neighbor_idx] = (int *)HECMW_malloc(sizeof(int) * n_import_node);
  if (import_node[neighbor_idx] == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  rtc = create_comm_node_pre(global_mesh, node_flag, import_node, neighbor_idx,
                             neighbor_domain);
  HECMW_assert(rtc == n_import_node);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int create_export_info_nb(const struct hecmwST_local_mesh *global_mesh,
                                 struct hecmwST_local_mesh *local_mesh,
                                 const char *node_flag, int **export_node,
                                 int neighbor_idx, int current_domain,
                                 int neighbor_domain) {
  int n_export_node, rtc;

  n_export_node =
      count_masked_comm_node(global_mesh, node_flag, current_domain);
  HECMW_assert(n_export_node >= 0);

  local_mesh->export_index[neighbor_idx + 1] =
      local_mesh->export_index[neighbor_idx] + n_export_node;

  export_node[neighbor_idx] = (int *)HECMW_malloc(sizeof(int) * n_export_node);
  if (export_node[neighbor_idx] == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  rtc = create_comm_node_pre(global_mesh, node_flag, export_node, neighbor_idx,
                             current_domain);
  HECMW_assert(rtc == n_export_node);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int create_shared_info_nb(const struct hecmwST_local_mesh *global_mesh,
                                 struct hecmwST_local_mesh *local_mesh,
                                 const char *elem_flag, int **shared_elem,
                                 int neighbor_idx, int neighbor_domain) {
  int n_shared_elem, rtc;

  if (is_spdup_available(global_mesh)) {
    n_shared_elem =
        count_masked_shared_elem_mod(global_mesh, elem_flag, neighbor_domain);
  } else {
    n_shared_elem = count_masked_shared_elem(global_mesh, elem_flag);
  }

  HECMW_assert(n_shared_elem >= 0);

  local_mesh->shared_index[neighbor_idx + 1] =
      local_mesh->shared_index[neighbor_idx] + n_shared_elem;

  shared_elem[neighbor_idx] = (int *)HECMW_malloc(sizeof(int) * n_shared_elem);
  if (shared_elem[neighbor_idx] == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  if (is_spdup_available(global_mesh)) {
    rtc = create_shared_elem_pre_mod(global_mesh, elem_flag, shared_elem,
                                     neighbor_idx, neighbor_domain);
  } else {
    rtc = create_shared_elem_pre(global_mesh, elem_flag, shared_elem,
                                 neighbor_idx);
  }

  HECMW_assert(rtc == n_shared_elem);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int create_comm_info_nb(const struct hecmwST_local_mesh *global_mesh,
                               struct hecmwST_local_mesh *local_mesh,
                               char *node_flag, char *elem_flag,
                               char *node_flag_neighbor,
                               char *elem_flag_neighbor, int current_domain) {
  int **import_node = NULL;
  int **export_node = NULL;
  int **shared_elem = NULL;
  int neighbor_domain;
  int size;
  int rtc;
  int i, j;

  local_mesh->import_index = NULL;
  local_mesh->export_index = NULL;
  local_mesh->shared_index = NULL;
  local_mesh->import_item  = NULL;
  local_mesh->export_item  = NULL;
  local_mesh->shared_item  = NULL;

  import_node = (int **)HECMW_malloc(sizeof(int *) * local_mesh->n_neighbor_pe);
  if (import_node == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  } else {
    for (i = 0; i < local_mesh->n_neighbor_pe; i++) {
      import_node[i] = NULL;
    }
  }
  export_node = (int **)HECMW_malloc(sizeof(int *) * local_mesh->n_neighbor_pe);
  if (export_node == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  } else {
    for (i = 0; i < local_mesh->n_neighbor_pe; i++) {
      export_node[i] = NULL;
    }
  }
  shared_elem = (int **)HECMW_malloc(sizeof(int *) * local_mesh->n_neighbor_pe);
  if (shared_elem == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  } else {
    for (i = 0; i < local_mesh->n_neighbor_pe; i++) {
      shared_elem[i] = NULL;
    }
  }

  local_mesh->import_index =
      (int *)HECMW_calloc(local_mesh->n_neighbor_pe + 1, sizeof(int));
  if (local_mesh->import_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  local_mesh->export_index =
      (int *)HECMW_calloc(local_mesh->n_neighbor_pe + 1, sizeof(int));
  if (local_mesh->export_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  local_mesh->shared_index =
      (int *)HECMW_calloc(local_mesh->n_neighbor_pe + 1, sizeof(int));
  if (local_mesh->shared_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (i = 0; i < local_mesh->n_neighbor_pe; i++) {
    neighbor_domain = local_mesh->neighbor_pe[i];

    rtc = mask_mesh_status_nb(global_mesh, node_flag_neighbor,
                              elem_flag_neighbor, neighbor_domain);
    if (rtc != RTC_NORMAL) goto error;

    if (is_spdup_available(global_mesh)) {
      rtc = mask_comm_node_mod(global_mesh, node_flag, node_flag_neighbor,
                               current_domain);
    } else {
      rtc = mask_comm_node(global_mesh, node_flag, node_flag_neighbor);
    }

    if (rtc != RTC_NORMAL) goto error;

    if (is_spdup_available(global_mesh)) {
      rtc = mask_comm_elem_mod(global_mesh, elem_flag, elem_flag_neighbor,
                               current_domain);
    } else {
      rtc = mask_comm_elem(global_mesh, elem_flag, elem_flag_neighbor);
    }

    if (rtc != RTC_NORMAL) goto error;

    rtc = create_import_info_nb(global_mesh, local_mesh, node_flag, import_node,
                                i, neighbor_domain);
    if (rtc != RTC_NORMAL) goto error;

    rtc = create_export_info_nb(global_mesh, local_mesh, node_flag, export_node,
                                i, current_domain, neighbor_domain);
    if (rtc != RTC_NORMAL) goto error;

    rtc = create_shared_info_nb(global_mesh, local_mesh, elem_flag, shared_elem,
                                i, neighbor_domain);
    if (rtc != RTC_NORMAL) goto error;

    if (is_spdup_available(global_mesh)) {
      /*K. Inagaki */
      rtc = spdup_clear_IEB(node_flag_neighbor, elem_flag_neighbor,
                            neighbor_domain);
      if (rtc != RTC_NORMAL) goto error;

      rtc = spdup_clear_MMbnd(node_flag_neighbor, elem_flag_neighbor,
                              neighbor_domain);
      if (rtc != RTC_NORMAL) goto error;

      rtc = spdup_clear_MMbnd(node_flag, elem_flag, current_domain);
      if (rtc != RTC_NORMAL) goto error;
    } else {
      for (j = 0; j < global_mesh->n_node; j++) {
        CLEAR_MM(node_flag[j]);
      }
      for (j = 0; j < global_mesh->n_elem; j++) {
        CLEAR_MM(elem_flag[j]);
      }

      memset(node_flag_neighbor, 0, sizeof(char) * global_mesh->n_node);
      memset(elem_flag_neighbor, 0, sizeof(char) * global_mesh->n_elem);
    }
  }

  size = sizeof(int) * local_mesh->import_index[local_mesh->n_neighbor_pe];
  local_mesh->import_item = (int *)HECMW_malloc(size);
  if (local_mesh->import_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  rtc = create_comm_item(local_mesh->n_neighbor_pe, import_node,
                         local_mesh->import_index, local_mesh->import_item);
  if (rtc != RTC_NORMAL) goto error;

  for (i = 0; i < local_mesh->n_neighbor_pe; i++) {
    HECMW_free(import_node[i]);
  }
  HECMW_free(import_node);
  import_node = NULL;

  size = sizeof(int) * local_mesh->export_index[local_mesh->n_neighbor_pe];
  local_mesh->export_item = (int *)HECMW_malloc(size);
  if (local_mesh->export_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  rtc = create_comm_item(local_mesh->n_neighbor_pe, export_node,
                         local_mesh->export_index, local_mesh->export_item);
  if (rtc != RTC_NORMAL) goto error;

  for (i = 0; i < local_mesh->n_neighbor_pe; i++) {
    HECMW_free(export_node[i]);
  }
  HECMW_free(export_node);
  export_node = NULL;

  size = sizeof(int) * local_mesh->shared_index[local_mesh->n_neighbor_pe];
  local_mesh->shared_item = (int *)HECMW_malloc(size);
  if (local_mesh->shared_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  rtc = create_comm_item(local_mesh->n_neighbor_pe, shared_elem,
                         local_mesh->shared_index, local_mesh->shared_item);
  if (rtc != RTC_NORMAL) goto error;

  for (i = 0; i < local_mesh->n_neighbor_pe; i++) {
    HECMW_free(shared_elem[i]);
  }
  HECMW_free(shared_elem);
  shared_elem = NULL;

  return RTC_NORMAL;

error:
  if (import_node) {
    int i;
    for (i = 0; i < local_mesh->n_neighbor_pe; i++) {
      HECMW_free(import_node[i]);
    }
    HECMW_free(import_node);
  }
  if (export_node) {
    int i;
    for (i = 0; i < local_mesh->n_neighbor_pe; i++) {
      HECMW_free(export_node[i]);
    }
    HECMW_free(export_node);
  }
  if (shared_elem) {
    int i;
    for (i = 0; i < local_mesh->n_neighbor_pe; i++) {
      HECMW_free(shared_elem[i]);
    }
    HECMW_free(shared_elem);
  }

  HECMW_free(local_mesh->import_index);
  HECMW_free(local_mesh->export_index);
  HECMW_free(local_mesh->shared_index);
  HECMW_free(local_mesh->import_item);
  HECMW_free(local_mesh->export_item);
  HECMW_free(local_mesh->shared_item);

  local_mesh->import_index = NULL;
  local_mesh->export_index = NULL;
  local_mesh->shared_index = NULL;
  local_mesh->import_item  = NULL;
  local_mesh->export_item  = NULL;
  local_mesh->shared_item  = NULL;

  return RTC_ERROR;
}

/*------------------------------------------------------------------------------------------------*/

static int create_import_info_eb(const struct hecmwST_local_mesh *global_mesh,
                                 struct hecmwST_local_mesh *local_mesh,
                                 const char *elem_flag, int **import_elem,
                                 int neighbor_idx, int neighbor_domain) {
  int n_import_elem, rtc;

  n_import_elem =
      count_masked_comm_elem(global_mesh, elem_flag, neighbor_domain);
  HECMW_assert(n_import_elem >= 0);

  local_mesh->import_index[neighbor_idx + 1] =
      local_mesh->import_index[neighbor_idx] + n_import_elem;

  import_elem[neighbor_idx] = (int *)HECMW_malloc(sizeof(int) * n_import_elem);
  if (import_elem[neighbor_idx] == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  rtc = create_comm_elem_pre(global_mesh, elem_flag, import_elem, neighbor_idx,
                             neighbor_domain);
  HECMW_assert(rtc == n_import_elem);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int create_export_info_eb(const struct hecmwST_local_mesh *global_mesh,
                                 struct hecmwST_local_mesh *local_mesh,
                                 const char *elem_flag, int **export_elem,
                                 int neighbor_idx, int current_domain,
                                 int neighbor_domain) {
  int n_export_elem, rtc;

  n_export_elem =
      count_masked_comm_elem(global_mesh, elem_flag, current_domain);
  HECMW_assert(n_export_elem >= 0);

  local_mesh->export_index[neighbor_idx + 1] =
      local_mesh->export_index[neighbor_idx] + n_export_elem;

  export_elem[neighbor_idx] = (int *)HECMW_malloc(sizeof(int) * n_export_elem);
  if (export_elem[neighbor_idx] == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  rtc = create_comm_elem_pre(global_mesh, elem_flag, export_elem, neighbor_idx,
                             current_domain);
  HECMW_assert(rtc == n_export_elem);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int create_shared_info_eb(const struct hecmwST_local_mesh *global_mesh,
                                 struct hecmwST_local_mesh *local_mesh,
                                 const char *node_flag, int **shared_node,
                                 int neighbor_idx, int neighbor_domain) {
  int n_shared_node, rtc;

  n_shared_node = count_masked_shared_node(global_mesh, node_flag);
  HECMW_assert(n_shared_node >= 0);

  local_mesh->shared_index[neighbor_idx + 1] =
      local_mesh->shared_index[neighbor_idx] + n_shared_node;

  shared_node[neighbor_idx] = (int *)HECMW_malloc(sizeof(int) * n_shared_node);
  if (shared_node[neighbor_idx] == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  rtc =
      create_shared_node_pre(global_mesh, node_flag, shared_node, neighbor_idx);
  HECMW_assert(rtc == n_shared_node);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*------------------------------------------------------------------------------------------------*/

static int create_comm_info_eb(const struct hecmwST_local_mesh *global_mesh,
                               struct hecmwST_local_mesh *local_mesh,
                               char *node_flag, char *elem_flag,
                               char *node_flag_neighbor,
                               char *elem_flag_neighbor, int current_domain) {
  int **import_elem = NULL;
  int **export_elem = NULL;
  int **shared_node = NULL;
  int neighbor_domain;
  int size;
  int rtc;
  int i, j;

  /* allocation */
  local_mesh->import_index = NULL;
  local_mesh->export_index = NULL;
  local_mesh->shared_index = NULL;
  local_mesh->import_item  = NULL;
  local_mesh->export_item  = NULL;
  local_mesh->shared_item  = NULL;

  import_elem = (int **)HECMW_malloc(sizeof(int *) * local_mesh->n_neighbor_pe);
  if (import_elem == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  } else {
    for (i = 0; i < local_mesh->n_neighbor_pe; i++) {
      import_elem[i] = NULL;
    }
  }
  export_elem = (int **)HECMW_malloc(sizeof(int *) * local_mesh->n_neighbor_pe);
  if (export_elem == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  } else {
    for (i = 0; i < local_mesh->n_neighbor_pe; i++) {
      export_elem[i] = NULL;
    }
  }
  shared_node = (int **)HECMW_malloc(sizeof(int *) * local_mesh->n_neighbor_pe);
  if (shared_node == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  } else {
    for (i = 0; i < local_mesh->n_neighbor_pe; i++) {
      shared_node[i] = NULL;
    }
  }

  local_mesh->import_index =
      (int *)HECMW_calloc(local_mesh->n_neighbor_pe + 1, sizeof(int));
  if (local_mesh->import_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  local_mesh->export_index =
      (int *)HECMW_calloc(local_mesh->n_neighbor_pe + 1, sizeof(int));
  if (local_mesh->export_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  local_mesh->shared_index =
      (int *)HECMW_calloc(local_mesh->n_neighbor_pe + 1, sizeof(int));
  if (local_mesh->shared_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  /* create communication table */
  for (i = 0; i < local_mesh->n_neighbor_pe; i++) {
    neighbor_domain = local_mesh->neighbor_pe[i];

    for (j = 0; j < global_mesh->n_node; j++) {
      CLEAR_BIT(node_flag[j], MASK);
      CLEAR_BIT(node_flag[j], MARK);
    }
    for (j = 0; j < global_mesh->n_elem; j++) {
      CLEAR_BIT(elem_flag[j], MASK);
      CLEAR_BIT(elem_flag[j], MARK);
    }

    memset(node_flag_neighbor, 0, sizeof(char) * global_mesh->n_node);
    memset(elem_flag_neighbor, 0, sizeof(char) * global_mesh->n_elem);

    /* mask boundary node & element */
    rtc = mask_mesh_status_eb(global_mesh, node_flag_neighbor,
                              elem_flag_neighbor, neighbor_domain);
    if (rtc != RTC_NORMAL) goto error;

    rtc = mask_comm_node(global_mesh, node_flag, node_flag_neighbor);
    if (rtc != RTC_NORMAL) goto error;

    rtc = mask_comm_elem(global_mesh, elem_flag, elem_flag_neighbor);
    if (rtc != RTC_NORMAL) goto error;

    /* create import element information (preliminary) */
    rtc = create_import_info_eb(global_mesh, local_mesh, elem_flag, import_elem,
                                i, neighbor_domain);
    if (rtc != RTC_NORMAL) goto error;

    /* create export element information (preliminary) */
    rtc = create_export_info_eb(global_mesh, local_mesh, elem_flag, export_elem,
                                i, current_domain, neighbor_domain);
    if (rtc != RTC_NORMAL) goto error;

    /* create shared node information (preliminary) */
    rtc = create_shared_info_eb(global_mesh, local_mesh, node_flag, shared_node,
                                i, neighbor_domain);
    if (rtc != RTC_NORMAL) goto error;
  }

  /* create import element information */
  size = sizeof(int) * local_mesh->import_index[local_mesh->n_neighbor_pe];
  local_mesh->import_item = (int *)HECMW_malloc(size);
  if (local_mesh->import_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  rtc = create_comm_item(local_mesh->n_neighbor_pe, import_elem,
                         local_mesh->import_index, local_mesh->import_item);
  if (rtc != RTC_NORMAL) goto error;

  for (i = 0; i < local_mesh->n_neighbor_pe; i++) {
    HECMW_free(import_elem[i]);
  }
  HECMW_free(import_elem);
  import_elem = NULL;

  /* create export node information */
  size = sizeof(int) * local_mesh->export_index[local_mesh->n_neighbor_pe];
  local_mesh->export_item = (int *)HECMW_malloc(size);
  if (local_mesh->export_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  rtc = create_comm_item(local_mesh->n_neighbor_pe, export_elem,
                         local_mesh->export_index, local_mesh->export_item);
  if (rtc != RTC_NORMAL) goto error;

  for (i = 0; i < local_mesh->n_neighbor_pe; i++) {
    HECMW_free(export_elem[i]);
  }
  HECMW_free(export_elem);
  export_elem = NULL;

  /* create shared element information */
  size = sizeof(int) * local_mesh->shared_index[local_mesh->n_neighbor_pe];
  local_mesh->shared_item = (int *)HECMW_malloc(size);
  if (local_mesh->shared_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  rtc = create_comm_item(local_mesh->n_neighbor_pe, shared_node,
                         local_mesh->shared_index, local_mesh->shared_item);
  if (rtc != RTC_NORMAL) goto error;

  for (i = 0; i < local_mesh->n_neighbor_pe; i++) {
    HECMW_free(shared_node[i]);
  }
  HECMW_free(shared_node);
  shared_node = NULL;

  return RTC_NORMAL;

error:
  if (import_elem) {
    int i;
    for (i = 0; i < local_mesh->n_neighbor_pe; i++) {
      HECMW_free(import_elem[i]);
    }
    HECMW_free(import_elem);
  }
  if (export_elem) {
    int i;
    for (i = 0; i < local_mesh->n_neighbor_pe; i++) {
      HECMW_free(export_elem[i]);
    }
    HECMW_free(export_elem);
  }
  if (shared_node) {
    int i;
    for (i = 0; i < local_mesh->n_neighbor_pe; i++) {
      HECMW_free(shared_node[i]);
    }
    HECMW_free(shared_node);
  }
  HECMW_free(local_mesh->import_index);
  HECMW_free(local_mesh->export_index);
  HECMW_free(local_mesh->shared_index);
  HECMW_free(local_mesh->import_item);
  HECMW_free(local_mesh->export_item);
  HECMW_free(local_mesh->shared_item);

  local_mesh->import_index = NULL;
  local_mesh->export_index = NULL;
  local_mesh->shared_index = NULL;
  local_mesh->import_item  = NULL;
  local_mesh->export_item  = NULL;
  local_mesh->shared_item  = NULL;

  return RTC_ERROR;
}

/*================================================================================================*/

static int create_comm_info(const struct hecmwST_local_mesh *global_mesh,
                            struct hecmwST_local_mesh *local_mesh,
                            char *node_flag, char *elem_flag,
                            char *node_flag_neighbor, char *elem_flag_neighbor,
                            int current_domain) {
  int rtc;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(node_flag);
  HECMW_assert(elem_flag);

  HECMW_log(HECMW_LOG_DEBUG, "Starting creation of interface table...");

  switch (global_mesh->hecmw_flag_parttype) {
    case HECMW_FLAG_PARTTYPE_NODEBASED: /* for node-based partitioning */
      rtc = create_comm_info_nb(global_mesh, local_mesh, node_flag, elem_flag,
                                node_flag_neighbor, elem_flag_neighbor,
                                current_domain);
      if (rtc != RTC_NORMAL) goto error;

      break;

    case HECMW_FLAG_PARTTYPE_ELEMBASED: /* for element-based partitioning */
      rtc = create_comm_info_eb(global_mesh, local_mesh, node_flag, elem_flag,
                                node_flag_neighbor, elem_flag_neighbor,
                                current_domain);
      if (rtc != RTC_NORMAL) goto error;

      break;

    default:
      HECMW_set_error(HECMW_PART_E_INVALID_PTYPE, "");
      goto error;
  }

  HECMW_log(HECMW_LOG_DEBUG, "Creation of interface table done");

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*==================================================================================================

   create distributed mesh information

==================================================================================================*/

/*K. Inagaki */
static int set_node_global2local_internal(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh, int *node_global2local,
    const char *node_flag, int domain) {
  int counter;
  int i, node;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(node_global2local);
  HECMW_assert(node_flag);
  HECMW_assert(global_mesh->n_node > 0);

  for (counter = 0, i = 0; i < n_int_nlist[domain]; i++) {
    node                        = int_nlist[domain][i];
    node_global2local[node - 1] = ++counter;
  }
  local_mesh->nn_internal = counter;

  return RTC_NORMAL;
}

static int set_node_global2local_external(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh, int *node_global2local,
    const char *node_flag) {
  int counter;
  int i;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(node_global2local);
  HECMW_assert(node_flag);
  HECMW_assert(global_mesh->n_node > 0);

  /* ordinary external nodes are marked as BOUNDARY && OVERLAP */
  for (counter = local_mesh->nn_internal, i = 0; i < global_mesh->n_node; i++) {
    if (!EVAL_BIT(node_flag[i], INTERNAL) && EVAL_BIT(node_flag[i], BOUNDARY) &&
        EVAL_BIT(node_flag[i], OVERLAP)) {
      node_global2local[i] = ++counter;
    }
  }
  local_mesh->nn_middle = counter;

  /* added external contact slave nodes are marked as BOUNDARY but not OVERLAP
   */
  for (i = 0; i < global_mesh->n_node; i++) {
    if (!EVAL_BIT(node_flag[i], INTERNAL) && EVAL_BIT(node_flag[i], BOUNDARY) &&
        !EVAL_BIT(node_flag[i], OVERLAP)) {
      node_global2local[i] = ++counter;
    }
  }
  local_mesh->n_node       = counter;
  local_mesh->n_node_gross = counter;

  HECMW_assert(local_mesh->n_node > 0);

  return RTC_NORMAL;
}

/*K. Inagaki */
static int set_node_global2local_external_mod(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh, int *node_global2local,
    const char *node_flag, int domain) {
  int counter;
  int i, node;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(node_global2local);
  HECMW_assert(node_flag);
  HECMW_assert(global_mesh->n_node > 0);

  for (counter = local_mesh->nn_internal, i = n_bnd_nlist[2 * domain];
       i < n_bnd_nlist[2 * domain + 1]; i++) {
    node                        = bnd_nlist[domain][i];
    node_global2local[node - 1] = ++counter;
  }
  local_mesh->n_node       = counter;
  local_mesh->n_node_gross = counter;

  HECMW_assert(local_mesh->n_node > 0);

  return RTC_NORMAL;
}

static int set_node_global2local_all(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh, int *node_global2local,
    const char *node_flag) {
  int counter;
  int i;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(node_global2local);
  HECMW_assert(node_flag);
  HECMW_assert(global_mesh->n_node > 0);

  for (counter = 0, i = 0; i < global_mesh->n_node; i++) {
    if (EVAL_BIT(node_flag[i], INTERNAL) || EVAL_BIT(node_flag[i], BOUNDARY)) {
      node_global2local[i] = ++counter;
    }
  }
  local_mesh->n_node       = counter;
  local_mesh->n_node_gross = counter;

  HECMW_assert(local_mesh->n_node > 0);

  return RTC_NORMAL;
}

static int const_nn_internal(const struct hecmwST_local_mesh *global_mesh,
                             struct hecmwST_local_mesh *local_mesh,
                             const char *node_flag) {
  int counter;
  int i;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(node_flag);
  HECMW_assert(global_mesh->n_node > 0);

  for (counter = 0, i = 0; i < global_mesh->n_node; i++) {
    if (EVAL_BIT(node_flag[i], INTERNAL)) counter++;
  }
  local_mesh->nn_internal = counter;

  return 0;
}

static int const_node_internal_list(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh, int *node_global2local,
    const char *node_flag) {
  int counter;
  int i;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(node_global2local);
  HECMW_assert(node_flag);
  HECMW_assert(global_mesh->n_node > 0);

  if (local_mesh->nn_internal == 0) {
    local_mesh->node_internal_list = NULL;
    return RTC_NORMAL;
  }

  local_mesh->node_internal_list =
      (int *)HECMW_malloc(sizeof(int) * local_mesh->nn_internal);
  if (local_mesh->node_internal_list == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (counter = 0, i = 0; i < global_mesh->n_node; i++) {
    if (EVAL_BIT(node_flag[i], INTERNAL)) {
      local_mesh->node_internal_list[counter++] = node_global2local[i];
    }
  }
  HECMW_assert(counter == local_mesh->nn_internal);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*K. Inagaki */
static int set_node_global2local(const struct hecmwST_local_mesh *global_mesh,
                                 struct hecmwST_local_mesh *local_mesh,
                                 int *node_global2local, const char *node_flag,
                                 int current_domain) {
  int rtc;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(node_global2local);
  HECMW_assert(node_flag);

  switch (global_mesh->hecmw_flag_parttype) {
    case HECMW_FLAG_PARTTYPE_NODEBASED:

      rtc = set_node_global2local_internal(global_mesh, local_mesh,
                                           node_global2local, node_flag,
                                           current_domain);
      if (rtc != RTC_NORMAL) goto error;

      if (is_spdup_available(global_mesh)) {
        rtc = set_node_global2local_external_mod(global_mesh, local_mesh,
                                                 node_global2local, node_flag,
                                                 current_domain);
      } else {
        rtc = set_node_global2local_external(global_mesh, local_mesh,
                                             node_global2local, node_flag);
      }

      if (rtc != RTC_NORMAL) goto error;

      local_mesh->node_internal_list = NULL;

      break;

    case HECMW_FLAG_PARTTYPE_ELEMBASED:

      rtc = const_nn_internal(global_mesh, local_mesh, node_flag);
      if (rtc != RTC_NORMAL) goto error;

      rtc = set_node_global2local_all(global_mesh, local_mesh,
                                      node_global2local, node_flag);
      if (rtc != RTC_NORMAL) goto error;

      rtc = const_node_internal_list(global_mesh, local_mesh, node_global2local,
                                     node_flag);
      if (rtc != RTC_NORMAL) goto error;

      break;

    default:
      HECMW_set_error(HECMW_PART_E_INVALID_PTYPE, "%d",
                      global_mesh->hecmw_flag_parttype);
      goto error;
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*K. Inagaki */
static int clear_node_global2local(const struct hecmwST_local_mesh *global_mesh,
                                   struct hecmwST_local_mesh *local_mesh,
                                   int *node_global2local, int domain) {
  int rtc;
  int i, node;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(node_global2local);

  if (is_spdup_available(global_mesh)) {
    for (i = 0; i < n_int_nlist[domain]; i++) {
      node                        = int_nlist[domain][i];
      node_global2local[node - 1] = 0;
    }
    for (i = n_bnd_nlist[2 * domain]; i < n_bnd_nlist[2 * domain + 1]; i++) {
      node                        = bnd_nlist[domain][i];
      node_global2local[node - 1] = 0;
    }
  } else {
    for (i = 0; i < global_mesh->n_node; i++) {
      node_global2local[i] = 0;
    }
  }

  return RTC_NORMAL;
}

/*------------------------------------------------------------------------------------------------*/

static int set_node_local2global(const struct hecmwST_local_mesh *global_mesh,
                                 struct hecmwST_local_mesh *local_mesh,
                                 const int *node_global2local,
                                 int *node_local2global) {
  int counter;
  int i;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(node_global2local);
  HECMW_assert(node_local2global);
  HECMW_assert(global_mesh->n_node > 0);

  for (counter = 0, i = 0; i < global_mesh->n_node; i++) {
    if (node_global2local[i]) {
      node_local2global[node_global2local[i] - 1] = i + 1;
      counter++;
    }
  }
  HECMW_assert(counter == local_mesh->n_node);

  return RTC_NORMAL;
}

/*K. Inagaki */
static int set_node_local2global_mod(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh, const int *node_global2local,
    int *node_local2global, int domain) {
  int counter;
  int i, idx1, idx2, node1, node2, n_int, n_bnd, n_out, maxn;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(node_global2local);
  HECMW_assert(node_local2global);
  HECMW_assert(global_mesh->n_node > 0);

  n_int = n_int_nlist[domain];
  n_bnd = n_bnd_nlist[2 * domain];
  n_out = n_bnd_nlist[2 * domain + 1] - n_bnd_nlist[2 * domain];
  maxn  = global_mesh->n_node + 1;

  node1 = (n_int == 0) ? maxn : int_nlist[domain][0];
  node2 = (n_out == 0) ? maxn : bnd_nlist[domain][n_bnd];
  for (counter = 0, idx1 = 0, idx2 = 0, i = 0; i < n_int + n_out; i++) {
    if (node1 < node2) {
      node_local2global[node_global2local[node1 - 1] - 1] = node1;
      idx1++;
      node1 = (idx1 == n_int) ? maxn : int_nlist[domain][idx1];
    } else {
      node_local2global[node_global2local[node2 - 1] - 1] = node2;
      idx2++;
      node2 = (idx2 == n_out) ? maxn : bnd_nlist[domain][idx2 + n_bnd];
    }
    counter++;
  }

  HECMW_assert(counter == local_mesh->n_node);

  return RTC_NORMAL;
}

/*------------------------------------------------------------------------------------------------*/

static int set_elem_global2local_internal(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh, int *elem_global2local,
    const char *elem_flag) {
  int counter;
  int i;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(elem_global2local);
  HECMW_assert(elem_flag);
  HECMW_assert(global_mesh->n_elem);

  for (counter = 0, i = 0; i < global_mesh->n_elem; i++) {
    if (EVAL_BIT(elem_flag[i], INTERNAL)) {
      elem_global2local[i] = ++counter;
    }
  }
  local_mesh->ne_internal = counter;

  return RTC_NORMAL;
}

static int set_elem_global2local_external(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh, int *elem_global2local,
    const char *elem_flag) {
  int counter;
  int i;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(elem_global2local);
  HECMW_assert(elem_flag);
  HECMW_assert(global_mesh->n_elem);

  for (counter = local_mesh->ne_internal, i = 0; i < global_mesh->n_elem; i++) {
    if (!EVAL_BIT(elem_flag[i], INTERNAL) && EVAL_BIT(elem_flag[i], BOUNDARY)) {
      elem_global2local[i] = ++counter;
    }
  }
  local_mesh->n_elem       = counter;
  local_mesh->n_elem_gross = counter;

  HECMW_assert(local_mesh->n_elem > 0);

  return RTC_NORMAL;
}

static int set_elem_global2local_all(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh, int *elem_global2local,
    const char *elem_flag) {
  int counter;
  int i;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(elem_global2local);
  HECMW_assert(elem_flag);
  HECMW_assert(global_mesh->n_elem > 0);

  for (counter = 0, i = 0; i < global_mesh->n_elem; i++) {
    if (EVAL_BIT(elem_flag[i], INTERNAL) || EVAL_BIT(elem_flag[i], BOUNDARY)) {
      elem_global2local[i] = ++counter;
    }
  }
  local_mesh->n_elem       = counter;
  local_mesh->n_elem_gross = counter;

  HECMW_assert(local_mesh->n_elem > 0);

  return RTC_NORMAL;
}

/*K. Inagaki */
static int set_elem_global2local_all_mod(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh, int *elem_global2local,
    const char *elem_flag, int domain) {
  int counter;
  int i, idx1, idx2, elem1, elem2, n_int, n_bnd, n_out, maxe;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(elem_global2local);
  HECMW_assert(elem_flag);
  HECMW_assert(global_mesh->n_elem > 0);

  n_int = n_int_elist[domain];
  n_bnd = n_bnd_elist[2 * domain];
  n_out = n_bnd_elist[2 * domain + 1] - n_bnd_elist[2 * domain];
  maxe  = global_mesh->n_elem + 1;

  elem1 = (n_int == 0) ? maxe : int_elist[domain][0];
  elem2 = (n_out == 0) ? maxe : bnd_elist[domain][n_bnd];
  for (counter = 0, idx1 = 0, idx2 = 0, i = 0; i < n_int + n_out; i++) {
    if (elem1 < elem2) {
      elem_global2local[elem1 - 1] = ++counter;
      idx1++;
      elem1 = (idx1 == n_int) ? maxe : int_elist[domain][idx1];
    } else {
      elem_global2local[elem2 - 1] = ++counter;
      idx2++;
      elem2 = (idx2 == n_out) ? maxe : bnd_elist[domain][idx2 + n_bnd];
    }
  }

  local_mesh->n_elem       = counter;
  local_mesh->n_elem_gross = counter;

  HECMW_assert(local_mesh->n_elem > 0);

  return RTC_NORMAL;
}

static int const_ne_internal(const struct hecmwST_local_mesh *global_mesh,
                             struct hecmwST_local_mesh *local_mesh,
                             const char *elem_flag) {
  int counter;
  int i;

  HECMW_assert(global_mesh->n_elem > 0);

  for (counter = 0, i = 0; i < global_mesh->n_elem; i++) {
    if (EVAL_BIT(elem_flag[i], INTERNAL)) counter++;
  }
  local_mesh->ne_internal = counter;

  return RTC_NORMAL;
}

/*K. Inagaki */
static int const_elem_internal_list(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh, int *elem_global2local,
    const char *elem_flag, int domain) {
  int counter;
  int i, elem;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(elem_global2local);
  HECMW_assert(elem_flag);
  HECMW_assert(global_mesh->n_elem > 0);

  if (local_mesh->ne_internal == 0) {
    local_mesh->elem_internal_list = NULL;
    return RTC_NORMAL;
  }

  local_mesh->elem_internal_list =
      (int *)HECMW_malloc(sizeof(int) * local_mesh->ne_internal);
  if (local_mesh->elem_internal_list == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (counter = 0, i = 0; i < n_int_elist[domain]; i++) {
    elem                                      = int_elist[domain][i];
    local_mesh->elem_internal_list[counter++] = elem_global2local[elem - 1];
  }

  HECMW_assert(counter == local_mesh->ne_internal);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int set_elem_global2local(const struct hecmwST_local_mesh *global_mesh,
                                 struct hecmwST_local_mesh *local_mesh,
                                 int *elem_global2local, const char *elem_flag,
                                 int current_domain) {
  int rtc;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(elem_global2local);
  HECMW_assert(elem_flag);

  switch (global_mesh->hecmw_flag_parttype) {
    case HECMW_FLAG_PARTTYPE_NODEBASED: /* for node-based partitioning */

      local_mesh->ne_internal = n_int_elist[current_domain];

      if (is_spdup_available(global_mesh)) {
        rtc = set_elem_global2local_all_mod(global_mesh, local_mesh,
                                            elem_global2local, elem_flag,
                                            current_domain);
      } else {
        rtc = set_elem_global2local_all(global_mesh, local_mesh,
                                        elem_global2local, elem_flag);
      }

      if (rtc != RTC_NORMAL) goto error;

      rtc = const_elem_internal_list(global_mesh, local_mesh, elem_global2local,
                                     elem_flag, current_domain);

      if (rtc != RTC_NORMAL) goto error;

      break;

    case HECMW_FLAG_PARTTYPE_ELEMBASED: /* for element-based partitioning */

      rtc = set_elem_global2local_internal(global_mesh, local_mesh,
                                           elem_global2local, elem_flag);
      if (rtc != RTC_NORMAL) goto error;

      rtc = set_elem_global2local_external(global_mesh, local_mesh,
                                           elem_global2local, elem_flag);
      if (rtc != RTC_NORMAL) goto error;

      local_mesh->elem_internal_list = NULL;

      break;

    default:
      HECMW_set_error(HECMW_PART_E_INVALID_PTYPE, "%d",
                      global_mesh->hecmw_flag_parttype);
      goto error;
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int clear_elem_global2local(const struct hecmwST_local_mesh *global_mesh,
                                   struct hecmwST_local_mesh *local_mesh,
                                   int *elem_global2local, int domain) {
  int rtc;
  int i, elem;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(elem_global2local);

  if (is_spdup_available(global_mesh)) {
    for (i = 0; i < n_int_elist[domain]; i++) {
      elem                        = int_elist[domain][i];
      elem_global2local[elem - 1] = 0;
    }
    for (i = n_bnd_elist[2 * domain]; i < n_bnd_elist[2 * domain + 1]; i++) {
      elem                        = bnd_elist[domain][i];
      elem_global2local[elem - 1] = 0;
    }

  } else {
    for (i = 0; i < global_mesh->n_elem; i++) {
      elem_global2local[i] = 0;
    }
  }

  return RTC_NORMAL;
}

/*------------------------------------------------------------------------------------------------*/

static int set_elem_local2global(const struct hecmwST_local_mesh *global_mesh,
                                 struct hecmwST_local_mesh *local_mesh,
                                 const int *elem_global2local,
                                 int *elem_local2global) {
  int counter;
  int i;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(elem_global2local);
  HECMW_assert(elem_local2global);
  HECMW_assert(global_mesh->n_elem > 0);

  for (counter = 0, i = 0; i < global_mesh->n_elem; i++) {
    if (elem_global2local[i]) {
      elem_local2global[elem_global2local[i] - 1] = i + 1;
      counter++;
    }
  }
  HECMW_assert(counter == local_mesh->n_elem);

  return RTC_NORMAL;
}

/*K. Inagaki */
static int set_elem_local2global_mod(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh, const int *elem_global2local,
    int *elem_local2global, int domain) {
  int counter;
  int i, idx1, idx2, elem1, elem2, n_int, n_bnd, n_out, maxe;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(elem_global2local);
  HECMW_assert(elem_local2global);
  HECMW_assert(global_mesh->n_elem > 0);

  n_int = n_int_elist[domain];
  n_bnd = n_bnd_elist[2 * domain];
  n_out = n_bnd_elist[2 * domain + 1] - n_bnd_elist[2 * domain];
  maxe  = global_mesh->n_elem + 1;

  elem1 = (n_int == 0) ? maxe : int_elist[domain][0];
  elem2 = (n_out == 0) ? maxe : bnd_elist[domain][n_bnd];
  for (counter = 0, idx1 = 0, idx2 = 0, i = 0; i < n_int + n_out; i++) {
    if (elem1 < elem2) {
      elem_local2global[elem_global2local[elem1 - 1] - 1] = elem1;
      idx1++;
      elem1 = (idx1 == n_int) ? maxe : int_elist[domain][idx1];
    } else {
      elem_local2global[elem_global2local[elem2 - 1] - 1] = elem2;
      idx2++;
      elem2 = (idx2 == n_out) ? maxe : bnd_elist[domain][idx2 + n_bnd];
    }
    counter++;
  }

  HECMW_assert(counter == local_mesh->n_elem);

  return RTC_NORMAL;
}

/*================================================================================================*/

static int const_gridfile(const struct hecmwST_local_mesh *global_mesh,
                          struct hecmwST_local_mesh *local_mesh) {
  strcpy(local_mesh->gridfile, global_mesh->gridfile);

  return RTC_NORMAL;
}

static int const_hecmw_n_file(const struct hecmwST_local_mesh *global_mesh,
                              struct hecmwST_local_mesh *local_mesh) {
  local_mesh->hecmw_n_file = global_mesh->hecmw_n_file;

  return RTC_NORMAL;
}

static int const_files(const struct hecmwST_local_mesh *global_mesh,
                       struct hecmwST_local_mesh *local_mesh) {
  local_mesh->files = global_mesh->files;

  return RTC_NORMAL;
}

static int const_header(const struct hecmwST_local_mesh *global_mesh,
                        struct hecmwST_local_mesh *local_mesh) {
  strcpy(local_mesh->header, global_mesh->header);

  return RTC_NORMAL;
}

static int const_hecmw_flag_adapt(const struct hecmwST_local_mesh *global_mesh,
                                  struct hecmwST_local_mesh *local_mesh) {
  local_mesh->hecmw_flag_adapt = global_mesh->hecmw_flag_adapt;

  return RTC_NORMAL;
}

static int const_hecmw_flag_initcon(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh) {
  local_mesh->hecmw_flag_initcon = global_mesh->hecmw_flag_initcon;

  return RTC_NORMAL;
}

static int const_hecmw_flag_parttype(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh) {
  local_mesh->hecmw_flag_parttype = global_mesh->hecmw_flag_parttype;

  return RTC_NORMAL;
}

static int const_hecmw_flag_partdepth(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh) {
  local_mesh->hecmw_flag_partdepth = global_mesh->hecmw_flag_partdepth;

  return RTC_NORMAL;
}

static int const_hecmw_flag_version(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh) {
  local_mesh->hecmw_flag_version = global_mesh->hecmw_flag_version;

  return RTC_NORMAL;
}

static int const_hecmw_flag_partcontact(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh) {
  local_mesh->hecmw_flag_partcontact = global_mesh->hecmw_flag_partcontact;

  return RTC_NORMAL;
}

static int const_zero_temp(const struct hecmwST_local_mesh *global_mesh,
                           struct hecmwST_local_mesh *local_mesh) {
  local_mesh->zero_temp = global_mesh->zero_temp;

  return RTC_NORMAL;
}

static int const_global_info(const struct hecmwST_local_mesh *global_mesh,
                             struct hecmwST_local_mesh *local_mesh) {
  int rtc;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);

  rtc = const_gridfile(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_hecmw_n_file(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_files(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_header(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_hecmw_flag_adapt(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_hecmw_flag_initcon(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_hecmw_flag_parttype(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_hecmw_flag_partdepth(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_hecmw_flag_version(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_hecmw_flag_partcontact(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_zero_temp(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * - - - - - - - - - */

static int const_n_dof(const struct hecmwST_local_mesh *global_mesh,
                       struct hecmwST_local_mesh *local_mesh) {
  HECMW_assert(global_mesh->n_dof > 0);

  local_mesh->n_dof = global_mesh->n_dof;

  HECMW_assert(local_mesh->n_dof > 0);

  return RTC_NORMAL;
}

static int const_n_dof_grp(const struct hecmwST_local_mesh *global_mesh,
                           struct hecmwST_local_mesh *local_mesh) {
  HECMW_assert(global_mesh->n_dof_grp);

  local_mesh->n_dof_grp = global_mesh->n_dof_grp;

  HECMW_assert(global_mesh->n_dof_grp);

  return RTC_NORMAL;
}

static int const_node_dof_index(const struct hecmwST_local_mesh *global_mesh,
                                struct hecmwST_local_mesh *local_mesh,
                                const char *node_flag) {
  int counter;
  int i, j;

  HECMW_assert(local_mesh->n_dof_grp > 0);
  HECMW_assert(global_mesh->node_dof_index);

  local_mesh->node_dof_index =
      (int *)HECMW_calloc(local_mesh->n_dof_grp + 1, sizeof(int));
  if (local_mesh->node_dof_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (counter = 0, i = 0; i < global_mesh->n_dof_grp; i++) {
    for (j = global_mesh->node_dof_index[i];
         j < global_mesh->node_dof_index[i + 1]; j++) {
      if (EVAL_BIT(node_flag[j], INTERNAL)) counter++;
    }
    local_mesh->node_dof_index[i + 1] = counter;
  }
  HECMW_assert(local_mesh->node_dof_index[local_mesh->n_dof_grp] ==
               local_mesh->nn_internal);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*K. Inagaki */
static int const_node_dof_index_mod(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh, const char *node_flag, int domain) {
  int counter;
  int i, j, node;

  HECMW_assert(local_mesh->n_dof_grp > 0);
  HECMW_assert(global_mesh->node_dof_index);

  local_mesh->node_dof_index =
      (int *)HECMW_calloc(local_mesh->n_dof_grp + 1, sizeof(int));
  if (local_mesh->node_dof_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (counter = 0, i = 0; i < global_mesh->n_dof_grp; i++) {
    for (j = 0; j < n_int_nlist[domain]; j++) {
      node = int_nlist[domain][j];
      if (node <= global_mesh->node_dof_index[i]) continue;
      if (node > global_mesh->node_dof_index[i + 1]) continue;
      counter++;
    }
    local_mesh->node_dof_index[i + 1] = counter;
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_node_dof_item(const struct hecmwST_local_mesh *global_mesh,
                               struct hecmwST_local_mesh *local_mesh) {
  HECMW_assert(global_mesh->node_dof_item);

  local_mesh->node_dof_item = global_mesh->node_dof_item;

  return 0;
}

static int const_node(const struct hecmwST_local_mesh *global_mesh,
                      struct hecmwST_local_mesh *local_mesh,
                      const int *node_local2global) {
  int i;

  HECMW_assert(local_mesh->n_node > 0);
  HECMW_assert(global_mesh->node);

  local_mesh->node =
      (double *)HECMW_malloc(sizeof(double) * local_mesh->n_node * 3);
  if (local_mesh->node == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (i = 0; i < local_mesh->n_node; i++) {
    local_mesh->node[3 * i] = global_mesh->node[3 * (node_local2global[i] - 1)];
    local_mesh->node[3 * i + 1] =
        global_mesh->node[3 * (node_local2global[i] - 1) + 1];
    local_mesh->node[3 * i + 2] =
        global_mesh->node[3 * (node_local2global[i] - 1) + 2];
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_node_id(const struct hecmwST_local_mesh *global_mesh,
                         struct hecmwST_local_mesh *local_mesh,
                         const int *node_local2global) {
  int i;

  HECMW_assert(local_mesh->n_node > 0);
  HECMW_assert(global_mesh->node_ID);

  local_mesh->node_ID =
      (int *)HECMW_malloc(sizeof(int) * local_mesh->n_node * 2);
  if (local_mesh->node_ID == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (i = 0; i < local_mesh->n_node; i++) {
    local_mesh->node_ID[2 * i] =
        global_mesh->node_ID[2 * (node_local2global[i] - 1)];
    local_mesh->node_ID[2 * i + 1] =
        global_mesh->node_ID[2 * (node_local2global[i] - 1) + 1];
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_global_node_id(const struct hecmwST_local_mesh *global_mesh,
                                struct hecmwST_local_mesh *local_mesh,
                                const int *node_local2global) {
  int i;

  HECMW_assert(local_mesh->n_node > 0);
  HECMW_assert(global_mesh->global_node_ID);

  local_mesh->global_node_ID =
      (int *)HECMW_malloc(sizeof(int) * local_mesh->n_node);
  if (local_mesh->global_node_ID == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (i = 0; i < local_mesh->n_node; i++) {
    local_mesh->global_node_ID[i] =
        global_mesh->global_node_ID[node_local2global[i] - 1];
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_node_init_val_index(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh, const int *node_local2global) {
  int old_idx;
  int i;

  HECMW_assert(local_mesh->hecmw_flag_initcon);
  HECMW_assert(local_mesh->n_node > 0);
  HECMW_assert(global_mesh->node_init_val_index);

  local_mesh->node_init_val_index =
      (int *)HECMW_calloc(local_mesh->n_node + 1, sizeof(int));
  if (local_mesh->node_init_val_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (i = 0; i < local_mesh->n_node; i++) {
    old_idx = node_local2global[i] - 1;

    local_mesh->node_init_val_index[i + 1] =
        local_mesh->node_init_val_index[i] +
        global_mesh->node_init_val_index[old_idx + 1] -
        global_mesh->node_init_val_index[old_idx];
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_node_init_val_item(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh, const int *node_local2global) {
  int size;
  int counter;
  int i, j, gstart, gend, lstart, lend;

  HECMW_assert(local_mesh->hecmw_flag_initcon);
  HECMW_assert(local_mesh->n_node > 0);
  HECMW_assert(local_mesh->node_init_val_index);
  HECMW_assert(global_mesh->node_init_val_item);

  if (local_mesh->node_init_val_index[local_mesh->n_node] == 0) {
    local_mesh->node_init_val_item = NULL;
    return 0;
  }

  size = sizeof(double) * local_mesh->node_init_val_index[local_mesh->n_node];
  local_mesh->node_init_val_item = (double *)HECMW_malloc(size);
  if (local_mesh->node_init_val_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (counter = 0, i = 0; i < local_mesh->n_node; i++) {
    gstart = global_mesh->node_init_val_index[node_local2global[i] - 1];
    gend   = global_mesh->node_init_val_index[node_local2global[i]];
    lstart = local_mesh->node_init_val_index[i];
    lend   = local_mesh->node_init_val_index[i + 1];

    HECMW_assert(gend - gstart == lend - lstart);

    for (j = 0; j < lend - lstart; j++) {
      local_mesh->node_init_val_item[lstart + j] =
          global_mesh->node_init_val_item[gstart + j];
      counter++;
    }
    HECMW_assert(counter == local_mesh->node_init_val_index[i + 1]);
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_node_info(const struct hecmwST_local_mesh *global_mesh,
                           struct hecmwST_local_mesh *local_mesh,
                           const int *node_local2global, const char *node_flag,
                           int current_domain) {
  int rtc;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(node_local2global);
  HECMW_assert(node_flag);

  rtc = const_n_dof(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_n_dof_grp(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  switch (global_mesh->hecmw_flag_parttype) {
    case HECMW_FLAG_PARTTYPE_NODEBASED:
      rtc = const_node_dof_index_mod(global_mesh, local_mesh, node_flag,
                                     current_domain);
      break;
    case HECMW_FLAG_PARTTYPE_ELEMBASED:
      rtc = const_node_dof_index(global_mesh, local_mesh, node_flag);
      break;
    default:
      HECMW_set_error(errno, "");
      goto error;
  }

  if (rtc != RTC_NORMAL) goto error;

  rtc = const_node_dof_item(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_node(global_mesh, local_mesh, node_local2global);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_node_id(global_mesh, local_mesh, node_local2global);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_global_node_id(global_mesh, local_mesh, node_local2global);
  if (rtc != RTC_NORMAL) goto error;

  if (local_mesh->hecmw_flag_initcon) {
    rtc = const_node_init_val_index(global_mesh, local_mesh, node_local2global);
    if (rtc != RTC_NORMAL) goto error;

    rtc = const_node_init_val_item(global_mesh, local_mesh, node_local2global);
    if (rtc != RTC_NORMAL) goto error;
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * - - - - - - - - - */

static int const_n_elem_type(const struct hecmwST_local_mesh *global_mesh,
                             struct hecmwST_local_mesh *local_mesh) {
  HECMW_assert(global_mesh->n_elem_type > 0);

  local_mesh->n_elem_type = global_mesh->n_elem_type;

  HECMW_assert(local_mesh->n_elem_type > 0);

  return RTC_NORMAL;
}

static int const_elem_type(const struct hecmwST_local_mesh *global_mesh,
                           struct hecmwST_local_mesh *local_mesh,
                           const int *elem_local2global) {
  int i;

  HECMW_assert(local_mesh->n_elem > 0);
  HECMW_assert(global_mesh->elem_type);

  local_mesh->elem_type = (int *)HECMW_malloc(sizeof(int) * local_mesh->n_elem);
  if (local_mesh->elem_type == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (i = 0; i < local_mesh->n_elem; i++) {
    local_mesh->elem_type[i] = global_mesh->elem_type[elem_local2global[i] - 1];
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_elem_type_index(const struct hecmwST_local_mesh *global_mesh,
                                 struct hecmwST_local_mesh *local_mesh,
                                 const int *elem_global2local) {
  int counter;
  int i, j;

  HECMW_assert(local_mesh->n_elem_type > 0);
  HECMW_assert(global_mesh->n_elem_type > 0);
  HECMW_assert(global_mesh->elem_type_index);

  local_mesh->elem_type_index =
      (int *)HECMW_calloc(local_mesh->n_elem_type + 1, sizeof(int));
  if (local_mesh->elem_type_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (counter = 0, i = 0; i < global_mesh->n_elem_type; i++) {
    for (j = global_mesh->elem_type_index[i];
         j < global_mesh->elem_type_index[i + 1]; j++) {
      if (elem_global2local[j]) counter++;
    }
    local_mesh->elem_type_index[i + 1] = counter;
  }
  HECMW_assert(local_mesh->elem_type_index[local_mesh->n_elem_type] ==
               local_mesh->n_elem);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*K. Inagaki */
static int const_elem_type_index_mod(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh, const int *elem_global2local,
    int domain) {
  int counter;
  int i, j, idx1, idx2, elem_tmp, elem1, elem2, n_int, n_bnd, n_out, maxe;

  HECMW_assert(local_mesh->n_elem_type > 0);
  HECMW_assert(global_mesh->n_elem_type > 0);
  HECMW_assert(global_mesh->elem_type_index);

  local_mesh->elem_type_index =
      (int *)HECMW_calloc(local_mesh->n_elem_type + 1, sizeof(int));
  if (local_mesh->elem_type_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  n_int = n_int_elist[domain];
  n_bnd = n_bnd_elist[2 * domain];
  n_out = n_bnd_elist[2 * domain + 1] - n_bnd_elist[2 * domain];
  maxe  = global_mesh->n_elem + 1;

  for (counter = 0, i = 0; i < global_mesh->n_elem_type; i++) {
    elem1 = (n_int == 0) ? maxe : int_elist[domain][0];
    elem2 = (n_out == 0) ? maxe : bnd_elist[domain][n_bnd];
    for (idx1 = 0, idx2 = 0, j = 0; j < n_int + n_out; j++) {
      if (elem1 < elem2) {
        elem_tmp = elem1 - 1;
        idx1++;
        elem1 = (idx1 == n_int) ? maxe : int_elist[domain][idx1];
      } else {
        elem_tmp = elem2 - 1;
        idx2++;
        elem2 = (idx2 == n_out) ? maxe : bnd_elist[domain][idx2 + n_bnd];
      }
      if (elem_tmp >= global_mesh->elem_type_index[i] &&
          elem_tmp < global_mesh->elem_type_index[i + 1]) {
        counter++;
      }
    }
    local_mesh->elem_type_index[i + 1] = counter;
  }

  HECMW_assert(local_mesh->elem_type_index[local_mesh->n_elem_type] ==
               local_mesh->n_elem);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_elem_type_item(const struct hecmwST_local_mesh *global_mesh,
                                struct hecmwST_local_mesh *local_mesh) {
  HECMW_assert(global_mesh->elem_type_item);

  local_mesh->elem_type_item = global_mesh->elem_type_item;

  return RTC_NORMAL;
}

static int const_elem_node_index(const struct hecmwST_local_mesh *global_mesh,
                                 struct hecmwST_local_mesh *local_mesh,
                                 const int *elem_local2global) {
  int old_idx;
  int i;

  HECMW_assert(local_mesh->n_elem > 0);
  HECMW_assert(global_mesh->elem_node_index);

  local_mesh->elem_node_index =
      (int *)HECMW_calloc(local_mesh->n_elem + 1, sizeof(int));
  if (local_mesh->elem_node_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (i = 0; i < local_mesh->n_elem; i++) {
    old_idx = elem_local2global[i] - 1;

    local_mesh->elem_node_index[i + 1] =
        local_mesh->elem_node_index[i] +
        global_mesh->elem_node_index[old_idx + 1] -
        global_mesh->elem_node_index[old_idx];
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_elem_node_item(const struct hecmwST_local_mesh *global_mesh,
                                struct hecmwST_local_mesh *local_mesh,
                                const int *node_global2local,
                                const int *elem_local2global) {
  int node;
  int size;
  int counter;
  int i, j, gstart, gend, lstart, lend;

  HECMW_assert(local_mesh->n_elem > 0);
  HECMW_assert(local_mesh->elem_node_index);
  HECMW_assert(local_mesh->elem_node_index[local_mesh->n_elem] > 0);
  HECMW_assert(global_mesh->elem_node_item);

  size = sizeof(int) * local_mesh->elem_node_index[local_mesh->n_elem];
  local_mesh->elem_node_item = (int *)HECMW_malloc(size);
  if (local_mesh->elem_node_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (counter = 0, i = 0; i < local_mesh->n_elem; i++) {
    gstart = global_mesh->elem_node_index[elem_local2global[i] - 1];
    gend   = global_mesh->elem_node_index[elem_local2global[i]];
    lstart = local_mesh->elem_node_index[i];
    lend   = local_mesh->elem_node_index[i + 1];

    for (j = 0; j < lend - lstart; j++) {
      node = global_mesh->elem_node_item[gstart + j];
      local_mesh->elem_node_item[lstart + j] = node_global2local[node - 1];
      counter++;
    }
    HECMW_assert(counter == local_mesh->elem_node_index[i + 1]);
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_elem_id(const struct hecmwST_local_mesh *global_mesh,
                         struct hecmwST_local_mesh *local_mesh,
                         const int *elem_local2global) {
  int i;

  HECMW_assert(local_mesh->n_elem > 0);
  HECMW_assert(global_mesh->elem_ID);

  local_mesh->elem_ID =
      (int *)HECMW_malloc(sizeof(int) * local_mesh->n_elem * 2);
  if (local_mesh->elem_ID == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (i = 0; i < local_mesh->n_elem; i++) {
    local_mesh->elem_ID[2 * i] =
        global_mesh->elem_ID[2 * (elem_local2global[i] - 1)];
    local_mesh->elem_ID[2 * i + 1] =
        global_mesh->elem_ID[2 * (elem_local2global[i] - 1) + 1];
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_global_elem_id(const struct hecmwST_local_mesh *global_mesh,
                                struct hecmwST_local_mesh *local_mesh,
                                const int *elem_local2global) {
  int i;

  HECMW_assert(local_mesh->n_elem);
  HECMW_assert(global_mesh->global_elem_ID);

  local_mesh->global_elem_ID =
      (int *)HECMW_malloc(sizeof(int) * local_mesh->n_elem);
  if (local_mesh->global_elem_ID == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (i = 0; i < local_mesh->n_elem; i++) {
    local_mesh->global_elem_ID[i] =
        global_mesh->global_elem_ID[elem_local2global[i] - 1];
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_section_id(const struct hecmwST_local_mesh *global_mesh,
                            struct hecmwST_local_mesh *local_mesh,
                            const int *elem_local2global) {
  int i;

  HECMW_assert(local_mesh->n_elem);
  HECMW_assert(global_mesh->section_ID);

  local_mesh->section_ID =
      (int *)HECMW_malloc(sizeof(int) * local_mesh->n_elem);
  if (local_mesh->section_ID == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (i = 0; i < local_mesh->n_elem; i++) {
    local_mesh->section_ID[i] =
        global_mesh->section_ID[elem_local2global[i] - 1];
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_elem_mat_id_index(const struct hecmwST_local_mesh *global_mesh,
                                   struct hecmwST_local_mesh *local_mesh,
                                   const int *elem_local2global) {
  int old_idx;
  int i;

  HECMW_assert(local_mesh->n_elem > 0);
  HECMW_assert(global_mesh->elem_mat_ID_index);

  local_mesh->elem_mat_ID_index =
      (int *)HECMW_calloc(local_mesh->n_elem + 1, sizeof(int));
  if (local_mesh->elem_mat_ID_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (i = 0; i < local_mesh->n_elem; i++) {
    old_idx = elem_local2global[i] - 1;

    local_mesh->elem_mat_ID_index[i + 1] =
        local_mesh->elem_mat_ID_index[i] +
        global_mesh->elem_mat_ID_index[old_idx + 1] -
        global_mesh->elem_mat_ID_index[old_idx];
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_n_elem_mat_id(struct hecmwST_local_mesh *local_mesh) {
  HECMW_assert(local_mesh->n_elem > 0);
  HECMW_assert(local_mesh->elem_mat_ID_index);

  local_mesh->n_elem_mat_ID = local_mesh->elem_mat_ID_index[local_mesh->n_elem];

  return RTC_NORMAL;
}

static int const_elem_mat_id_item(const struct hecmwST_local_mesh *global_mesh,
                                  struct hecmwST_local_mesh *local_mesh,
                                  const int *elem_local2global) {
  int size;
  int counter;
  int i, j, gstart, gend, lstart, lend;

  HECMW_assert(local_mesh->n_elem > 0);
  HECMW_assert(local_mesh->elem_mat_ID_index[local_mesh->n_elem] >= 0);

  if (local_mesh->elem_mat_ID_index[local_mesh->n_elem] == 0) {
    local_mesh->elem_mat_ID_item = NULL;
    return RTC_NORMAL;
  }

  size = sizeof(int) * local_mesh->elem_mat_ID_index[local_mesh->n_elem];
  local_mesh->elem_mat_ID_item = (int *)HECMW_malloc(size);
  if (local_mesh->elem_mat_ID_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (counter = 0, i = 0; i < local_mesh->n_elem; i++) {
    gstart = global_mesh->elem_mat_ID_index[elem_local2global[i] - 1];
    gend   = global_mesh->elem_mat_ID_index[elem_local2global[i]];
    lstart = local_mesh->elem_mat_ID_index[i];
    lend   = local_mesh->elem_mat_ID_index[i + 1];

    HECMW_assert(lend - lstart == gend - gstart);

    for (j = 0; j < lend - lstart; j++) {
      local_mesh->elem_mat_ID_item[lstart + j] =
          global_mesh->elem_mat_ID_item[gstart + j];
      counter++;
    }
    HECMW_assert(counter == local_mesh->elem_mat_ID_index[i + 1]);
  }
  HECMW_assert(counter == local_mesh->n_elem_mat_ID);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_elem_info(const struct hecmwST_local_mesh *global_mesh,
                           struct hecmwST_local_mesh *local_mesh,
                           const int *node_global2local,
                           const int *elem_global2local,
                           const int *elem_local2global, int current_domain) {
  int rtc;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(node_global2local);
  HECMW_assert(elem_global2local);
  HECMW_assert(elem_local2global);

  rtc = const_n_elem_type(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_elem_type(global_mesh, local_mesh, elem_local2global);
  if (rtc != RTC_NORMAL) goto error;

  if (is_spdup_available(global_mesh)) {
    rtc = const_elem_type_index_mod(global_mesh, local_mesh, elem_global2local,
                                    current_domain);
  } else {
    rtc = const_elem_type_index(global_mesh, local_mesh, elem_global2local);
  }

  if (rtc != RTC_NORMAL) goto error;

  rtc = const_elem_type_item(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_elem_node_index(global_mesh, local_mesh, elem_local2global);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_elem_node_item(global_mesh, local_mesh, node_global2local,
                             elem_local2global);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_elem_id(global_mesh, local_mesh, elem_local2global);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_global_elem_id(global_mesh, local_mesh, elem_local2global);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_section_id(global_mesh, local_mesh, elem_local2global);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_elem_mat_id_index(global_mesh, local_mesh, elem_local2global);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_n_elem_mat_id(local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_elem_mat_id_item(global_mesh, local_mesh, elem_local2global);
  if (rtc != RTC_NORMAL) goto error;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * - - - - - - - - - */

static int const_hecmw_comm(const struct hecmwST_local_mesh *global_mesh,
                            struct hecmwST_local_mesh *local_mesh) {
  local_mesh->HECMW_COMM = global_mesh->HECMW_COMM;

  return RTC_NORMAL;
}

static int const_zero(struct hecmwST_local_mesh *local_mesh,
                      int current_domain) {
  local_mesh->zero = (current_domain == 0) ? 1 : 0;

  return RTC_NORMAL;
}

static int const_petot(const struct hecmwST_local_mesh *global_mesh,
                       struct hecmwST_local_mesh *local_mesh) {
  local_mesh->PETOT = global_mesh->n_subdomain;

  return RTC_NORMAL;
}

static int const_pesmptot(const struct hecmwST_local_mesh *global_mesh,
                          struct hecmwST_local_mesh *local_mesh) {
  local_mesh->PEsmpTOT = global_mesh->PEsmpTOT;

  return RTC_NORMAL;
}

static int const_my_rank(struct hecmwST_local_mesh *local_mesh,
                         int current_domain) {
  local_mesh->my_rank = current_domain;

  return RTC_NORMAL;
}

static int const_errnof(const struct hecmwST_local_mesh *global_mesh,
                        struct hecmwST_local_mesh *local_mesh) {
  local_mesh->errnof = global_mesh->errnof;

  return RTC_NORMAL;
}

static int const_n_subdomain(const struct hecmwST_local_mesh *global_mesh,
                             struct hecmwST_local_mesh *local_mesh) {
  local_mesh->n_subdomain = global_mesh->n_subdomain;

  return RTC_NORMAL;
}

static int const_import_item(struct hecmwST_local_mesh *local_mesh,
                             const int *global2local) {
  int new_id;
  int i;

  if (local_mesh->n_neighbor_pe == 0) {
    local_mesh->import_item = NULL;
    return RTC_NORMAL;
  }

  HECMW_assert(local_mesh->n_neighbor_pe > 0);
  HECMW_assert(local_mesh->import_index);
  HECMW_assert(local_mesh->import_index[local_mesh->n_neighbor_pe] >= 0);
  HECMW_assert(local_mesh->import_item);

  for (i = 0; i < local_mesh->import_index[local_mesh->n_neighbor_pe]; i++) {
    new_id                     = global2local[local_mesh->import_item[i] - 1];
    local_mesh->import_item[i] = new_id;
  }

  return RTC_NORMAL;
}

static int const_export_item(struct hecmwST_local_mesh *local_mesh,
                             const int *global2local) {
  int new_id;
  int i;

  if (local_mesh->n_neighbor_pe == 0) {
    local_mesh->export_item = NULL;
    return RTC_NORMAL;
  }

  HECMW_assert(local_mesh->n_neighbor_pe > 0);
  HECMW_assert(local_mesh->export_index);
  HECMW_assert(local_mesh->export_index[local_mesh->n_neighbor_pe] >= 0);
  HECMW_assert(local_mesh->export_item);

  for (i = 0; i < local_mesh->export_index[local_mesh->n_neighbor_pe]; i++) {
    new_id                     = global2local[local_mesh->export_item[i] - 1];
    local_mesh->export_item[i] = new_id;
  }

  return RTC_NORMAL;
}

static int const_shared_item(struct hecmwST_local_mesh *local_mesh,
                             const int *global2local) {
  int new_id;
  int i;

  if (local_mesh->n_neighbor_pe == 0) {
    local_mesh->shared_item = NULL;
    return RTC_NORMAL;
  }

  HECMW_assert(local_mesh->n_neighbor_pe > 0);
  HECMW_assert(local_mesh->shared_index);
  HECMW_assert(local_mesh->shared_index[local_mesh->n_neighbor_pe] >= 0);
  HECMW_assert(local_mesh->shared_item);

  for (i = 0; i < local_mesh->shared_index[local_mesh->n_neighbor_pe]; i++) {
    new_id                     = global2local[local_mesh->shared_item[i] - 1];
    local_mesh->shared_item[i] = new_id;
  }

  return RTC_NORMAL;
}

static int const_comm_info(const struct hecmwST_local_mesh *global_mesh,
                           struct hecmwST_local_mesh *local_mesh,
                           const int *node_global2local,
                           const int *elem_global2local, int current_domain) {
  int rtc;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(node_global2local);
  HECMW_assert(elem_global2local);

  rtc = const_hecmw_comm(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_zero(local_mesh, current_domain);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_petot(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_pesmptot(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_my_rank(local_mesh, current_domain);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_errnof(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_n_subdomain(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  switch (global_mesh->hecmw_flag_parttype) {
    case HECMW_FLAG_PARTTYPE_NODEBASED:
      rtc = const_import_item(local_mesh, node_global2local);
      if (rtc != RTC_NORMAL) goto error;

      rtc = const_export_item(local_mesh, node_global2local);
      if (rtc != RTC_NORMAL) goto error;

      rtc = const_shared_item(local_mesh, elem_global2local);
      if (rtc != RTC_NORMAL) goto error;

      break;

    case HECMW_FLAG_PARTTYPE_ELEMBASED:
      rtc = const_import_item(local_mesh, elem_global2local);
      if (rtc != RTC_NORMAL) goto error;

      rtc = const_export_item(local_mesh, elem_global2local);
      if (rtc != RTC_NORMAL) goto error;

      rtc = const_shared_item(local_mesh, node_global2local);
      if (rtc != RTC_NORMAL) goto error;

      break;

    default:
      HECMW_set_error(HECMW_PART_E_INVALID_PTYPE, "%d",
                      global_mesh->hecmw_flag_parttype);
      goto error;
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * - - - - - - - - - */

static int const_n_adapt(const struct hecmwST_local_mesh *global_mesh,
                         struct hecmwST_local_mesh *local_mesh) {
  local_mesh->n_adapt = global_mesh->n_adapt;

  return RTC_NORMAL;
}

static int const_coarse_grid_level(const struct hecmwST_local_mesh *global_mesh,
                                   struct hecmwST_local_mesh *local_mesh) {
  local_mesh->coarse_grid_level = global_mesh->coarse_grid_level;

  return RTC_NORMAL;
}

static int const_when_i_was_refined_node(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh) {
  local_mesh->when_i_was_refined_node = global_mesh->when_i_was_refined_node;

  return RTC_NORMAL;
}

static int const_when_i_was_refined_elem(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh) {
  local_mesh->when_i_was_refined_elem = global_mesh->when_i_was_refined_elem;

  return RTC_NORMAL;
}

static int const_adapt_parent_type(const struct hecmwST_local_mesh *global_mesh,
                                   struct hecmwST_local_mesh *local_mesh) {
  local_mesh->adapt_parent_type = global_mesh->adapt_parent_type;

  return RTC_NORMAL;
}

static int const_adapt_type(const struct hecmwST_local_mesh *global_mesh,
                            struct hecmwST_local_mesh *local_mesh) {
  local_mesh->adapt_type = global_mesh->adapt_type;

  return RTC_NORMAL;
}

static int const_adapt_level(const struct hecmwST_local_mesh *global_mesh,
                             struct hecmwST_local_mesh *local_mesh) {
  local_mesh->adapt_level = global_mesh->adapt_level;

  return RTC_NORMAL;
}

static int const_adapt_parent(const struct hecmwST_local_mesh *global_mesh,
                              struct hecmwST_local_mesh *local_mesh) {
  local_mesh->adapt_parent = global_mesh->adapt_parent;

  return RTC_NORMAL;
}

static int const_adapt_children_index(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh) {
  local_mesh->adapt_children_index = global_mesh->adapt_children_index;

  return RTC_NORMAL;
}

static int const_adapt_children_item(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh) {
  local_mesh->adapt_children_item = global_mesh->adapt_children_item;

  return RTC_NORMAL;
}

static int const_adapt_info(const struct hecmwST_local_mesh *global_mesh,
                            struct hecmwST_local_mesh *local_mesh) {
  int rtc;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);

  rtc = const_n_adapt(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_coarse_grid_level(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_when_i_was_refined_node(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_when_i_was_refined_elem(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_adapt_parent_type(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_adapt_type(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_adapt_level(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_adapt_parent(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_adapt_children_index(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_adapt_children_item(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * - - - - - - - - - */

static int const_n_sect(const struct hecmwST_local_mesh *global_mesh,
                        struct hecmwST_local_mesh *local_mesh) {
  local_mesh->section->n_sect = global_mesh->section->n_sect;

  return RTC_NORMAL;
}

static int const_sect_type(const struct hecmwST_local_mesh *global_mesh,
                           struct hecmwST_local_mesh *local_mesh) {
  local_mesh->section->sect_type = global_mesh->section->sect_type;

  return RTC_NORMAL;
}

static int const_sect_opt(const struct hecmwST_local_mesh *global_mesh,
                          struct hecmwST_local_mesh *local_mesh) {
  local_mesh->section->sect_opt = global_mesh->section->sect_opt;

  return RTC_NORMAL;
}

static int const_sect_mat_id_index(const struct hecmwST_local_mesh *global_mesh,
                                   struct hecmwST_local_mesh *local_mesh) {
  local_mesh->section->sect_mat_ID_index =
      global_mesh->section->sect_mat_ID_index;

  return RTC_NORMAL;
}

static int const_sect_mat_id_item(const struct hecmwST_local_mesh *global_mesh,
                                  struct hecmwST_local_mesh *local_mesh) {
  local_mesh->section->sect_mat_ID_item =
      global_mesh->section->sect_mat_ID_item;

  return RTC_NORMAL;
}

static int const_sect_i_index(const struct hecmwST_local_mesh *global_mesh,
                              struct hecmwST_local_mesh *local_mesh) {
  local_mesh->section->sect_I_index = global_mesh->section->sect_I_index;

  return RTC_NORMAL;
}

static int const_sect_i_item(const struct hecmwST_local_mesh *global_mesh,
                             struct hecmwST_local_mesh *local_mesh) {
  local_mesh->section->sect_I_item = global_mesh->section->sect_I_item;

  return RTC_NORMAL;
}

static int const_sect_r_index(const struct hecmwST_local_mesh *global_mesh,
                              struct hecmwST_local_mesh *local_mesh) {
  local_mesh->section->sect_R_index = global_mesh->section->sect_R_index;

  return RTC_NORMAL;
}

static int const_sect_r_item(const struct hecmwST_local_mesh *global_mesh,
                             struct hecmwST_local_mesh *local_mesh) {
  local_mesh->section->sect_R_item = global_mesh->section->sect_R_item;

  return RTC_NORMAL;
}

static int const_sect_info(const struct hecmwST_local_mesh *global_mesh,
                           struct hecmwST_local_mesh *local_mesh) {
  int rtc;

  HECMW_assert(global_mesh);
  HECMW_assert(local_mesh);
  HECMW_assert(global_mesh->section);
  HECMW_assert(local_mesh->section);

  rtc = const_n_sect(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_sect_type(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_sect_opt(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_sect_mat_id_index(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_sect_mat_id_item(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_sect_i_index(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_sect_i_item(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_sect_r_index(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_sect_r_item(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * - - - - - - - - - */

static int const_n_mat(const struct hecmwST_local_mesh *global_mesh,
                       struct hecmwST_local_mesh *local_mesh) {
  local_mesh->material->n_mat = global_mesh->material->n_mat;

  return RTC_NORMAL;
}

static int const_n_mat_item(const struct hecmwST_local_mesh *global_mesh,
                            struct hecmwST_local_mesh *local_mesh) {
  local_mesh->material->n_mat_item = global_mesh->material->n_mat_item;

  return RTC_NORMAL;
}

static int const_n_mat_subitem(const struct hecmwST_local_mesh *global_mesh,
                               struct hecmwST_local_mesh *local_mesh) {
  local_mesh->material->n_mat_subitem = global_mesh->material->n_mat_subitem;

  return RTC_NORMAL;
}

static int const_n_mat_table(const struct hecmwST_local_mesh *global_mesh,
                             struct hecmwST_local_mesh *local_mesh) {
  local_mesh->material->n_mat_table = global_mesh->material->n_mat_table;

  return RTC_NORMAL;
}

static int const_mat_name(const struct hecmwST_local_mesh *global_mesh,
                          struct hecmwST_local_mesh *local_mesh) {
  local_mesh->material->mat_name = global_mesh->material->mat_name;

  return RTC_NORMAL;
}

static int const_mat_item_index(const struct hecmwST_local_mesh *global_mesh,
                                struct hecmwST_local_mesh *local_mesh) {
  local_mesh->material->mat_item_index = global_mesh->material->mat_item_index;

  return RTC_NORMAL;
}

static int const_mat_subitem_index(const struct hecmwST_local_mesh *global_mesh,
                                   struct hecmwST_local_mesh *local_mesh) {
  local_mesh->material->mat_subitem_index =
      global_mesh->material->mat_subitem_index;

  return RTC_NORMAL;
}

static int const_mat_table_index(const struct hecmwST_local_mesh *global_mesh,
                                 struct hecmwST_local_mesh *local_mesh) {
  local_mesh->material->mat_table_index =
      global_mesh->material->mat_table_index;

  return RTC_NORMAL;
}

static int const_mat_val(const struct hecmwST_local_mesh *global_mesh,
                         struct hecmwST_local_mesh *local_mesh) {
  local_mesh->material->mat_val = global_mesh->material->mat_val;

  return RTC_NORMAL;
}

static int const_mat_temp(const struct hecmwST_local_mesh *global_mesh,
                          struct hecmwST_local_mesh *local_mesh) {
  local_mesh->material->mat_temp = global_mesh->material->mat_temp;

  return RTC_NORMAL;
}

static int const_mat_info(const struct hecmwST_local_mesh *global_mesh,
                          struct hecmwST_local_mesh *local_mesh) {
  int rtc;

  HECMW_assert(global_mesh);
  HECMW_assert(global_mesh->material);
  HECMW_assert(local_mesh);
  HECMW_assert(local_mesh->material);

  rtc = const_n_mat(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_n_mat_item(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_n_mat_subitem(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_n_mat_table(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_mat_name(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_mat_item_index(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_mat_subitem_index(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_mat_table_index(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_mat_val(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_mat_temp(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * - - - - - - - - - */

static int const_n_mpc(const struct hecmwST_local_mesh *global_mesh,
                       struct hecmwST_local_mesh *local_mesh,
                       const int *node_global2local, char *mpc_flag) {
  struct hecmwST_mpc *mpc_global = global_mesh->mpc;
  struct hecmwST_mpc *mpc_local  = local_mesh->mpc;
  int node, diff, evalsum, counter;
  int i, j;

  for (counter = 0, i = 0; i < mpc_global->n_mpc; i++) {
    diff = mpc_global->mpc_index[i + 1] - mpc_global->mpc_index[i];

    evalsum = 0;

    for (j = mpc_global->mpc_index[i]; j < mpc_global->mpc_index[i + 1]; j++) {
      node = mpc_global->mpc_item[j];
      if (node_global2local[node - 1] > 0) evalsum++;
    }

    if (evalsum == diff) {
      MASK_BIT(mpc_flag[i], MASK);
      counter++;
    }
  }
  mpc_local->n_mpc = counter;

  return RTC_NORMAL;
}

static int const_mpc_index(const struct hecmwST_local_mesh *global_mesh,
                           struct hecmwST_local_mesh *local_mesh,
                           const char *mpc_flag) {
  struct hecmwST_mpc *mpc_global = global_mesh->mpc;
  struct hecmwST_mpc *mpc_local  = local_mesh->mpc;
  int counter;
  int i;

  mpc_local->mpc_index = (int *)HECMW_calloc(mpc_local->n_mpc + 1, sizeof(int));
  if (local_mesh->mpc->mpc_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (counter = 0, i = 0; i < mpc_global->n_mpc; i++) {
    if (EVAL_BIT(mpc_flag[i], MASK)) {
      mpc_local->mpc_index[counter + 1] = mpc_local->mpc_index[counter] +
                                          mpc_global->mpc_index[i + 1] -
                                          mpc_global->mpc_index[i];
      counter++;
    }
  }
  HECMW_assert(counter == mpc_local->n_mpc);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_mpc_item(const struct hecmwST_local_mesh *global_mesh,
                          struct hecmwST_local_mesh *local_mesh,
                          const int *node_global2local, const char *mpc_flag) {
  struct hecmwST_mpc *mpc_global = global_mesh->mpc;
  struct hecmwST_mpc *mpc_local  = local_mesh->mpc;
  int mcounter, icounter;
  int i, j;

  mpc_local->mpc_item =
      (int *)HECMW_malloc(sizeof(int) * mpc_local->mpc_index[mpc_local->n_mpc]);
  if (mpc_local->mpc_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (mcounter = 0, icounter = 0, i = 0; i < mpc_global->n_mpc; i++) {
    if (EVAL_BIT(mpc_flag[i], MASK)) {
      for (j = mpc_global->mpc_index[i]; j < mpc_global->mpc_index[i + 1];
           j++) {
        mpc_local->mpc_item[mcounter++] =
            node_global2local[mpc_global->mpc_item[j] - 1];
      }
      HECMW_assert(mcounter == mpc_local->mpc_index[++icounter]);
    }
  }
  HECMW_assert(icounter == mpc_local->n_mpc);
  HECMW_assert(mcounter == mpc_local->mpc_index[mpc_local->n_mpc]);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_mpc_dof(const struct hecmwST_local_mesh *global_mesh,
                         struct hecmwST_local_mesh *local_mesh,
                         const char *mpc_flag) {
  struct hecmwST_mpc *mpc_global = global_mesh->mpc;
  struct hecmwST_mpc *mpc_local  = local_mesh->mpc;
  int mcounter, icounter;
  int i, j;

  mpc_local->mpc_dof =
      (int *)HECMW_malloc(sizeof(int) * mpc_local->mpc_index[mpc_local->n_mpc]);
  if (local_mesh->mpc->mpc_dof == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (mcounter = 0, icounter = 0, i = 0; i < mpc_global->n_mpc; i++) {
    if (EVAL_BIT(mpc_flag[i], MASK)) {
      for (j = mpc_global->mpc_index[i]; j < mpc_global->mpc_index[i + 1];
           j++) {
        mpc_local->mpc_dof[mcounter++] = mpc_global->mpc_dof[j];
      }
      HECMW_assert(mcounter == mpc_local->mpc_index[++icounter]);
    }
  }
  HECMW_assert(icounter == mpc_local->n_mpc);
  HECMW_assert(mcounter == mpc_local->mpc_index[mpc_local->n_mpc]);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_mpc_val(const struct hecmwST_local_mesh *global_mesh,
                         struct hecmwST_local_mesh *local_mesh,
                         const char *mpc_flag) {
  struct hecmwST_mpc *mpc_global = global_mesh->mpc;
  struct hecmwST_mpc *mpc_local  = local_mesh->mpc;
  int size;
  int mcounter, icounter;
  int i, j;

  size               = sizeof(double) * mpc_local->mpc_index[mpc_local->n_mpc];
  mpc_local->mpc_val = (double *)HECMW_malloc(size);
  if (local_mesh->mpc->mpc_val == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (mcounter = 0, icounter = 0, i = 0; i < mpc_global->n_mpc; i++) {
    if (EVAL_BIT(mpc_flag[i], MASK)) {
      for (j = mpc_global->mpc_index[i]; j < mpc_global->mpc_index[i + 1];
           j++) {
        mpc_local->mpc_val[mcounter++] = mpc_global->mpc_val[j];
      }
      HECMW_assert(mcounter == mpc_local->mpc_index[++icounter]);
    }
  }
  HECMW_assert(icounter == local_mesh->mpc->n_mpc);
  HECMW_assert(mcounter == local_mesh->mpc->mpc_index[local_mesh->mpc->n_mpc]);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_mpc_const(const struct hecmwST_local_mesh *global_mesh,
                           struct hecmwST_local_mesh *local_mesh,
                           const char *mpc_flag) {
  struct hecmwST_mpc *mpc_global = global_mesh->mpc;
  struct hecmwST_mpc *mpc_local  = local_mesh->mpc;
  int size;
  int icounter;
  int i;

  size                 = sizeof(double) * mpc_local->n_mpc;
  mpc_local->mpc_const = (double *)HECMW_malloc(size);
  if (local_mesh->mpc->mpc_const == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (icounter = 0, i = 0; i < mpc_global->n_mpc; i++) {
    if (EVAL_BIT(mpc_flag[i], MASK)) {
      mpc_local->mpc_const[icounter] = mpc_global->mpc_const[i];
      icounter++;
    }
  }
  HECMW_assert(icounter == local_mesh->mpc->n_mpc);

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_mpc_info(const struct hecmwST_local_mesh *global_mesh,
                          struct hecmwST_local_mesh *local_mesh,
                          const int *node_global2local) {
  char *mpc_flag = NULL;
  int rtc;

  HECMW_assert(global_mesh);
  HECMW_assert(global_mesh->mpc);
  HECMW_assert(local_mesh);
  HECMW_assert(local_mesh->mpc);
  HECMW_assert(node_global2local);

  if (global_mesh->mpc->n_mpc == 0) {
    init_struct_mpc(local_mesh);
    return RTC_NORMAL;
  }

  mpc_flag = (char *)HECMW_calloc(global_mesh->mpc->n_mpc, sizeof(char));
  if (mpc_flag == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  rtc = const_n_mpc(global_mesh, local_mesh, node_global2local, mpc_flag);
  if (rtc != RTC_NORMAL) goto error;

  if (local_mesh->mpc->n_mpc == 0) {
    init_struct_mpc(local_mesh);
    HECMW_free(mpc_flag);
    return RTC_NORMAL;
  }

  rtc = const_mpc_index(global_mesh, local_mesh, mpc_flag);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_mpc_item(global_mesh, local_mesh, node_global2local, mpc_flag);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_mpc_dof(global_mesh, local_mesh, mpc_flag);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_mpc_val(global_mesh, local_mesh, mpc_flag);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_mpc_const(global_mesh, local_mesh, mpc_flag);
  if (rtc != RTC_NORMAL) goto error;

  HECMW_free(mpc_flag);

  return RTC_NORMAL;

error:
  HECMW_free(mpc_flag);

  return RTC_ERROR;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * - - - - - - - - - */

static int const_n_amp(const struct hecmwST_local_mesh *global_mesh,
                       struct hecmwST_local_mesh *local_mesh) {
  local_mesh->amp->n_amp = global_mesh->amp->n_amp;

  return RTC_NORMAL;
}

static int const_amp_name(const struct hecmwST_local_mesh *global_mesh,
                          struct hecmwST_local_mesh *local_mesh) {
  local_mesh->amp->amp_name = global_mesh->amp->amp_name;

  return RTC_NORMAL;
}

static int const_amp_type_definition(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh) {
  local_mesh->amp->amp_type_definition = global_mesh->amp->amp_type_definition;

  return RTC_NORMAL;
}

static int const_amp_type_time(const struct hecmwST_local_mesh *global_mesh,
                               struct hecmwST_local_mesh *local_mesh) {
  local_mesh->amp->amp_type_time = global_mesh->amp->amp_type_time;

  return RTC_NORMAL;
}

static int const_amp_type_value(const struct hecmwST_local_mesh *global_mesh,
                                struct hecmwST_local_mesh *local_mesh) {
  local_mesh->amp->amp_type_value = global_mesh->amp->amp_type_value;

  return RTC_NORMAL;
}

static int const_amp_index(const struct hecmwST_local_mesh *global_mesh,
                           struct hecmwST_local_mesh *local_mesh) {
  local_mesh->amp->amp_index = global_mesh->amp->amp_index;

  return RTC_NORMAL;
}

static int const_amp_val(const struct hecmwST_local_mesh *global_mesh,
                         struct hecmwST_local_mesh *local_mesh) {
  local_mesh->amp->amp_val = global_mesh->amp->amp_val;

  return RTC_NORMAL;
}

static int const_amp_table(const struct hecmwST_local_mesh *global_mesh,
                           struct hecmwST_local_mesh *local_mesh) {
  local_mesh->amp->amp_table = global_mesh->amp->amp_table;

  return RTC_NORMAL;
}

static int const_amp_info(const struct hecmwST_local_mesh *global_mesh,
                          struct hecmwST_local_mesh *local_mesh) {
  int rtc;

  HECMW_assert(global_mesh);
  HECMW_assert(global_mesh->amp);
  HECMW_assert(local_mesh);
  HECMW_assert(local_mesh->amp);

  if (global_mesh->amp->n_amp == 0) {
    init_struct_amp(local_mesh);
    return RTC_NORMAL;
  }

  rtc = const_n_amp(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_amp_name(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_amp_type_definition(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_amp_type_time(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_amp_type_value(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_amp_index(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_amp_val(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_amp_table(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * - - - - - - - - - */

static int *const_node_grp_mask_eqn(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh, const int *node_global2local,
    int eqn_block_idx) {
  struct hecmwST_node_grp *node_group_global = global_mesh->node_group;
  int *n_eqn_item                            = NULL;
  int diff, evalsum;
  int i, j, is, ie, js;

  is = node_group_global->grp_index[eqn_block_idx];
  ie = node_group_global->grp_index[eqn_block_idx + 1];

  n_eqn_item = (int *)HECMW_malloc(sizeof(int) * (ie - is));
  if (n_eqn_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (js = 0, i = 0; i < ie - is; i++) {
    diff = node_group_global->grp_item[is + i] - js;
    for (evalsum = 0, j = js; j < node_group_global->grp_item[is + i]; j++) {
      if (node_global2local[j] > 0 &&
          node_global2local[j] <= local_mesh->nn_internal)
        evalsum++;
    }

    if (evalsum) {
      HECMW_assert(evalsum == diff);
      n_eqn_item[i] = diff;
    } else {
      n_eqn_item[i] = 0;
    }

    js = node_group_global->grp_item[is + i];
  }

  return n_eqn_item;

error:
  return NULL;
}

static int const_node_n_grp(const struct hecmwST_local_mesh *global_mesh,
                            struct hecmwST_local_mesh *local_mesh) {
  local_mesh->node_group->n_grp = global_mesh->node_group->n_grp;

  return RTC_NORMAL;
}

static int const_node_grp_name(const struct hecmwST_local_mesh *global_mesh,
                               struct hecmwST_local_mesh *local_mesh) {
  local_mesh->node_group->grp_name = global_mesh->node_group->grp_name;

  return RTC_NORMAL;
}

static int const_node_grp_index(const struct hecmwST_local_mesh *global_mesh,
                                struct hecmwST_local_mesh *local_mesh,
                                const int *node_global2local,
                                const int *n_eqn_item, int eqn_block_idx) {
  struct hecmwST_node_grp *node_group_global = global_mesh->node_group;
  struct hecmwST_node_grp *node_group_local  = local_mesh->node_group;
  int node;
  int counter, diff;
  int i, j;

  node_group_local->grp_index =
      (int *)HECMW_calloc(node_group_local->n_grp + 1, sizeof(int));
  if (node_group_local->grp_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (counter = 0, i = 0; i < node_group_global->n_grp; i++) {
    if (i != eqn_block_idx) {
      for (j = node_group_global->grp_index[i];
           j < node_group_global->grp_index[i + 1]; j++) {
        node = node_group_global->grp_item[j];
        if (node_global2local[node - 1]) counter++;
      }

    } else {
      diff =
          node_group_global->grp_index[i + 1] - node_group_global->grp_index[i];
      for (j = 0; j < diff; j++) {
        if (n_eqn_item[j] > 0) counter++;
      }
    }

    node_group_local->grp_index[i + 1] = counter;
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*K. Inagaki */
static int const_node_grp_index_mod(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh, const int *node_global2local,
    const int *n_eqn_item, int eqn_block_idx, int domain) {
  struct hecmwST_node_grp *node_group_global = global_mesh->node_group;
  struct hecmwST_node_grp *node_group_local  = local_mesh->node_group;
  int node;
  int counter, diff;
  int i, j;

  node_group_local->grp_index =
      (int *)HECMW_calloc(node_group_local->n_grp + 1, sizeof(int));
  if (node_group_local->grp_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (counter = 0, i = 0; i < node_group_global->n_grp; i++) {
    if (i != eqn_block_idx) {
      if (node_group_global->grp_index[i + 1] -
              node_group_global->grp_index[i] ==
          global_mesh->n_node) {
        counter += n_int_nlist[domain];
        counter += n_bnd_nlist[2 * domain + 1] - n_bnd_nlist[2 * domain];
      } else {
        counter += ngrp_idx[domain][i + 1] - ngrp_idx[domain][i];
        /*
        for( j=node_group_global->grp_index[i];
        j<node_group_global->grp_index[i+1]; j++ ) {
            node = node_group_global->grp_item[j];
            if( node_global2local[node-1] )  counter++;
        }
        */
      }

    } else {
      diff =
          node_group_global->grp_index[i + 1] - node_group_global->grp_index[i];
      for (j = 0; j < diff; j++) {
        if (n_eqn_item[j] > 0) counter++;
      }
    }

    node_group_local->grp_index[i + 1] = counter;
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_node_grp_item(const struct hecmwST_local_mesh *global_mesh,
                               struct hecmwST_local_mesh *local_mesh,
                               const int *node_global2local,
                               const int *n_eqn_item, int eqn_block_idx) {
  struct hecmwST_node_grp *node_group_global = global_mesh->node_group;
  struct hecmwST_node_grp *node_group_local  = local_mesh->node_group;
  int node;
  int size;
  int counter;
  int i, j, k, js, je, ks, ls;

  size = sizeof(int) * node_group_local->grp_index[node_group_local->n_grp];
  node_group_local->grp_item = (int *)HECMW_malloc(size);
  if (node_group_local->grp_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (counter = 0, i = 0; i < node_group_global->n_grp; i++) {
    if (i != eqn_block_idx) {
      for (j = node_group_global->grp_index[i];
           j < node_group_global->grp_index[i + 1]; j++) {
        node = node_group_global->grp_item[j];
        if (node_global2local[node - 1]) {
          node_group_local->grp_item[counter++] = node_global2local[node - 1];
        }
      }

    } else {
      js = node_group_global->grp_index[i];
      je = node_group_global->grp_index[i + 1];
      for (ks = 0, ls = 0, j = js; j < je; j++) {
        if (n_eqn_item[j - js]) {
          HECMW_assert(n_eqn_item[j - js] ==
                       node_group_global->grp_item[j] - ks);
          node_group_local->grp_item[counter] = ls + n_eqn_item[j - js];

          for (k = ks; k < node_group_global->grp_item[j]; k++) {
            HECMW_assert(ls < node_global2local[k] &&
                         node_global2local[k] <=
                             node_group_local->grp_item[counter]);
          }
          ls = node_group_local->grp_item[counter];
          counter++;
        }
        ks = node_group_global->grp_item[j];
      }
    }
    HECMW_assert(counter == node_group_local->grp_index[i + 1]);
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*K. Inagaki */
static int const_node_grp_item_mod(const struct hecmwST_local_mesh *global_mesh,
                                   struct hecmwST_local_mesh *local_mesh,
                                   const int *node_global2local,
                                   const int *n_eqn_item, int eqn_block_idx,
                                   int domain) {
  struct hecmwST_node_grp *node_group_global = global_mesh->node_group;
  struct hecmwST_node_grp *node_group_local  = local_mesh->node_group;
  int node;
  int size;
  int counter;
  int i, j, k, js, je, ks, ls;
  int idx1, idx2, node1, node2, n_int, n_bnd, n_out, maxn;

  size = sizeof(int) * node_group_local->grp_index[node_group_local->n_grp];
  node_group_local->grp_item = (int *)HECMW_malloc(size);
  if (node_group_local->grp_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  n_int = n_int_nlist[domain];
  n_bnd = n_bnd_nlist[2 * domain];
  n_out = n_bnd_nlist[2 * domain + 1] - n_bnd_nlist[2 * domain];
  maxn  = global_mesh->n_node + 1;

  for (counter = 0, i = 0; i < node_group_global->n_grp; i++) {
    if (i != eqn_block_idx) {
      if (node_group_global->grp_index[i + 1] -
              node_group_global->grp_index[i] ==
          global_mesh->n_node) {
        idx1  = 0;
        idx2  = 0;
        node1 = (n_int == 0) ? maxn : int_nlist[domain][0];
        node2 = (n_out == 0) ? maxn : bnd_nlist[domain][n_bnd];
        for (j = 0; j < n_int + n_out; j++) {
          if (node1 < node2) {
            node_group_local->grp_item[counter++] =
                node_global2local[node1 - 1];
            idx1++;
            node1 = (idx1 == n_int) ? maxn : int_nlist[domain][idx1];
          } else {
            node_group_local->grp_item[counter++] =
                node_global2local[node2 - 1];
            idx2++;
            node2 = (idx2 == n_out) ? maxn : bnd_nlist[domain][idx2 + n_bnd];
          }
        }
      } else {
        if (ngrp_idx[domain][i + 1] - ngrp_idx[domain][i] == 0) continue;
        for (j = ngrp_idx[domain][i]; j < ngrp_idx[domain][i + 1]; j++) {
          node                                  = ngrp_item[domain][j];
          node_group_local->grp_item[counter++] = node_global2local[node - 1];
        }
      }
    } else {
      js = node_group_global->grp_index[i];
      je = node_group_global->grp_index[i + 1];
      for (ks = 0, ls = 0, j = js; j < je; j++) {
        if (n_eqn_item[j - js]) {
          HECMW_assert(n_eqn_item[j - js] ==
                       node_group_global->grp_item[j] - ks);
          node_group_local->grp_item[counter] = ls + n_eqn_item[j - js];

          for (k = ks; k < node_group_global->grp_item[j]; k++) {
            HECMW_assert(ls < node_global2local[k] &&
                         node_global2local[k] <=
                             node_group_local->grp_item[counter]);
          }
          ls = node_group_local->grp_item[counter];
          counter++;
        }
        ks = node_group_global->grp_item[j];
      }
    }
    HECMW_assert(counter == node_group_local->grp_index[i + 1]);
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_node_grp_info(const struct hecmwST_local_mesh *global_mesh,
                               struct hecmwST_local_mesh *local_mesh,
                               const int *node_global2local,
                               int current_domain) {
  int *n_eqn_item = NULL;
  int eqn_block_idx;
  int rtc;

  HECMW_assert(global_mesh);
  HECMW_assert(global_mesh->node_group);
  HECMW_assert(local_mesh);
  HECMW_assert(local_mesh->node_group);
  HECMW_assert(node_global2local);

  if (global_mesh->node_group->n_grp == 0) {
    init_struct_node_grp(local_mesh);
    return RTC_NORMAL;
  }

  eqn_block_idx = search_eqn_block_idx(global_mesh);

  if (eqn_block_idx >= 0) {
    n_eqn_item = const_node_grp_mask_eqn(global_mesh, local_mesh,
                                         node_global2local, eqn_block_idx);
    if (n_eqn_item == NULL) goto error;
  }

  rtc = const_node_n_grp(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_node_grp_name(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  if (is_spdup_available(global_mesh)) {
    rtc = const_node_grp_index_mod(global_mesh, local_mesh, node_global2local,
                                   n_eqn_item, eqn_block_idx, current_domain);
    if (rtc != RTC_NORMAL) goto error;
    rtc = const_node_grp_item_mod(global_mesh, local_mesh, node_global2local,
                                  n_eqn_item, eqn_block_idx, current_domain);
    if (rtc != RTC_NORMAL) goto error;

  } else {
    rtc = const_node_grp_index(global_mesh, local_mesh, node_global2local,
                               n_eqn_item, eqn_block_idx);
    if (rtc != RTC_NORMAL) goto error;
    rtc = const_node_grp_item(global_mesh, local_mesh, node_global2local,
                              n_eqn_item, eqn_block_idx);
    if (rtc != RTC_NORMAL) goto error;
  }

  HECMW_free(n_eqn_item);

  return RTC_NORMAL;

error:
  HECMW_free(n_eqn_item);

  return RTC_ERROR;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * - - - - - - - - - */

static int const_elem_n_grp(const struct hecmwST_local_mesh *global_mesh,
                            struct hecmwST_local_mesh *local_mesh) {
  local_mesh->elem_group->n_grp = global_mesh->elem_group->n_grp;

  return RTC_NORMAL;
}

static int const_elem_grp_name(const struct hecmwST_local_mesh *global_mesh,
                               struct hecmwST_local_mesh *local_mesh) {
  local_mesh->elem_group->grp_name = global_mesh->elem_group->grp_name;

  return RTC_NORMAL;
}

static int const_elem_grp_index(const struct hecmwST_local_mesh *global_mesh,
                                struct hecmwST_local_mesh *local_mesh,
                                const int *elem_global2local) {
  struct hecmwST_elem_grp *elem_group_global = global_mesh->elem_group;
  struct hecmwST_elem_grp *elem_group_local  = local_mesh->elem_group;
  int elem;
  int counter;
  int i, j;

  elem_group_local->grp_index =
      (int *)HECMW_calloc(elem_group_local->n_grp + 1, sizeof(int));
  if (elem_group_local->grp_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (counter = 0, i = 0; i < elem_group_global->n_grp; i++) {
    for (j = elem_group_global->grp_index[i];
         j < elem_group_global->grp_index[i + 1]; j++) {
      elem = elem_group_global->grp_item[j];
      if (elem_global2local[elem - 1]) counter++;
    }
    elem_group_local->grp_index[i + 1] = counter;
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*K. Inagaki */
static int const_elem_grp_index_mod(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh, const int *elem_global2local,
    int domain) {
  struct hecmwST_elem_grp *elem_group_global = global_mesh->elem_group;
  struct hecmwST_elem_grp *elem_group_local  = local_mesh->elem_group;
  int elem;
  int counter;
  int i, j, idx1, idx2, elem1, elem2;

  elem_group_local->grp_index =
      (int *)HECMW_calloc(elem_group_local->n_grp + 1, sizeof(int));
  if (elem_group_local->grp_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (counter = 0, i = 0; i < elem_group_global->n_grp; i++) {
    if (elem_group_global->grp_index[i + 1] - elem_group_global->grp_index[i] ==
        global_mesh->n_elem) {
      counter += n_int_elist[domain];
      counter += n_bnd_elist[2 * domain + 1] - n_bnd_elist[2 * domain];
    } else {
      counter += egrp_idx[domain][i + 1] - egrp_idx[domain][i];
    }
    elem_group_local->grp_index[i + 1] = counter;
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_elem_grp_item(const struct hecmwST_local_mesh *global_mesh,
                               struct hecmwST_local_mesh *local_mesh,
                               const int *elem_global2local) {
  struct hecmwST_elem_grp *elem_group_global = global_mesh->elem_group;
  struct hecmwST_elem_grp *elem_group_local  = local_mesh->elem_group;
  int elem;
  int size;
  int counter;
  int i, j;

  size = sizeof(int) * elem_group_local->grp_index[elem_group_local->n_grp];
  elem_group_local->grp_item = (int *)HECMW_malloc(size);
  if (local_mesh->elem_group->grp_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (counter = 0, i = 0; i < elem_group_global->n_grp; i++) {
    for (j = elem_group_global->grp_index[i];
         j < elem_group_global->grp_index[i + 1]; j++) {
      elem = elem_group_global->grp_item[j];
      if (elem_global2local[elem - 1]) {
        elem_group_local->grp_item[counter++] = elem_global2local[elem - 1];
      }
    }
    HECMW_assert(counter == elem_group_local->grp_index[i + 1]);
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*K. Inagaki */
static int const_elem_grp_item_mod(const struct hecmwST_local_mesh *global_mesh,
                                   struct hecmwST_local_mesh *local_mesh,
                                   const int *elem_global2local, int domain) {
  struct hecmwST_elem_grp *elem_group_global = global_mesh->elem_group;
  struct hecmwST_elem_grp *elem_group_local  = local_mesh->elem_group;
  int elem;
  int size;
  int counter;
  int i, j, idx1, idx2, elem1, elem2, n_int, n_bnd, n_out, maxe;

  size = sizeof(int) * elem_group_local->grp_index[elem_group_local->n_grp];
  elem_group_local->grp_item = (int *)HECMW_malloc(size);
  if (local_mesh->elem_group->grp_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  n_int = n_int_elist[domain];
  n_bnd = n_bnd_elist[2 * domain];
  n_out = n_bnd_elist[2 * domain + 1] - n_bnd_elist[2 * domain];
  maxe  = global_mesh->n_elem + 1;

  for (counter = 0, i = 0; i < elem_group_global->n_grp; i++) {
    if (elem_group_global->grp_index[i + 1] - elem_group_global->grp_index[i] ==
        global_mesh->n_elem) {
      elem1 = (n_int == 0) ? maxe : int_elist[domain][0];
      elem2 = (n_out == 0) ? maxe : bnd_elist[domain][n_bnd];
      for (idx1 = 0, idx2 = 0, j = 0; j < n_int + n_out; j++) {
        if (elem1 < elem2) {
          elem_group_local->grp_item[counter++] = elem_global2local[elem1 - 1];
          idx1++;
          elem1 = (idx1 == n_int) ? maxe : int_elist[domain][idx1];
        } else {
          elem_group_local->grp_item[counter++] = elem_global2local[elem2 - 1];
          idx2++;
          elem2 = (idx2 == n_out) ? maxe : bnd_elist[domain][idx2 + n_bnd];
        }
      }
    } else {
      if (egrp_idx[domain][i + 1] - egrp_idx[domain][i] == 0) continue;
      for (j = egrp_idx[domain][i]; j < egrp_idx[domain][i + 1]; j++) {
        elem                                  = egrp_item[domain][j];
        elem_group_local->grp_item[counter++] = elem_global2local[elem - 1];
      }
    }
    HECMW_assert(counter == elem_group_local->grp_index[i + 1]);
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_elem_grp_info(const struct hecmwST_local_mesh *global_mesh,
                               struct hecmwST_local_mesh *local_mesh,
                               const int *elem_global2local,
                               int current_domain) {
  int rtc;

  HECMW_assert(global_mesh);
  HECMW_assert(global_mesh->elem_group);
  HECMW_assert(local_mesh);
  HECMW_assert(local_mesh->elem_group);
  HECMW_assert(elem_global2local);

  if (global_mesh->elem_group->n_grp == 0) {
    init_struct_elem_grp(local_mesh);
    return RTC_NORMAL;
  }

  rtc = const_elem_n_grp(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_elem_grp_name(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  if (is_spdup_available(global_mesh)) {
    rtc = const_elem_grp_index_mod(global_mesh, local_mesh, elem_global2local,
                                   current_domain);
    if (rtc != RTC_NORMAL) goto error;
    rtc = const_elem_grp_item_mod(global_mesh, local_mesh, elem_global2local,
                                  current_domain);
    if (rtc != RTC_NORMAL) goto error;

  } else {
    rtc = const_elem_grp_index(global_mesh, local_mesh, elem_global2local);
    if (rtc != RTC_NORMAL) goto error;
    rtc = const_elem_grp_item(global_mesh, local_mesh, elem_global2local);
    if (rtc != RTC_NORMAL) goto error;
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * - - - - - - - - - */

static int const_surf_n_grp(const struct hecmwST_local_mesh *global_mesh,
                            struct hecmwST_local_mesh *local_mesh) {
  local_mesh->surf_group->n_grp = global_mesh->surf_group->n_grp;

  return RTC_NORMAL;
}

static int const_surf_grp_name(const struct hecmwST_local_mesh *global_mesh,
                               struct hecmwST_local_mesh *local_mesh) {
  local_mesh->surf_group->grp_name = global_mesh->surf_group->grp_name;

  return RTC_NORMAL;
}

static int const_surf_grp_index(const struct hecmwST_local_mesh *global_mesh,
                                struct hecmwST_local_mesh *local_mesh,
                                const int *elem_global2local) {
  struct hecmwST_surf_grp *surf_group_global = global_mesh->surf_group;
  struct hecmwST_surf_grp *surf_group_local  = local_mesh->surf_group;
  int elem;
  int counter;
  int i, j;

  surf_group_local->grp_index =
      (int *)HECMW_calloc(surf_group_local->n_grp + 1, sizeof(int));
  if (surf_group_local->grp_index == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (counter = 0, i = 0; i < surf_group_global->n_grp; i++) {
    for (j = surf_group_global->grp_index[i];
         j < surf_group_global->grp_index[i + 1]; j++) {
      elem = surf_group_global->grp_item[2 * j];
      if (elem_global2local[elem - 1]) counter++;
    }
    surf_group_local->grp_index[i + 1] = counter;
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_surf_grp_item(const struct hecmwST_local_mesh *global_mesh,
                               struct hecmwST_local_mesh *local_mesh,
                               const int *elem_global2local) {
  struct hecmwST_surf_grp *surf_group_global = global_mesh->surf_group;
  struct hecmwST_surf_grp *surf_group_local  = local_mesh->surf_group;
  int elem, surf;
  int size;
  int counter;
  int i, j;

  size = sizeof(int) * surf_group_local->grp_index[surf_group_local->n_grp] * 2;
  surf_group_local->grp_item = (int *)HECMW_malloc(size);
  if (surf_group_local->grp_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (counter = 0, i = 0; i < surf_group_global->n_grp; i++) {
    for (j = surf_group_global->grp_index[i];
         j < surf_group_global->grp_index[i + 1]; j++) {
      elem = surf_group_global->grp_item[2 * j];
      surf = surf_group_global->grp_item[2 * j + 1];
      if (elem_global2local[elem - 1]) {
        surf_group_local->grp_item[2 * counter] = elem_global2local[elem - 1];
        surf_group_local->grp_item[2 * counter + 1] = surf;
        counter++;
      }
    }
    HECMW_assert(counter == surf_group_local->grp_index[i + 1]);
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_surf_grp_info(const struct hecmwST_local_mesh *global_mesh,
                               struct hecmwST_local_mesh *local_mesh,
                               const int *elem_global2local) {
  int rtc;

  HECMW_assert(global_mesh);
  HECMW_assert(global_mesh->surf_group);
  HECMW_assert(local_mesh);
  HECMW_assert(local_mesh->surf_group);
  HECMW_assert(elem_global2local);

  if (global_mesh->surf_group->n_grp == 0) {
    init_struct_surf_grp(local_mesh);
    return RTC_NORMAL;
  }

  rtc = const_surf_n_grp(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_surf_grp_name(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_surf_grp_index(global_mesh, local_mesh, elem_global2local);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_surf_grp_item(global_mesh, local_mesh, elem_global2local);
  if (rtc != RTC_NORMAL) goto error;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * - - - - - - - - - */

static int const_contact_pair_n_pair(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh) {
  local_mesh->contact_pair->n_pair = global_mesh->contact_pair->n_pair;

  return RTC_NORMAL;
}

static int const_contact_pair_name(const struct hecmwST_local_mesh *global_mesh,
                                   struct hecmwST_local_mesh *local_mesh) {
  local_mesh->contact_pair->name = global_mesh->contact_pair->name;

  return RTC_NORMAL;
}

static int const_contact_pair_type(const struct hecmwST_local_mesh *global_mesh,
                                   struct hecmwST_local_mesh *local_mesh) {
  struct hecmwST_contact_pair *cpair_global = global_mesh->contact_pair;
  struct hecmwST_contact_pair *cpair_local  = local_mesh->contact_pair;
  int i;

  cpair_local->type = (int *)HECMW_calloc(cpair_local->n_pair, sizeof(int));
  if (cpair_local->type == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (i = 0; i < cpair_global->n_pair; i++) {
    cpair_local->type[i] = cpair_global->type[i];
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_contact_pair_slave_grp_id(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh) {
  struct hecmwST_contact_pair *cpair_global = global_mesh->contact_pair;
  struct hecmwST_contact_pair *cpair_local  = local_mesh->contact_pair;
  int i;

  cpair_local->slave_grp_id =
      (int *)HECMW_calloc(cpair_local->n_pair, sizeof(int));
  if (cpair_local->slave_grp_id == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (i = 0; i < cpair_global->n_pair; i++) {
    cpair_local->slave_grp_id[i] = cpair_global->slave_grp_id[i];
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_contact_pair_slave_orisgrp_id(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh) {
  struct hecmwST_contact_pair *cpair_global = global_mesh->contact_pair;
  struct hecmwST_contact_pair *cpair_local  = local_mesh->contact_pair;
  int i;

  cpair_local->slave_orisgrp_id =
      (int *)HECMW_calloc(cpair_local->n_pair, sizeof(int));
  if (cpair_local->slave_orisgrp_id == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (i = 0; i < cpair_global->n_pair; i++) {
    cpair_local->slave_orisgrp_id[i] = cpair_global->slave_orisgrp_id[i];
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_contact_pair_master_grp_id(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_local_mesh *local_mesh) {
  struct hecmwST_contact_pair *cpair_global = global_mesh->contact_pair;
  struct hecmwST_contact_pair *cpair_local  = local_mesh->contact_pair;
  int i;

  cpair_local->master_grp_id =
      (int *)HECMW_calloc(cpair_local->n_pair, sizeof(int));
  if (cpair_local->master_grp_id == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (i = 0; i < cpair_global->n_pair; i++) {
    cpair_local->master_grp_id[i] = cpair_global->master_grp_id[i];
  }

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

static int const_contact_pair_info(const struct hecmwST_local_mesh *global_mesh,
                                   struct hecmwST_local_mesh *local_mesh) {
  int rtc;

  HECMW_assert(global_mesh);
  HECMW_assert(global_mesh->contact_pair);
  HECMW_assert(local_mesh);
  HECMW_assert(local_mesh->contact_pair);

  if (global_mesh->contact_pair->n_pair == 0) {
    init_struct_contact_pair(local_mesh);
    return RTC_NORMAL;
  }

  rtc = const_contact_pair_n_pair(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_contact_pair_name(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_contact_pair_type(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_contact_pair_slave_grp_id(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_contact_pair_slave_orisgrp_id(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_contact_pair_master_grp_id(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  return RTC_NORMAL;

error:
  return RTC_ERROR;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * - - - - - - - - - */

static int const_local_data(const struct hecmwST_local_mesh *global_mesh,
                            struct hecmwST_local_mesh *local_mesh,
                            const struct hecmw_part_cont_data *cont_data,
                            const char *node_flag, const char *elem_flag,
                            int *node_global2local, int *elem_global2local,
                            int current_domain) {
  int *node_local2global = NULL;
  int *elem_local2global = NULL;
  int rtc, i;

  HECMW_log(HECMW_LOG_DEBUG, "Starting creation of local mesh data...\n");

  rtc = set_node_global2local(global_mesh, local_mesh, node_global2local,
                              node_flag, current_domain);
  if (rtc != RTC_NORMAL) goto error;

  node_local2global = (int *)HECMW_calloc(local_mesh->n_node, sizeof(int));
  if (node_local2global == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  if (is_spdup_available(global_mesh)) {
    rtc = set_node_local2global_mod(global_mesh, local_mesh, node_global2local,
                                    node_local2global, current_domain);
  } else {
    rtc = set_node_local2global(global_mesh, local_mesh, node_global2local,
                                node_local2global);
  }

  if (rtc != RTC_NORMAL) goto error;

  rtc = set_elem_global2local(global_mesh, local_mesh, elem_global2local,
                              elem_flag, current_domain);

  if (rtc != RTC_NORMAL) goto error;

  elem_local2global = (int *)HECMW_calloc(local_mesh->n_elem, sizeof(int));
  if (elem_local2global == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  if (is_spdup_available(global_mesh)) {
    rtc = set_elem_local2global_mod(global_mesh, local_mesh, elem_global2local,
                                    elem_local2global, current_domain);
  } else {
    rtc = set_elem_local2global(global_mesh, local_mesh, elem_global2local,
                                elem_local2global);
  }

  if (rtc != RTC_NORMAL) goto error;

  rtc = const_global_info(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_node_info(global_mesh, local_mesh, node_local2global, node_flag,
                        current_domain);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_elem_info(global_mesh, local_mesh, node_global2local,
                        elem_global2local, elem_local2global, current_domain);

  if (rtc != RTC_NORMAL) goto error;
  rtc = const_comm_info(global_mesh, local_mesh, node_global2local,
                        elem_global2local, current_domain);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_adapt_info(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_sect_info(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_mat_info(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_mpc_info(global_mesh, local_mesh, node_global2local);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_amp_info(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_node_grp_info(global_mesh, local_mesh, node_global2local,
                            current_domain);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_elem_grp_info(global_mesh, local_mesh, elem_global2local,
                            current_domain);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_surf_grp_info(global_mesh, local_mesh, elem_global2local);
  if (rtc != RTC_NORMAL) goto error;

  rtc = const_contact_pair_info(global_mesh, local_mesh);
  if (rtc != RTC_NORMAL) goto error;

  rtc = clear_node_global2local(global_mesh, local_mesh, node_global2local,
                                current_domain);
  rtc = clear_elem_global2local(global_mesh, local_mesh, elem_global2local,
                                current_domain);

  HECMW_free(node_local2global);
  HECMW_free(elem_local2global);

  HECMW_log(HECMW_LOG_DEBUG, "Creation of local mesh data done\n");

  return RTC_NORMAL;

error:
  HECMW_free(node_local2global);
  HECMW_free(elem_local2global);
  clean_struct_local_mesh(local_mesh);

  return RTC_ERROR;
}

/*==================================================================================================

  print UCD format data

==================================================================================================*/

static int print_ucd_entire_set_node_data(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_result_data *result_data, const char *node_flag) {
  int size;
  int nn_item;
  int i;

  result_data->nn_component = 1;

  result_data->nn_dof =
      (int *)HECMW_malloc(sizeof(int) * result_data->nn_component);
  if (result_data->nn_dof == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  result_data->nn_dof[0] = 1;

  result_data->node_label =
      (char **)HECMW_malloc(sizeof(char *) * result_data->nn_component);
  if (result_data->node_label == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  } else {
    for (i = 0; i < result_data->nn_component; i++) {
      result_data->node_label[i] = NULL;
    }
  }
  for (i = 0; i < result_data->nn_component; i++) {
    result_data->node_label[i] =
        (char *)HECMW_malloc(sizeof(char) * (HECMW_NAME_LEN + 1));
    if (result_data->node_label[i] == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    }
  }
  strcpy(result_data->node_label[0], "rank_of_node");

  for (nn_item = 0, i = 0; i < result_data->nn_component; i++) {
    nn_item += result_data->nn_dof[i];
  }

  size                       = sizeof(double) * nn_item * global_mesh->n_node;
  result_data->node_val_item = (double *)HECMW_malloc(size);
  if (result_data->node_val_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  switch (global_mesh->hecmw_flag_parttype) {
    case HECMW_FLAG_PARTTYPE_NODEBASED:
      for (i = 0; i < global_mesh->n_node; i++) {
        result_data->node_val_item[i] = (double)global_mesh->node_ID[2 * i + 1];
      }
      break;

    case HECMW_FLAG_PARTTYPE_ELEMBASED:
      for (i = 0; i < global_mesh->n_node; i++) {
        if (EVAL_BIT(node_flag[i], OVERLAP)) {
          result_data->node_val_item[i] =
              (double)global_mesh->n_subdomain + 2.0;
        } else {
          result_data->node_val_item[i] =
              (double)global_mesh->node_ID[2 * i + 1];
        }
      }
      break;

    default:
      HECMW_set_error(HECMW_PART_E_INVALID_PTYPE, "%d",
                      global_mesh->hecmw_flag_parttype);
      goto error;
  }

  return RTC_NORMAL;

error:
  free_struct_result_data(result_data);

  return RTC_ERROR;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * - - - - - - - - - */

static int print_ucd_entire_set_elem_data(
    const struct hecmwST_local_mesh *global_mesh,
    struct hecmwST_result_data *result_data, const char *elem_flag) {
  int size;
  int ne_item;
  int i;

  result_data->ne_component = 1;

  result_data->ne_dof =
      (int *)HECMW_malloc(sizeof(int) * result_data->ne_component);
  if (result_data->ne_dof == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  result_data->ne_dof[0] = 1;

  result_data->elem_label =
      (char **)HECMW_malloc(sizeof(char *) * result_data->ne_component);
  if (result_data->elem_label == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  } else {
    for (i = 0; i < result_data->ne_component; i++) {
      result_data->elem_label[i] = NULL;
    }
  }
  for (i = 0; i < result_data->ne_component; i++) {
    result_data->elem_label[i] =
        (char *)HECMW_malloc(sizeof(char) * (HECMW_NAME_LEN + 1));
    if (result_data->elem_label[i] == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    }
  }
  strcpy(result_data->elem_label[0], "partitioning_image");

  /* modify element information*/
  for (i = 0; i < global_mesh->n_elem; i++) {
    switch (global_mesh->elem_type[i]) {
      case HECMW_ETYPE_SHT6:
        global_mesh->elem_type[i] = HECMW_ETYPE_SHT1;
        break;

      case HECMW_ETYPE_SHQ8:
        global_mesh->elem_type[i] = HECMW_ETYPE_SHQ1;
        break;

      case HECMW_ETYPE_BEM3:
        global_mesh->elem_type[i] = HECMW_ETYPE_ROD1;
        break;

      case HECMW_ETYPE_ROD31:
        global_mesh->elem_type[i] = HECMW_ETYPE_ROD1;
        break;
    }
  }

  for (ne_item = 0, i = 0; i < result_data->ne_component; i++) {
    ne_item += result_data->ne_dof[i];
  }

  size                       = sizeof(double) * ne_item * global_mesh->n_elem;
  result_data->elem_val_item = (double *)HECMW_malloc(size);
  if (result_data->elem_val_item == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  switch (global_mesh->hecmw_flag_parttype) {
    case HECMW_FLAG_PARTTYPE_NODEBASED:
      for (i = 0; i < global_mesh->n_elem; i++) {
        if (EVAL_BIT(elem_flag[i], OVERLAP)) {
          result_data->elem_val_item[i] =
              (double)global_mesh->n_subdomain + 2.0;
        } else {
          result_data->elem_val_item[i] =
              (double)global_mesh->elem_ID[2 * i + 1];
        }
      }
      break;

    case HECMW_FLAG_PARTTYPE_ELEMBASED:
      for (i = 0; i < global_mesh->n_elem; i++) {
        result_data->elem_val_item[i] = (double)global_mesh->elem_ID[2 * i + 1];
      }
      break;

    default:
      HECMW_set_error(HECMW_PART_E_INVALID_PTYPE, "%d",
                      global_mesh->hecmw_flag_parttype);
      goto error;
  }

  return RTC_NORMAL;

error:
  free_struct_result_data(result_data);

  return RTC_ERROR;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

static int print_ucd_entire(const struct hecmwST_local_mesh *global_mesh,
                            const char *node_flag, const char *elem_flag,
                            const char *ofname) {
  struct hecmwST_result_data *result_data;

  result_data = (struct hecmwST_result_data *)HECMW_malloc(
      sizeof(struct hecmwST_result_data));
  if (result_data == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  } else {
    init_struct_result_data(result_data);
  }

  if (print_ucd_entire_set_node_data(global_mesh, result_data, node_flag)) {
    goto error;
  }

  if (print_ucd_entire_set_elem_data(global_mesh, result_data, elem_flag)) {
    goto error;
  }

  if (HECMW_ucd_legacy_print(global_mesh, result_data, ofname)) {
    goto error;
  }

  free_struct_result_data(result_data);

  return RTC_NORMAL;

error:
  free_struct_result_data(result_data);

  return RTC_ERROR;
}

static int init_partition(struct hecmwST_local_mesh *global_mesh,
                          struct hecmw_part_cont_data *cont_data) {
  HECMW_log(HECMW_LOG_DEBUG, "Starting initialization for partitioner...");

  /* global_mesh->n_subdomain */
  global_mesh->n_subdomain = cont_data->n_domain;

  /* global_mesh->hecmw_flag_parttype */
  switch (cont_data->type) {
    case HECMW_PART_TYPE_NODE_BASED: /* for node-based partitioning */
      global_mesh->hecmw_flag_parttype = HECMW_FLAG_PARTTYPE_NODEBASED;
      break;

    case HECMW_PART_TYPE_ELEMENT_BASED: /* for element-based partitioning */
      global_mesh->hecmw_flag_parttype = HECMW_FLAG_PARTTYPE_ELEMBASED;
      break;

    default:
      HECMW_set_error(HECMW_PART_E_INVALID_PTYPE, "%d", cont_data->type);
      goto error;
  }

  /* global_mesh->hecmw_flag_partdepth */
  global_mesh->hecmw_flag_partdepth = cont_data->depth;

  /* global_mesh->hecmw_flag_partcontact */
  if (global_mesh->contact_pair->n_pair > 0) {
    switch (cont_data->contact) {
      case HECMW_PART_CONTACT_AGGREGATE:
        global_mesh->hecmw_flag_partcontact = HECMW_FLAG_PARTCONTACT_AGGREGATE;
        break;

      case HECMW_PART_CONTACT_DISTRIBUTE:
        global_mesh->hecmw_flag_partcontact = HECMW_FLAG_PARTCONTACT_DISTRIBUTE;
        break;

      case HECMW_PART_CONTACT_SIMPLE:
        global_mesh->hecmw_flag_partcontact = HECMW_FLAG_PARTCONTACT_SIMPLE;
        break;

      case HECMW_PART_CONTACT_DEFAULT:
      default:
        cont_data->contact = HECMW_PART_CONTACT_SIMPLE;
        global_mesh->hecmw_flag_partcontact = HECMW_FLAG_PARTCONTACT_SIMPLE;
        break;
    }
  }

  HECMW_log(HECMW_LOG_DEBUG, "Initialization for partitioner done");

  return RTC_NORMAL;

error:
  return RTC_ERROR;
  ;
}

/*==================================================================================================

  main function

==================================================================================================*/

extern struct hecmwST_local_mesh *HECMW_partition_inner(
    struct hecmwST_local_mesh *global_mesh,
    struct hecmw_part_cont_data *cont_data) {
  struct hecmwST_local_mesh *local_mesh = NULL;
  struct hecmw_ctrl_meshfiles *ofheader = NULL;
  char *node_flag                       = NULL;
  char *elem_flag                       = NULL;
  char *node_flag_neighbor              = NULL;
  char *elem_flag_neighbor              = NULL;
  int *node_global2local                = NULL;
  int *elem_global2local                = NULL;
  char ofname[HECMW_FILENAME_LEN + 1];
  int *num_elem, *num_node, *num_ielem, *num_inode, *num_nbpe;
  int *sum_elem, *sum_node, *sum_ielem, *sum_inode, *sum_nbpe;
  int current_domain, nrank, iS, iE;
  int rtc;
  int i;
  int error_in_ompsection = 0;

  if (global_mesh == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'global_mesh\' is NULL");
    goto error;
  }
  if (cont_data == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'cont_data\' is NULL");
    goto error;
  }

  rtc = init_partition(global_mesh, cont_data);
  if (rtc != RTC_NORMAL) goto error;

  rtc = HECMW_part_init_log(global_mesh->n_subdomain);
  if (rtc != RTC_NORMAL) goto error;

  if (global_mesh->my_rank == 0) {
    rtc = HECMW_part_set_log_part_type(cont_data->type);
    if (rtc != RTC_NORMAL) goto error;
    rtc = HECMW_part_set_log_part_method(cont_data->method);
    if (rtc != RTC_NORMAL) goto error;
    rtc = HECMW_part_set_log_part_depth(cont_data->depth);
    if (rtc != RTC_NORMAL) goto error;
    rtc = HECMW_part_set_log_part_contact(cont_data->contact);
    if (rtc != RTC_NORMAL) goto error;

    rtc = HECMW_part_set_log_n_node_g(global_mesh->n_node);
    if (rtc != RTC_NORMAL) goto error;
    rtc = HECMW_part_set_log_n_elem_g(global_mesh->n_elem);
    if (rtc != RTC_NORMAL) goto error;
  }

  if (global_mesh->n_subdomain == 1) {
    current_domain = 0;

    if (global_mesh->my_rank == 0) {
      HECMW_log(HECMW_LOG_INFO, "Creating local mesh for domain #%d ...",
                current_domain);

      ofheader = HECMW_ctrl_get_meshfiles_header_sub(
          "part_out", global_mesh->n_subdomain, current_domain);
      if (ofheader == NULL) {
        HECMW_log(HECMW_LOG_ERROR, "not set output file header");
        error_in_ompsection = 1;
        goto error;
      }
      if (ofheader->n_mesh == 0) {
        HECMW_log(HECMW_LOG_ERROR, "output file name is not set");
        error_in_ompsection = 1;
        goto error;
      }

      get_dist_file_name(ofheader->meshfiles[0].filename, current_domain,
                         ofname);
      HECMW_assert(ofname != NULL);

      HECMW_log(HECMW_LOG_DEBUG,
                "Starting writing local mesh for domain #%d...",
                current_domain);

      rtc = HECMW_put_dist_mesh(global_mesh, ofname);
      if (rtc != 0) {
        HECMW_log(HECMW_LOG_ERROR, "Failed to write local mesh for domain #%d",
                  current_domain);
        goto error;
      }

      HECMW_log(HECMW_LOG_DEBUG, "Writing local mesh for domain #%d done",
                current_domain);

      rtc = HECMW_part_set_log_n_elem(0, global_mesh->n_elem);
      if (rtc != 0) goto error;
      rtc = HECMW_part_set_log_n_node(0, global_mesh->n_node);
      if (rtc != 0) goto error;
      rtc = HECMW_part_set_log_ne_internal(0, global_mesh->ne_internal);
      if (rtc != 0) goto error;
      rtc = HECMW_part_set_log_nn_internal(0, global_mesh->nn_internal);
      if (rtc != 0) goto error;

      rtc = HECMW_part_print_log();
      if (rtc) goto error;
    }
    HECMW_part_finalize_log();

    return global_mesh;
  }

  num_elem = (int *)HECMW_calloc(global_mesh->n_subdomain, sizeof(int));
  if (num_elem == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  num_node = (int *)HECMW_calloc(global_mesh->n_subdomain, sizeof(int));
  if (num_node == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  num_ielem = (int *)HECMW_calloc(global_mesh->n_subdomain, sizeof(int));
  if (num_ielem == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  num_inode = (int *)HECMW_calloc(global_mesh->n_subdomain, sizeof(int));
  if (num_inode == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  num_nbpe = (int *)HECMW_calloc(global_mesh->n_subdomain, sizeof(int));
  if (num_nbpe == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  sum_elem = (int *)HECMW_calloc(global_mesh->n_subdomain, sizeof(int));
  if (sum_elem == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  sum_node = (int *)HECMW_calloc(global_mesh->n_subdomain, sizeof(int));
  if (sum_node == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  sum_ielem = (int *)HECMW_calloc(global_mesh->n_subdomain, sizeof(int));
  if (sum_ielem == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  sum_inode = (int *)HECMW_calloc(global_mesh->n_subdomain, sizeof(int));
  if (sum_inode == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }
  sum_nbpe = (int *)HECMW_calloc(global_mesh->n_subdomain, sizeof(int));
  if (sum_nbpe == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  rtc = wnumbering(global_mesh, cont_data);
  if (rtc != RTC_NORMAL) goto error;

  if (cont_data->is_print_part == 1) {
    if (global_mesh->my_rank == 0) {
      print_part(global_mesh, cont_data->part_file_name);
    }
  }

  /*K. Inagaki */
  rtc = spdup_makelist_main(global_mesh);
  if (rtc != RTC_NORMAL) goto error;

#ifdef _OPENMP
#pragma omp parallel default(none), \
    private(node_flag, elem_flag, local_mesh, nrank, iS, iE, i,         \
            current_domain, rtc, ofheader, ofname),                     \
    private(node_global2local, elem_global2local,                       \
            node_flag_neighbor, elem_flag_neighbor),                    \
    shared(global_mesh, cont_data, num_elem, num_node,                  \
           num_ielem, num_inode, num_nbpe, error_in_ompsection)
  {
#endif /* _OPENMP */

    node_flag = (char *)HECMW_calloc(global_mesh->n_node, sizeof(char));
    if (node_flag == NULL) {
      HECMW_set_error(errno, "");
      error_in_ompsection = 1;
      goto error_omp;
    }
    elem_flag = (char *)HECMW_calloc(global_mesh->n_elem, sizeof(char));
    if (elem_flag == NULL) {
      HECMW_set_error(errno, "");
      error_in_ompsection = 1;
      goto error_omp;
    }

    /*K. Inagaki */
    node_global2local = (int *)HECMW_calloc(global_mesh->n_node, sizeof(int));
    if (node_global2local == NULL) {
      HECMW_set_error(errno, "");
      error_in_ompsection = 1;
      goto error_omp;
    }
    elem_global2local = (int *)HECMW_calloc(global_mesh->n_elem, sizeof(int));
    if (elem_global2local == NULL) {
      HECMW_set_error(errno, "");
      error_in_ompsection = 1;
      goto error_omp;
    }
    node_flag_neighbor =
        (char *)HECMW_malloc(sizeof(char) * global_mesh->n_node);
    if (node_flag_neighbor == NULL) {
      HECMW_set_error(errno, "");
      error_in_ompsection = 1;
      goto error_omp;
    }
    elem_flag_neighbor =
        (char *)HECMW_malloc(sizeof(char) * global_mesh->n_elem);
    if (elem_flag_neighbor == NULL) {
      HECMW_set_error(errno, "");
      error_in_ompsection = 1;
      goto error_omp;
    }
    memset(node_flag_neighbor, 0, sizeof(char) * global_mesh->n_node);
    memset(elem_flag_neighbor, 0, sizeof(char) * global_mesh->n_elem);

    local_mesh = HECMW_dist_alloc();
    if (local_mesh == NULL) {
      error_in_ompsection = 1;
      goto error_omp;
    }

    nrank = global_mesh->n_subdomain / HECMW_comm_get_size();
    iS    = HECMW_comm_get_rank() * nrank;
    iE    = iS + nrank;
    if (HECMW_comm_get_rank() == HECMW_comm_get_size() - 1)
      iE = global_mesh->n_subdomain;

#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1), reduction(+ : error_in_ompsection)
#endif
    for (i = iS; i < iE; i++) {
      if (error_in_ompsection) continue;

      current_domain = i;

      HECMW_log(HECMW_LOG_INFO, "Creating local mesh for domain #%d ...",
                current_domain);

      rtc = create_neighbor_info(global_mesh, local_mesh, node_flag, elem_flag,
                                 current_domain);
      if (rtc != RTC_NORMAL) {
        error_in_ompsection = 1;
        continue;
      }

      if (global_mesh->n_subdomain > 1) {
        rtc = create_comm_info(global_mesh, local_mesh, node_flag, elem_flag,
                               node_flag_neighbor, elem_flag_neighbor,
                               current_domain);
        if (rtc != RTC_NORMAL) {
          error_in_ompsection = 1;
          continue;
        }
      }

      rtc = const_local_data(global_mesh, local_mesh, cont_data, node_flag,
                             elem_flag, node_global2local, elem_global2local,
                             current_domain);
      if (rtc != RTC_NORMAL) {
        error_in_ompsection = 1;
        continue;
      }

      num_elem[i]  = local_mesh->n_elem;
      num_node[i]  = local_mesh->n_node;
      num_ielem[i] = local_mesh->ne_internal;
      num_inode[i] = local_mesh->nn_internal;
      num_nbpe[i]  = local_mesh->n_neighbor_pe;

      ofheader = HECMW_ctrl_get_meshfiles_header_sub(
          "part_out", global_mesh->n_subdomain, current_domain);
      if (ofheader == NULL) {
        HECMW_log(HECMW_LOG_ERROR, "not set output file header");
        error_in_ompsection = 1;
        continue;
      }
      if (ofheader->n_mesh == 0) {
        HECMW_log(HECMW_LOG_ERROR, "output file name is not set");
        error_in_ompsection = 1;
        continue;
      }

      get_dist_file_name(ofheader->meshfiles[0].filename, current_domain,
                         ofname);
      HECMW_assert(ofname != NULL);

      HECMW_log(HECMW_LOG_DEBUG,
                "Starting writing local mesh for domain #%d...",
                current_domain);

      rtc = HECMW_put_dist_mesh(local_mesh, ofname);
      if (rtc != 0) {
        HECMW_log(HECMW_LOG_ERROR, "Failed to write local mesh for domain #%d",
                  current_domain);
        error_in_ompsection = 1;
      } else {
        HECMW_log(HECMW_LOG_DEBUG, "Writing local mesh for domain #%d done",
                  current_domain);
      }

      clean_struct_local_mesh(local_mesh);

      HECMW_ctrl_free_meshfiles(ofheader);
      ofheader = NULL;

      if (is_spdup_available(global_mesh)) {
        /*K. Inagaki */
        spdup_clear_IEB(node_flag, elem_flag, current_domain);
      } else {
        int j;
        for (j = 0; j < global_mesh->n_node; j++) {
          CLEAR_IEB(node_flag[j]);
        }
        for (j = 0; j < global_mesh->n_elem; j++) {
          CLEAR_IEB(elem_flag[j]);
        }
      }
    }
#ifdef _OPENMP
    if (error_in_ompsection) goto error_omp;

#pragma omp single
#endif
    if (cont_data->is_print_ucd == 1) {
      if (global_mesh->my_rank == 0) {
        print_ucd_entire(global_mesh, node_flag, elem_flag,
                         cont_data->ucd_file_name);
      }
    }

  error_omp:
    HECMW_dist_free(local_mesh);
    HECMW_free(node_flag);
    HECMW_free(elem_flag);
    /*K. Inagaki */
    HECMW_free(node_global2local);
    HECMW_free(elem_global2local);
    HECMW_free(node_flag_neighbor);
    HECMW_free(elem_flag_neighbor);

#ifdef _OPENMP
  } /* omp end parallel */
  if (error_in_ompsection) goto error;
#endif

  rtc = HECMW_Allreduce(num_elem, sum_elem, global_mesh->n_subdomain, HECMW_INT,
                        HECMW_SUM, HECMW_comm_get_comm());
  if (rtc != 0) goto error;
  rtc = HECMW_Allreduce(num_node, sum_node, global_mesh->n_subdomain, HECMW_INT,
                        HECMW_SUM, HECMW_comm_get_comm());
  if (rtc != 0) goto error;
  rtc = HECMW_Allreduce(num_ielem, sum_ielem, global_mesh->n_subdomain,
                        HECMW_INT, HECMW_SUM, HECMW_comm_get_comm());
  if (rtc != 0) goto error;
  rtc = HECMW_Allreduce(num_inode, sum_inode, global_mesh->n_subdomain,
                        HECMW_INT, HECMW_SUM, HECMW_comm_get_comm());
  if (rtc != 0) goto error;
  rtc = HECMW_Allreduce(num_nbpe, sum_nbpe, global_mesh->n_subdomain,
                        HECMW_INT, HECMW_SUM, HECMW_comm_get_comm());
  if (rtc != 0) goto error;

  if (global_mesh->my_rank == 0) {
    for (i = 0; i < global_mesh->n_subdomain; i++) {
      rtc = HECMW_part_set_log_n_elem(i, sum_elem[i]);
      if (rtc != 0) goto error;
      rtc = HECMW_part_set_log_n_node(i, sum_node[i]);
      if (rtc != 0) goto error;
      rtc = HECMW_part_set_log_ne_internal(i, sum_ielem[i]);
      if (rtc != 0) goto error;
      rtc = HECMW_part_set_log_nn_internal(i, sum_inode[i]);
      if (rtc != 0) goto error;
      rtc = HECMW_part_set_log_n_neighbor_pe(i, sum_nbpe[i]);
      if (rtc != 0) goto error;
    }
    rtc = HECMW_part_print_log();
    if (rtc) goto error;
  }
  HECMW_part_finalize_log();

  HECMW_free(num_elem);
  HECMW_free(num_node);
  HECMW_free(num_ielem);
  HECMW_free(num_inode);
  HECMW_free(num_nbpe);
  HECMW_free(sum_elem);
  HECMW_free(sum_node);
  HECMW_free(sum_ielem);
  HECMW_free(sum_inode);
  HECMW_free(sum_nbpe);

  /*K. Inagaki */
  spdup_freelist(global_mesh);

  return global_mesh;

error:
  HECMW_free(node_flag);
  HECMW_free(elem_flag);
  HECMW_free(num_elem);
  HECMW_free(num_node);
  HECMW_free(num_ielem);
  HECMW_free(num_inode);
  HECMW_free(num_nbpe);
  HECMW_free(sum_elem);
  HECMW_free(sum_node);
  HECMW_free(sum_ielem);
  HECMW_free(sum_inode);
  HECMW_free(sum_nbpe);
  HECMW_dist_free(local_mesh);
  if (ofheader) {
    HECMW_ctrl_free_meshfiles(ofheader);
  }
  HECMW_part_finalize_log();

  return NULL;
}

extern struct hecmwST_local_mesh *HECMW_partition(
    struct hecmwST_local_mesh *global_mesh) {
  struct hecmwST_local_mesh *local_mesh;
  struct hecmw_part_cont_data *cont_data;

  HECMW_log(HECMW_LOG_INFO, "Starting domain decomposition...\n");

  if (global_mesh == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "\'global_mesh\' is NULL");
    goto error;
  }

  cont_data = HECMW_part_get_control(global_mesh);
  if (cont_data == NULL) goto error;

  local_mesh = HECMW_partition_inner(global_mesh, cont_data);
  if (local_mesh == NULL) goto error;

  HECMW_part_free_control(cont_data);

  HECMW_log(HECMW_LOG_INFO, "Domain decomposition done\n");

  return local_mesh;

error:
  return NULL;
}
