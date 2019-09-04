/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "hecmw_struct.h"
#include "hecmw_result.h"
#include "hecmw_util.h"
#include "hecmw_dist_copy_f2c.h"
#include "hecmw_result_copy_f2c.h"
#include "hecmw_dist_free.h"
#include "hecmw_visualizer.h"

static struct hecmwST_local_mesh *mesh;
static struct hecmwST_result_data *result;

/*----------------------------------------------------------------------------*/
static int alloc_local_mesh(void) {
  mesh = HECMW_malloc(sizeof(*mesh));
  if (mesh == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  mesh->section = HECMW_malloc(sizeof(*mesh->section));
  if (mesh->section == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  mesh->material = HECMW_malloc(sizeof(*mesh->material));
  if (mesh->material == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  mesh->mpc = HECMW_malloc(sizeof(*mesh->mpc));
  if (mesh->mpc == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  mesh->amp = HECMW_malloc(sizeof(*mesh->amp));
  if (mesh->amp == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  mesh->node_group = HECMW_malloc(sizeof(*mesh->node_group));
  if (mesh->node_group == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  mesh->elem_group = HECMW_malloc(sizeof(*mesh->elem_group));
  if (mesh->elem_group == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  mesh->surf_group = HECMW_malloc(sizeof(*mesh->surf_group));
  if (mesh->surf_group == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  mesh->contact_pair = HECMW_malloc(sizeof(*mesh->contact_pair));
  if (mesh->contact_pair == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
static int alloc_result(void) {
  result = HECMW_malloc(sizeof(*result));
  if (result == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }
  return 0;
}

/*----------------------------------------------------------------------------*/
void hecmw_visualize_if(int *step, int *err) {
  *err = 1;
  if (HECMW_visualize(mesh, result, *step)) return;
  *err = 0;
}

void hecmw_visualize_if_(int *step, int *err) {
  hecmw_visualize_if(step, err);
}

void hecmw_visualize_if__(int *step, int *err) {
  hecmw_visualize_if(step, err);
}

void HECMW_VISUALIZE_IF(int *step, int *err) {
  hecmw_visualize_if(step, err);
}

/*----------------------------------------------------------------------------*/
void hecmw_visualize_init_if(int *nnode, int *nelem, int *err) {
  *err = 1;

  if (alloc_local_mesh()) return;

  if (alloc_result()) return;

  if (HECMW_dist_copy_f2c_init(mesh)) return;

  if (HECMW_result_copy_f2c_init(result, *nnode, *nelem)) return;

  *err = 0;
}

void hecmw_visualize_init_if_(int *nnode, int *nelem, int *err) {
  hecmw_visualize_init_if(nnode, nelem, err);
}

void hecmw_visualize_init_if__(int *nnode, int *nelem, int *err) {
  hecmw_visualize_init_if(nnode, nelem, err);
}

void HECMW_VISUALIZE_INIT_IF(int *nnode, int *nelem, int *err) {
  hecmw_visualize_init_if(nnode, nelem, err);
}

/*----------------------------------------------------------------------------*/
void hecmw_visualize_finalize_if(int *err) {
  *err = 0;

  if (HECMW_dist_copy_f2c_finalize()) return;

  if (HECMW_result_copy_f2c_finalize()) return;

  HECMW_dist_free(mesh);
  mesh = NULL;

  HECMW_result_free(result);
  result = NULL;

  *err = 0;
}

void hecmw_visualize_finalize_if_(int *err) {
  hecmw_visualize_finalize_if(err);
}

void hecmw_visualize_finalize_if__(int *err) {
  hecmw_visualize_finalize_if(err);
}

void HECMW_VISUALIZE_FINALIZE_IF(int *err) { hecmw_visualize_finalize_if(err); }

/*----------------------------------------------------------------------------*/
void hecmw_init_for_visual_if(int *err) {
  *err = 1;
  if (HECMW_visualize_init()) return;
  *err = 0;
}

void hecmw_init_for_visual_if_(int *err) { hecmw_init_for_visual_if(err); }

void hecmw_init_for_visual_if__(int *err) { hecmw_init_for_visual_if(err); }

void HECMW_INIT_FOR_VISUAL_IF(int *err) { hecmw_init_for_visual_if(err); }

/*----------------------------------------------------------------------------*/
void hecmw_finalize_for_visual_if(int *err) {
  *err = 1;
  if (HECMW_visualize_finalize()) return;
  *err = 0;
}

void hecmw_finalize_for_visual_if_(int *err) {
  hecmw_finalize_for_visual_if(err);
}

void hecmw_finalize_for_visual_if__(int *err) {
  hecmw_finalize_for_visual_if(err);
}

void HECMW_FINALIZE_FOR_VISUAL_IF(int *err) {
  hecmw_finalize_for_visual_if(err);
}
