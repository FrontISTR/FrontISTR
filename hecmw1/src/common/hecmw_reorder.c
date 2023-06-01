/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "hecmw_struct.h"
#include "hecmw_util.h"
#include "hecmw_common_define.h"
#include "hecmw_etype.h"
#include "hecmw_reorder.h"

#define MASK_BIT(map, bit) ((map) |= (bit))
#define EVAL_BIT(map, bit) ((map) & (bit))
#define INV_BIT(map, bit) ((map) ^= (bit))
#define CLEAR_BIT(map, bit) \
  ((map) |= (bit));         \
  ((map) ^= (bit))

#define BIT_DOF_TWO 1
#define BIT_DOF_THREE 2
#define BIT_DOF_SIX 4
#define BIT_DOF_FOUR 8
#define BIT_DOF_ALL (BIT_DOF_TWO | BIT_DOF_THREE | BIT_DOF_SIX | BIT_DOF_FOUR)

#define HECMW_COMMON_EQUATION_BLOCK_NAME "EQUATION_BLOCK"
#define MY_RANK 1
#define NEIGHBOR_RANK 2
#define BOTH_RANK 3
#define MPC_BLOCK 4
#define CANDIDATE 8
#define ALL_BIT 255

#ifdef DEBUG
#define dw_node_flag(ptr, nmemb) dw_node_flag_(ptr, nmemb, __FILE__, __LINE__)
#else
#define dw_node_flag(ptr, nmemb) ((void)0)
#endif

struct equation_block {
  int n_eqn_block;
  int n_mpc_block;
  int *eqn_block_index;
};

/*----------------------------------------------------------------------------*/
#if 0
static int
get_eqn_block_idx( struct hecmwST_local_mesh *local_mesh )
{
  int i;

  for( i=0; i<local_mesh->node_group->n_grp; i++ )
    if( !strcmp( local_mesh->node_group->grp_name[i], HECMW_COMMON_EQUATION_BLOCK_NAME ) ) return i;

  return -1;
}
#endif

/*                                                                            */
/*  convert node id from old to new                                           */
/*                                                                            */

/*----------------------------------------------------------------------------*/
/*  nodal coordinates < node >                                                */
/*----------------------------------------------------------------------------*/
static int old2new_node(struct hecmwST_local_mesh *local_mesh,
                        int *node_new2old) {
  double *new, *old;
  int i;

  new = (double *)HECMW_malloc(sizeof(double) * local_mesh->n_node * 3);
  if (new == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  old = local_mesh->node;
  HECMW_assert(old);

  for (i = 0; i < local_mesh->n_node; i++) {
    new[3 * i]     = old[3 * (node_new2old[i] - 1)];
    new[3 * i + 1] = old[3 * (node_new2old[i] - 1) + 1];
    new[3 * i + 2] = old[3 * (node_new2old[i] - 1) + 2];
  }

  local_mesh->node = new;

  HECMW_free(old);

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  local node id & belonging domain of node < node_ID >                      */
/*----------------------------------------------------------------------------*/
static int old2new_node_ID(struct hecmwST_local_mesh *local_mesh,
                           int *node_new2old) {
  int *new, *old;
  int i;

  new = (int *)HECMW_malloc(sizeof(int) * local_mesh->n_node * 2);
  if (new == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  old = local_mesh->node_ID;
  HECMW_assert(old);

  for (i = 0; i < local_mesh->n_node; i++) {
    new[2 * i]     = old[2 * (node_new2old[i] - 1)];
    new[2 * i + 1] = old[2 * (node_new2old[i] - 1) + 1];
  }

  local_mesh->node_ID = new;

  HECMW_free(old);

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  global node id                                                            */
/*----------------------------------------------------------------------------*/
static int old2new_global_node_ID(struct hecmwST_local_mesh *local_mesh,
                                  int *node_new2old) {
  int *new, *old;
  int i;

  new = (int *)HECMW_malloc(sizeof(int) * local_mesh->n_node);
  if (new == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  old = local_mesh->global_node_ID;
  HECMW_assert(old);

  for (i = 0; i < local_mesh->n_node; i++) {
    new[i] = old[node_new2old[i] - 1];
  }

  local_mesh->global_node_ID = new;

  HECMW_free(old);

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  initial conditions of node                                                */
/*----------------------------------------------------------------------------*/
static int old2new_node_init_val(struct hecmwST_local_mesh *local_mesh,
                                 int *node_new2old) {
  int *new_index, *old_index;
  double *new_item, *old_item;
  int old_id;
  int counter;
  int i, j;

  new_index = (int *)HECMW_calloc(local_mesh->n_node + 1, sizeof(int));
  if (new_index == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }
  new_item = (double *)HECMW_malloc(
      sizeof(double) * local_mesh->node_init_val_index[local_mesh->n_node]);
  if (new_item == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  old_index = local_mesh->node_init_val_index;
  old_item  = local_mesh->node_init_val_item;
  HECMW_assert(old_index);
  HECMW_assert(old_item);

  for (counter = 0, i = 0; i < local_mesh->n_node; i++) {
    old_id = node_new2old[i];

    for (j = old_index[old_id - 1]; j < old_index[old_id]; j++) {
      new_item[counter++] = old_item[j];
    }
    new_index[i + 1] = counter;
  }

  local_mesh->node_init_val_index = new_index;
  local_mesh->node_init_val_item  = new_item;

  HECMW_free(old_item);
  HECMW_free(old_index);

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  component nodes of element < elem_node_item >                             */
/*----------------------------------------------------------------------------*/
static int old2new_elem_node_item(struct hecmwST_local_mesh *local_mesh,
                                  int *node_old2new) {
  int i;

  for (i = 0; i < local_mesh->elem_node_index[local_mesh->n_elem]; i++) {
    local_mesh->elem_node_item[i] =
        node_old2new[local_mesh->elem_node_item[i] - 1];
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  MPC group < local_mesh->mpc >                                             */
/*----------------------------------------------------------------------------*/
static int old2new_mpc_item(struct hecmwST_local_mesh *local_mesh,
                            int *node_old2new) {
  int i;

  if (!local_mesh->mpc->n_mpc) return 0;

  for (i = 0; i < local_mesh->mpc->mpc_index[local_mesh->mpc->n_mpc]; i++) {
    local_mesh->mpc->mpc_item[i] =
        node_old2new[local_mesh->mpc->mpc_item[i] - 1];
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  node group < local_mesh->node_group->grp_item >                           */
/*----------------------------------------------------------------------------*/
static int old2new_node_grp_item(struct hecmwST_local_mesh *local_mesh,
                                 int *node_old2new) {
  int i, j;

  if (!local_mesh->node_group->n_grp) return 0;

  for (i = 0; i < local_mesh->node_group->n_grp; i++) {
    if (strcmp(local_mesh->node_group->grp_name[i],
               HECMW_COMMON_EQUATION_BLOCK_NAME)) {
      for (j = local_mesh->node_group->grp_index[i];
           j < local_mesh->node_group->grp_index[i + 1]; j++) {
        local_mesh->node_group->grp_item[j] =
            node_old2new[local_mesh->node_group->grp_item[j] - 1];
      }
    }
  }

  return 0;
}

/*============================================================================*/
/*  convert node id from old to new                                           */
/*============================================================================*/
static int old2new_node_info(struct hecmwST_local_mesh *local_mesh,
                             int *node_new2old, int *node_old2new) {
  /* nodal coordinates */
  if (old2new_node(local_mesh, node_new2old)) return -1;

  /* local node id & belonging domain of node */
  if (old2new_node_ID(local_mesh, node_new2old)) return -1;

  /* global node id */
  if (old2new_global_node_ID(local_mesh, node_new2old)) return -1;

  /* initial conditions of node */
  if (local_mesh->hecmw_flag_initcon) {
    if (old2new_node_init_val(local_mesh, node_new2old)) return -1;
  }

  /* component nodes of element */
  if (old2new_elem_node_item(local_mesh, node_old2new)) return -1;

  /* MPC group */
  if (local_mesh->mpc->n_mpc) {
    if (old2new_mpc_item(local_mesh, node_old2new)) return -1;
  }

  /* node group */
  if (local_mesh->node_group->n_grp) {
    if (old2new_node_grp_item(local_mesh, node_old2new)) return -1;
  }

  return 0;
}

/*                                                                            */
/*  convert element id from old to new                                        */
/*                                                                            */

/*----------------------------------------------------------------------------*/
/*  finite element type < elem_type >                                         */
/*----------------------------------------------------------------------------*/
static int old2new_elem_type(struct hecmwST_local_mesh *local_mesh,
                             int *elem_new2old) {
  int *new, *old;
  int i;

  new = HECMW_malloc(sizeof(int) * local_mesh->n_elem);
  if (new == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  old = local_mesh->elem_type;
  HECMW_assert(old);

  for (i = 0; i < local_mesh->n_elem; i++) {
    new[i] = local_mesh->elem_type[elem_new2old[i] - 1];
  }

  local_mesh->elem_type = new;

  HECMW_free(old);

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  component nodes of element < elem_node_index, elem_node_item >            */
/*----------------------------------------------------------------------------*/
static int old2new_elem_node(struct hecmwST_local_mesh *local_mesh,
                             int *elem_new2old) {
  int *new_index, *old_index;
  int *new_item, *old_item;
  int old_id;
  int counter;
  int i, j;

  new_index = HECMW_calloc(local_mesh->n_elem + 1, sizeof(int));
  if (new_index == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }
  new_item = HECMW_malloc(sizeof(int) *
                          local_mesh->elem_node_index[local_mesh->n_elem]);
  if (new_item == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  old_index = local_mesh->elem_node_index;
  old_item  = local_mesh->elem_node_item;
  HECMW_assert(old_index);
  HECMW_assert(old_item);

  for (counter = 0, i = 0; i < local_mesh->n_elem; i++) {
    old_id = elem_new2old[i];

    for (j = old_index[old_id - 1]; j < old_index[old_id]; j++) {
      new_item[counter++] = old_item[j];
    }
    new_index[i + 1] = counter;
  }
  HECMW_assert(counter > local_mesh->n_elem);

  local_mesh->elem_node_index = new_index;
  local_mesh->elem_node_item  = new_item;

  HECMW_free(old_item);
  HECMW_free(old_index);

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  local element ID & belonging domain of element < elem_ID >                */
/*----------------------------------------------------------------------------*/
static int old2new_elem_ID(struct hecmwST_local_mesh *local_mesh,
                           int *elem_new2old) {
  int *new, *old;
  int i;

  new = (int *)HECMW_malloc(sizeof(int) * local_mesh->n_elem * 2);
  if (new == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  old = local_mesh->elem_ID;
  HECMW_assert(old);

  for (i = 0; i < local_mesh->n_elem; i++) {
    new[2 * i]     = old[2 * (elem_new2old[i] - 1)];
    new[2 * i + 1] = old[2 * (elem_new2old[i] - 1) + 1];
  }

  local_mesh->elem_ID = new;

  HECMW_free(old);

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  global element id < global_elem_ID >                                      */
/*----------------------------------------------------------------------------*/
static int old2new_global_elem_ID(struct hecmwST_local_mesh *local_mesh,
                                  int *elem_new2old) {
  int *new, *old;
  int i;

  new = (int *)HECMW_malloc(sizeof(int) * local_mesh->n_elem);
  if (new == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  old = local_mesh->global_elem_ID;
  HECMW_assert(old);

  for (i = 0; i < local_mesh->n_elem; i++) {
    new[i] = old[elem_new2old[i] - 1];
  }

  local_mesh->global_elem_ID = new;

  HECMW_free(old);

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  section id < section_ID >                                                 */
/*----------------------------------------------------------------------------*/
static int old2new_section_ID(struct hecmwST_local_mesh *local_mesh,
                              int *elem_new2old) {
  int *new, *old;
  int i;

  new = (int *)HECMW_malloc(sizeof(int) * local_mesh->n_elem);
  if (new == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  old = local_mesh->section_ID;
  HECMW_assert(old);

  for (i = 0; i < local_mesh->n_elem; i++) {
    new[i] = local_mesh->section_ID[elem_new2old[i] - 1];
  }

  local_mesh->section_ID = new;

  HECMW_free(old);

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  material id < elem_mat_ID_index, elem_mat_ID_item >                       */
/*----------------------------------------------------------------------------*/
static int old2new_mat_ID(struct hecmwST_local_mesh *local_mesh,
                          int *elem_new2old) {
  int *new_index, *old_index;
  int *new_item, *old_item;
  int old_id;
  int counter;
  int i, j;

  new_index = (int *)HECMW_calloc(local_mesh->n_elem + 1, sizeof(int));
  if (new_index == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }
  new_item = (int *)HECMW_malloc(
      sizeof(int) * local_mesh->elem_mat_ID_index[local_mesh->n_elem]);
  if (new_item == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  old_index = local_mesh->elem_mat_ID_index;
  old_item  = local_mesh->elem_mat_ID_item;
  HECMW_assert(old_index);
  HECMW_assert(old_item);

  for (counter = 0, i = 0; i < local_mesh->n_elem; i++) {
    old_id = elem_new2old[i];

    for (j = old_index[old_id - 1]; j < old_index[old_id]; j++) {
      new_item[counter++] = old_item[j];
    }
    new_index[i + 1] = counter;
  }
  HECMW_assert(counter >= local_mesh->n_elem);

  local_mesh->elem_mat_ID_index = new_index;
  local_mesh->elem_mat_ID_item  = new_item;

  HECMW_free(old_item);
  HECMW_free(old_index);

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  list of internal element < elem_internal_list >                           */
/*----------------------------------------------------------------------------*/
static int old2new_elem_internal_list(struct hecmwST_local_mesh *local_mesh,
                                      int *elem_old2new) {
  int i;

  for (i = 0; i < local_mesh->ne_internal; i++) {
    local_mesh->elem_internal_list[i] =
        elem_old2new[local_mesh->elem_internal_list[i] - 1];
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  element group < elem_group->grp_item >                                    */
/*----------------------------------------------------------------------------*/
static int old2new_elem_grp_item(struct hecmwST_local_mesh *local_mesh,
                                 int *elem_old2new) {
  int i;

  for (i = 0;
       i < local_mesh->elem_group->grp_index[local_mesh->elem_group->n_grp];
       i++) {
    local_mesh->elem_group->grp_item[i] =
        elem_old2new[local_mesh->elem_group->grp_item[i] - 1];
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  surface group < surf_group->grp_item >                                    */
/*----------------------------------------------------------------------------*/
static int old2new_surf_grp_item(struct hecmwST_local_mesh *local_mesh,
                                 int *elem_old2new) {
  int i;

  for (i = 0;
       i < local_mesh->surf_group->grp_index[local_mesh->surf_group->n_grp];
       i++) {
    local_mesh->surf_group->grp_item[2 * i] =
        elem_old2new[local_mesh->surf_group->grp_item[2 * i] - 1];
  }

  return 0;
}

/*============================================================================*/
/*  convert element id from old to new                                        */
/*============================================================================*/
static int old2new_elem_info(struct hecmwST_local_mesh *local_mesh,
                             int *elem_new2old, int *elem_old2new) {
  /* finite element type */
  if (old2new_elem_type(local_mesh, elem_new2old)) return -1;

  /* component nodes of element */
  if (old2new_elem_node(local_mesh, elem_new2old)) return -1;

  /* local element id & belonging domain of element */
  if (old2new_elem_ID(local_mesh, elem_new2old)) return -1;

  /* global element id */
  if (old2new_global_elem_ID(local_mesh, elem_new2old)) return -1;

  /* section id */
  if (old2new_section_ID(local_mesh, elem_new2old)) return -1;

  /* material id */
  if (old2new_mat_ID(local_mesh, elem_new2old)) return -1;

  /* list of internal element */
  if (old2new_elem_internal_list(local_mesh, elem_old2new)) return -1;

  /* element group */
  if (local_mesh->elem_group->n_grp) {
    if (old2new_elem_grp_item(local_mesh, elem_old2new)) return -1;
  }

  /* surface group */
  if (local_mesh->surf_group->n_grp) {
    if (old2new_surf_grp_item(local_mesh, elem_old2new)) return -1;
  }

  return 0;
}

/*                                                                            */
/*  reorder member of EQUATION_BLOCK                                          */
/*                                                                            */
#if 0
/*----------------------------------------------------------------------------*/
/*  node group ( only EQUATION_BLOCK ) < node_group->grp_item >               */
/*----------------------------------------------------------------------------*/
static int
old2new_eqn_block( struct hecmwST_local_mesh *local_mesh, int *eqn_block_old2new )
{
  struct hecmwST_node_grp *grp=local_mesh->node_group;
  int n_eqn_block, eqn_block_idx;
  int *new_item;
  int counter;
  int i, js;

  /* index of EQUATION_BLOCK */
  eqn_block_idx = get_eqn_block_idx( local_mesh );
  HECMW_assert( eqn_block_idx >= 0 );

  /* number of EQUATION_BLOCKs */
  n_eqn_block = grp->grp_index[eqn_block_idx+1] - grp->grp_index[eqn_block_idx];
  HECMW_assert( n_eqn_block > 0 );

  /* order EQUATION_BLOCK */
  new_item = (int *)HECMW_calloc( n_eqn_block, sizeof(int) );
  if( new_item == NULL ) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  for( js=0, i=grp->grp_index[eqn_block_idx]; i<grp->grp_index[eqn_block_idx+1]; i++ ) {
    new_item[eqn_block_old2new[i-grp->grp_index[eqn_block_idx]]] = grp->grp_item[i] - js;
    js = grp->grp_item[i];
  }

  for( counter=0, i=0; i<n_eqn_block; i++ ) {
    counter += new_item[i];
    grp->grp_item[grp->grp_index[eqn_block_idx]+i] = counter;
  }

  HECMW_free( new_item );

  return 0;
}

/*============================================================================*/
/*  reorder member of EQUATION_BLOCK                                          */
/*============================================================================*/
static int
old2new_eqn_block_info( struct hecmwST_local_mesh *local_mesh, int *eqn_block_old2new )
{
  /* node group ( only EQUATION_BLOCK ) */
  if(old2new_eqn_block( local_mesh, eqn_block_old2new )) return -1;

  return 0;
}
#endif
/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

/*                                                                            */
/*  reorder element according to finite element type                          */
/*                                                                            */

/*----------------------------------------------------------------------------*/
/*  count elements of each finite element type                                */
/*----------------------------------------------------------------------------*/
static int count_each_elem_type(struct hecmwST_local_mesh *local_mesh,
                                int *counter) {
  int etype;
  int i;

  /* count elements of each finite element type */
  for (i = 0; i < local_mesh->n_elem; i++) {
    etype = HECMW_get_etype_HECMW2UTIL(local_mesh->elem_type[i]);
    counter[etype]++;
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  count element types in mesh                                               */
/*----------------------------------------------------------------------------*/
static int set_n_elem_type(struct hecmwST_local_mesh *local_mesh,
                           int *counter) {
  int i;

  /* count element types */
  local_mesh->n_elem_type = 0;
  for (i = 0; i < HECMW_MESH_ETYPE_MAX + 1; i++) {
    if (counter[i]) {
      (local_mesh->n_elem_type)++;
    }
  }

  /* check data */
  if (local_mesh->n_elem_type <= 0) {
    HECMW_set_error(HECMW_COMMON_E_OUT_OF_RANGE, "");
    return -1;
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  create item and index for finite element type & set conversion table      */
/*  between old and new element id                                            */
/*----------------------------------------------------------------------------*/
static int set_elem_type_index(struct hecmwST_local_mesh *local_mesh,
                               int *counter, int *elem_new2old,
                               int *elem_old2new) {
  int etype;
  int types, elems;
  int i, j;

  /* allocation */
  if (local_mesh->elem_type_index) {
    HECMW_free(local_mesh->elem_type_index);
  }
  local_mesh->elem_type_index =
      (int *)HECMW_calloc(local_mesh->n_elem_type + 1, sizeof(int));
  if (local_mesh->elem_type_index == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  if (local_mesh->elem_type_item) {
    HECMW_free(local_mesh->elem_type_item);
  }
  local_mesh->elem_type_item =
      (int *)HECMW_malloc(sizeof(int) * local_mesh->n_elem_type);
  if (local_mesh->elem_type_item == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  /* create index for element type */
  for (elems = 0, types = 0, i = 0; i < HECMW_MESH_ETYPE_MAX + 1; i++) {
    if (counter[i]) {
      etype = HECMW_get_etype_UTIL2HECMW(i);

      for (j = 0; j < local_mesh->n_elem; j++) {
        if (local_mesh->elem_type[j] == etype) {
          elem_new2old[elems] = j + 1;
          elem_old2new[j]     = elems + 1;
          elems++;
        }
      }
      local_mesh->elem_type_index[types + 1] = elems;
      local_mesh->elem_type_item[types]      = etype;
      types++;
    }
  }

  return 0;
}

/*============================================================================*/
/*  reoreder element according to finite element type                         */
/*============================================================================*/
extern int HECMW_reorder_elem_type(struct hecmwST_local_mesh *local_mesh) {
  int *counter; /* counter of elements */
  int *elem_new2old,
      *elem_old2new; /* conversion table between old and new element id */

  HECMW_assert(local_mesh);

  /* allocation */
  counter = (int *)HECMW_calloc(HECMW_MESH_ETYPE_MAX + 1, sizeof(int));
  if (counter == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  elem_new2old = (int *)HECMW_calloc(local_mesh->n_elem, sizeof(int));
  if (elem_new2old == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  elem_old2new = (int *)HECMW_calloc(local_mesh->n_elem, sizeof(int));
  if (elem_old2new == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  /* count elements of each element type */
  if (count_each_elem_type(local_mesh, counter)) {
    return -1;
  }

  /* count finite element types in data */
  if (set_n_elem_type(local_mesh, counter)) {
    return -1;
  }

  /* create index for finite element type */
  if (set_elem_type_index(local_mesh, counter, elem_new2old, elem_old2new)) {
    return -1;
  }

  /* reorder relevant arrays */
  if (old2new_elem_info(local_mesh, elem_new2old, elem_old2new)) {
    return -1;
  }

  HECMW_free(counter);
  HECMW_free(elem_new2old);
  HECMW_free(elem_old2new);

  return 0;
}

/*                                                                            */
/*  reorder node according to DOF                                             */
/*                                                                            */

/*----------------------------------------------------------------------------*/
/*  mask node according to DOF                                                */
/*----------------------------------------------------------------------------*/
static int mask_node_dof_inner(struct hecmwST_local_mesh *local_mesh,
                               char *node_flag, const int is, const int ie,
                               const int n_comp, const int n_dof) {
  int nidx, node;
  int i, j;

  for (i = is; i < ie; i++) {
    nidx = local_mesh->elem_node_index[i];

    for (j = 0; j < n_comp; j++) {
      node = local_mesh->elem_node_item[nidx + j];
      MASK_BIT(node_flag[node - 1], n_dof);
    }
  }

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int mask_node_dof(struct hecmwST_local_mesh *local_mesh,
                         char *node_flag) {
  int n_comp, n_dof;
  int is, ie;
  int i;
  n_dof = -1;

  for (i = 0; i < local_mesh->n_elem_type; i++) {
    is     = local_mesh->elem_type_index[i];
    ie     = local_mesh->elem_type_index[i + 1];
    n_comp = HECMW_get_max_node(local_mesh->elem_type_item[i]);

    switch (local_mesh->elem_type_item[i]) {
      /* line element */
      case HECMW_ETYPE_ROD1:
      case HECMW_ETYPE_ROD2:
        n_dof = BIT_DOF_TWO;
        break;

      /* surface element */
      case HECMW_ETYPE_TRI1:
      case HECMW_ETYPE_TRI2:
      case HECMW_ETYPE_TRI22:
      case HECMW_ETYPE_QUA1:
      case HECMW_ETYPE_QUA2:
        n_dof = BIT_DOF_TWO;
        break;

      /* solid element */
      case HECMW_ETYPE_TET1:
      case HECMW_ETYPE_TET2:
      case HECMW_ETYPE_TET22:
      case HECMW_ETYPE_PRI1:
      case HECMW_ETYPE_PRI2:
      case HECMW_ETYPE_HEX1:
      case HECMW_ETYPE_HEX2:
      case HECMW_ETYPE_PYR1:
      case HECMW_ETYPE_PYR2:
      case HECMW_ETYPE_ROD31:
        n_dof = BIT_DOF_THREE;
        break;
      case HECMW_ETYPE_TET1_4:
      case HECMW_ETYPE_HEX1_4:
        n_dof = BIT_DOF_FOUR;
        break;

      /* master-slave type element */
      case HECMW_ETYPE_MST1:
      case HECMW_ETYPE_MST2:
      case HECMW_ETYPE_MSQ1:
      case HECMW_ETYPE_MSQ2:
        n_dof = BIT_DOF_THREE;
        break;

      /* interface element */
      case HECMW_ETYPE_JTB1:
      case HECMW_ETYPE_SPGDPT1:
      case HECMW_ETYPE_JTT1:
      case HECMW_ETYPE_JTT2:
      case HECMW_ETYPE_JTQ1:
      case HECMW_ETYPE_JTQ2:
        n_dof = BIT_DOF_THREE;
        break;

      /* beam element */
      case HECMW_ETYPE_BEM1:
      case HECMW_ETYPE_BEM2:
        n_dof = BIT_DOF_SIX;
        break;

      /* surface shell element */
      case HECMW_ETYPE_SHT1:
      case HECMW_ETYPE_SHT2:
      case HECMW_ETYPE_SHQ1:
      case HECMW_ETYPE_SHQ2:
      case HECMW_ETYPE_SHQ3:
        n_dof = BIT_DOF_SIX;
        break;

      /* surface shell element (351 361) */
      case HECMW_ETYPE_SHT6:
      case HECMW_ETYPE_SHQ8:
      case HECMW_ETYPE_BEM3:
        n_dof = BIT_DOF_THREE;
        break;

      /* patch element */
      case HECMW_ETYPE_PTT1:
      case HECMW_ETYPE_PTT2:
      case HECMW_ETYPE_PTQ1:
      case HECMW_ETYPE_PTQ2:
        n_dof = BIT_DOF_THREE;
        break;

      /* link element for MPC */
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
        n_dof = BIT_DOF_TWO;
        break;

      default:
        HECMW_assert(0);
    }

    if (mask_node_dof_inner(local_mesh, node_flag, is, ie, n_comp, n_dof))
      return -1;
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  reorder node according to DOF & create conversion table                   */
/*----------------------------------------------------------------------------*/
static int reorder_node_dof(struct hecmwST_local_mesh *local_mesh,
                            char *node_flag, int *node_new2old,
                            int *node_old2new, char *dof_flag, int *n_dof_tot) {
  int counter = 0;
  int i;

  /* six DOF */
  for (i = 0; i < local_mesh->n_node; i++) {
    if (EVAL_BIT(node_flag[i], BIT_DOF_SIX)) {
      node_old2new[i]       = counter + 1;
      node_new2old[counter] = i + 1;
      counter++;

      (n_dof_tot[HECMW_MESH_DOF_SIX])++;
      MASK_BIT(*dof_flag, BIT_DOF_SIX);
      CLEAR_BIT(node_flag[i], BIT_DOF_ALL);
    }
  }

  /* four DOF */
  for (i = 0; i < local_mesh->n_node; i++) {
    if (EVAL_BIT(node_flag[i], BIT_DOF_FOUR)) {
      node_old2new[i]       = counter + 1;
      node_new2old[counter] = i + 1;
      counter++;

      (n_dof_tot[HECMW_MESH_DOF_FOUR])++;
      MASK_BIT(*dof_flag, BIT_DOF_FOUR);
      CLEAR_BIT(node_flag[i], BIT_DOF_ALL);
    }
  }

  /* three DOF */
  for (i = 0; i < local_mesh->n_node; i++) {
    if (EVAL_BIT(node_flag[i], BIT_DOF_THREE)) {
      node_old2new[i]       = counter + 1;
      node_new2old[counter] = i + 1;
      counter++;

      (n_dof_tot[HECMW_MESH_DOF_THREE])++;
      MASK_BIT(*dof_flag, BIT_DOF_THREE);
      CLEAR_BIT(node_flag[i], BIT_DOF_ALL);
    }
  }

  /* two DOF */
  for (i = 0; i < local_mesh->n_node; i++) {
    if (EVAL_BIT(node_flag[i], BIT_DOF_TWO)) {
      node_old2new[i]       = counter + 1;
      node_new2old[counter] = i + 1;
      counter++;

      (n_dof_tot[HECMW_MESH_DOF_TWO])++;
      MASK_BIT(*dof_flag, BIT_DOF_TWO);
      CLEAR_BIT(node_flag[i], BIT_DOF_ALL);
    }
  }

  HECMW_assert(counter == local_mesh->n_node);

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  reorder node according to DOF with EQUATION_BLOCK                         */
/*----------------------------------------------------------------------------*/
#if 0  /* commented out by K.Goto; begin */
static int
mask_eqn_block( struct hecmwST_local_mesh *local_mesh,
                char *node_flag, char *block_flag, const int eqn_block_idx )
{
  struct hecmwST_node_grp *grp=local_mesh->node_group;
  int dof_max;
  int i, j, js;

  for( js=0, i=grp->grp_index[eqn_block_idx]; i<grp->grp_index[eqn_block_idx+1]; i++ ) {
    for( dof_max=0, j=js; j<grp->grp_item[i]; j++ ) {
      if( EVAL_BIT( node_flag[j], BIT_DOF_TWO ) )   MASK_BIT( block_flag[i-grp->grp_index[eqn_block_idx]], BIT_DOF_TWO );
      if( EVAL_BIT( node_flag[j], BIT_DOF_THREE ) ) MASK_BIT( block_flag[i-grp->grp_index[eqn_block_idx]], BIT_DOF_THREE );
      if( EVAL_BIT( node_flag[j], BIT_DOF_SIX ) )   MASK_BIT( block_flag[i-grp->grp_index[eqn_block_idx]], BIT_DOF_SIX );
      if( EVAL_BIT( node_flag[j], BIT_DOF_FOUR) )   MASK_BIT( block_flag[i-grp->grp_index[eqn_block_idx]], BIT_DOF_FOUR );
    }
    js = grp->grp_item[i];
  }

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int
reorder_node_dof_4mpc_inner( struct hecmwST_local_mesh *local_mesh,
                             char *node_flag, char *block_flag,
                             int *node_new2old, int *node_old2new, int *block_old2new,
                             const int n_eqn_block, const int eqn_block_idx,
                             char *dof_flag, int *n_dof_tot )
{
  struct hecmwST_node_grp *grp=local_mesh->node_group;
  int idx=grp->grp_index[eqn_block_idx];
  int counter=0, blocks=0;
  int i, j, js;

  /* six DOF */
  for( js=0, i=0; i<n_eqn_block; i++ ) {
    if( EVAL_BIT( block_flag[i], BIT_DOF_SIX ) ) {
      block_old2new[i] = blocks++;

      for( j=js; j<grp->grp_item[idx+i]; j++ ) {
        node_old2new[j]       = counter+1;
        node_new2old[counter] = j+1;
        counter++;

        (n_dof_tot[HECMW_MESH_DOF_SIX])++;
      }
      MASK_BIT( *dof_flag, BIT_DOF_SIX );
      CLEAR_BIT( block_flag[i], BIT_DOF_ALL );
    }
    js = grp->grp_item[idx+i];
  }

  /* four DOF */
  for( js=0, i=0; i<n_eqn_block; i++ ) {
    if( EVAL_BIT( block_flag[i], BIT_DOF_FOUR ) ) {
      block_old2new[i] = blocks++;

      for( j=js; j<grp->grp_item[idx+i]; j++ ) {
        node_old2new[j]       = counter+1;
        node_new2old[counter] = j+1;
        counter++;

        (n_dof_tot[HECMW_MESH_DOF_FOUR])++;
      }
      MASK_BIT( *dof_flag, BIT_DOF_FOUR );
      CLEAR_BIT( block_flag[i], BIT_DOF_ALL );
    }
    js = grp->grp_item[idx+i];
  }

  /* three DOF */
  for( js=0, i=0; i<n_eqn_block; i++ ) {
    if( EVAL_BIT( block_flag[i], BIT_DOF_THREE ) ) {
      block_old2new[i] = blocks++;

      for( j=js; j<grp->grp_item[idx+i]; j++ ) {
        node_old2new[j]       = counter+1;
        node_new2old[counter] = j+1;
        counter++;

        (n_dof_tot[HECMW_MESH_DOF_THREE])++;
      }
      MASK_BIT( *dof_flag, BIT_DOF_THREE );
      CLEAR_BIT( block_flag[i], BIT_DOF_ALL );
    }
    js = grp->grp_item[idx+i];
  }

  /* two DOF */
  for( js=0, i=0; i<n_eqn_block; i++ ) {
    if( EVAL_BIT( block_flag[i], BIT_DOF_TWO ) ) {
      block_old2new[i] = blocks++;

      for( j=js; j<grp->grp_item[idx+i]; j++ ) {
        node_old2new[j]       = counter+1;
        node_new2old[counter] = j+1;
        counter++;

        (n_dof_tot[HECMW_MESH_DOF_TWO])++;
      }
      MASK_BIT( *dof_flag, BIT_DOF_TWO );
      CLEAR_BIT( block_flag[i], BIT_DOF_ALL );
    }
    js = grp->grp_item[idx+i];
  }

  HECMW_assert( counter == local_mesh->n_node );

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int
reorder_node_dof_4mpc( struct hecmwST_local_mesh *local_mesh,
                       char *node_flag, int *node_new2old, int *node_old2new,
                       char *dof_flag, int *n_dof_tot )
{
  struct hecmwST_node_grp *grp=local_mesh->node_group;
  char *block_flag;
  int *block_old2new;
  int n_eqn_block, eqn_block_idx;

  /* group id of EQUATION_BLOCK */
  eqn_block_idx = get_eqn_block_idx( local_mesh );
  if( eqn_block_idx < 0 ) {
    HECMW_print_msg( HECMW_LOG_WARN, HECMW_COMMON_W_NO_EQN_BLOCK, "");
    return 1;
  }

  /* number of EQUATION_BLOCKs */
  n_eqn_block = grp->grp_index[eqn_block_idx+1] - grp->grp_index[eqn_block_idx];

  /* allocation */
  block_flag = (char *)HECMW_calloc( n_eqn_block, sizeof(char) );
  if( block_flag == NULL ) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  block_old2new = (int *)HECMW_malloc( sizeof(int)*n_eqn_block );
  if( block_old2new == NULL ) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  /* mask EQUATION_BLOCK according to DOF */
  if(mask_eqn_block( local_mesh, node_flag, block_flag, eqn_block_idx )) return -1;

  /* reorder nodes according to DOF with EQUATION_BLOCK */
  if(reorder_node_dof_4mpc_inner( local_mesh, node_flag, block_flag, node_new2old, node_old2new, block_old2new, n_eqn_block, eqn_block_idx, dof_flag, n_dof_tot )) return -1;

  /* reorder relevant arrays */
  if(old2new_eqn_block_info( local_mesh, block_old2new )) return -1;

  /* free */
  HECMW_free( block_flag );
  HECMW_free( block_old2new );

  return 0;
}
#endif /* commented out by K.Goto; end */

/*----------------------------------------------------------------------------*/
/*  count number of DOF groups & set maximal number of DOF                    */
/*----------------------------------------------------------------------------*/
static int count_n_dof_grp(struct hecmwST_local_mesh *local_mesh,
                           char *dof_flag) {
  local_mesh->n_dof_grp = 0;

  /* two DOF */
  if (EVAL_BIT(*dof_flag, BIT_DOF_TWO)) {
    (local_mesh->n_dof_grp)++;
    local_mesh->n_dof = HECMW_MESH_DOF_TWO;
  }

  /* three DOF */
  if (EVAL_BIT(*dof_flag, BIT_DOF_THREE)) {
    (local_mesh->n_dof_grp)++;
    local_mesh->n_dof = HECMW_MESH_DOF_THREE;
  }

  /* four DOF */
  if (EVAL_BIT(*dof_flag, BIT_DOF_FOUR)) {
    (local_mesh->n_dof_grp)++;
    local_mesh->n_dof = HECMW_MESH_DOF_FOUR;
  }

  /* six DOF */
  if (EVAL_BIT(*dof_flag, BIT_DOF_SIX)) {
    (local_mesh->n_dof_grp)++;
    local_mesh->n_dof = HECMW_MESH_DOF_SIX;
  }

  HECMW_assert(local_mesh->n_dof_grp > 0 &&
               local_mesh->n_dof_grp <= HECMW_MESH_DOF_TOT);

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  create index for DOF group & its value                                    */
/*----------------------------------------------------------------------------*/
static int create_node_dof_item(struct hecmwST_local_mesh *local_mesh,
                                char *dof_flag, int *n_dof_tot) {
  int counter = 0;

  /* allocation */
  if (local_mesh->node_dof_index) {
    HECMW_free(local_mesh->node_dof_index);
  }
  local_mesh->node_dof_index =
      (int *)HECMW_calloc(local_mesh->n_dof_grp + 1, sizeof(int));
  if (local_mesh->node_dof_index == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  if (local_mesh->node_dof_item) {
    HECMW_free(local_mesh->node_dof_item);
  }
  local_mesh->node_dof_item =
      (int *)HECMW_malloc(sizeof(int) * local_mesh->n_dof_grp);
  if (local_mesh->node_dof_item == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  /* six DOF */
  if (EVAL_BIT(*dof_flag, BIT_DOF_SIX)) {
    local_mesh->node_dof_index[counter + 1] =
        local_mesh->node_dof_index[counter] + n_dof_tot[HECMW_MESH_DOF_SIX];
    local_mesh->node_dof_item[counter] = HECMW_MESH_DOF_SIX;
    counter++;
  }

  /* four DOF */
  if (EVAL_BIT(*dof_flag, BIT_DOF_FOUR)) {
    local_mesh->node_dof_index[counter + 1] =
        local_mesh->node_dof_index[counter] + n_dof_tot[HECMW_MESH_DOF_FOUR];
    local_mesh->node_dof_item[counter] = HECMW_MESH_DOF_FOUR;
    counter++;
  }

  /* three DOF */
  if (EVAL_BIT(*dof_flag, BIT_DOF_THREE)) {
    local_mesh->node_dof_index[counter + 1] =
        local_mesh->node_dof_index[counter] + n_dof_tot[HECMW_MESH_DOF_THREE];
    local_mesh->node_dof_item[counter] = HECMW_MESH_DOF_THREE;
    counter++;
  }

  /* two DOF */
  if (EVAL_BIT(*dof_flag, BIT_DOF_TWO)) {
    local_mesh->node_dof_index[counter + 1] =
        local_mesh->node_dof_index[counter] + n_dof_tot[HECMW_MESH_DOF_TWO];
    local_mesh->node_dof_item[counter] = HECMW_MESH_DOF_TWO;
    counter++;
  }

  HECMW_assert(counter == local_mesh->n_dof_grp);
  HECMW_assert(local_mesh->node_dof_index[local_mesh->n_dof_grp] ==
               local_mesh->n_node);

  return 0;
}

/*============================================================================*/
/*  reorder nodes according to DOF                                            */
/*============================================================================*/
extern int HECMW_reorder_node_dof(struct hecmwST_local_mesh *local_mesh) {
  int *node_new2old,
      *node_old2new; /* conversion table between old and new node id */
  char *node_flag;   /* for masking */
  char dof_flag = '\0';
  int n_dof_tot[HECMW_MESH_DOF_MAX + 1];
  int i;

  HECMW_assert(local_mesh);

  /* allocation */
  node_flag = (char *)HECMW_calloc(local_mesh->n_node, sizeof(char));
  if (node_flag == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  node_new2old = (int *)HECMW_malloc(sizeof(int) * local_mesh->n_node);
  if (node_new2old == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  node_old2new = (int *)HECMW_malloc(sizeof(int) * local_mesh->n_node);
  if (node_old2new == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  for (i = 0; i < HECMW_MESH_DOF_MAX + 1; i++) {
    n_dof_tot[i] = 0;
  }

  /* mask node */
  if (mask_node_dof(local_mesh, node_flag)) {
    return -1;
  }

  /* reorder node */
  /* commented out by K.Goto; begin */
  /*
  if( local_mesh->mpc->n_mpc ) {
    int rtc;
    if((rtc = reorder_node_dof_4mpc( local_mesh, node_flag, node_new2old,
  node_old2new, &dof_flag, n_dof_tot )) < 0) {
      return -1;
        }
    if( rtc ) {
      if(reorder_node_dof( local_mesh, node_flag, node_new2old, node_old2new,
  &dof_flag, n_dof_tot )) {
        return -1;
          }
    }
  } else {
  */
  /* commented out by K.Goto; end */
  if (reorder_node_dof(local_mesh, node_flag, node_new2old, node_old2new,
                       &dof_flag, n_dof_tot)) {
    return -1;
  }
  /* commented out by K.Goto; begin */
  /*
  }
  */
  /* commented out by K.Goto; end */

  HECMW_free(node_flag);

  /* create DOF information */
  if (count_n_dof_grp(local_mesh, &dof_flag)) {
    return -1;
  }

  if (create_node_dof_item(local_mesh, &dof_flag, n_dof_tot)) {
    return -1;
  }

  /* reorder relevant arrays */
  if (old2new_node_info(local_mesh, node_new2old, node_old2new)) {
    return -1;
  }

  /* free */
  HECMW_free(node_new2old);
  HECMW_free(node_old2new);

  return 0;
}

/*                                                                            */
/*  reorder node according to MPC group                                       */
/*                                                                            */

/*----------------------------------------------------------------------------*/
/*  mask MPC group                                                            */
/*----------------------------------------------------------------------------*/
static int set_mpc_block(struct hecmwST_local_mesh *local_mesh,
                         int *mpc_node_flag, int *mpc_index_flag) {
  int node, min_group, group;
  int i, j, js, je;

  for (i = 0; i < local_mesh->mpc->n_mpc; i++) {
    js = local_mesh->mpc->mpc_index[i];
    je = local_mesh->mpc->mpc_index[i + 1];

    /* MPC???롼?פγ???????????륰?롼?פ?mpc_node_flag ?˳??   */
    /* ??PC???롼?פν????????å???mpc_index_flag ?˳?? (??? */
    min_group = i;
    for (j = js; j < je; j++) {
      node = local_mesh->mpc->mpc_item[j];

      /* ?????Ƥ?????????????Υ??롼?פ˽?????Ƥ????? ??group >= 0 */
      /*                       ¾?Υ??롼?פ˽?????Ƥ??ʤ?????group <  0 */
      group                   = mpc_node_flag[node - 1];
      mpc_node_flag[node - 1] = (group < 0) ? i : group;

      /* i????PC???롼?פγ????ν?????륰?롼?פκǾ????롼?פ?? */
      min_group = (mpc_index_flag[mpc_node_flag[node - 1]] < min_group)
                      ? mpc_index_flag[mpc_node_flag[node - 1]]
                      : min_group;
    }

    /* i????PC???롼?פγ???????????Ƥ??륰?롼?פϺǾ????롼?פ????*/
    for (j = js; j < je; j++) {
      node  = local_mesh->mpc->mpc_item[j];
      group = mpc_node_flag[node - 1];

      mpc_index_flag[group] = min_group;
    }
    mpc_index_flag[i] = min_group;
  }

  /* ???롼?פ?????????*/
  for (i = 0; i < local_mesh->mpc->n_mpc; i++) {
    group             = mpc_index_flag[i];
    mpc_index_flag[i] = mpc_index_flag[group];
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  count MPC blocks                                                          */
/*----------------------------------------------------------------------------*/
static int count_mpc_block(struct hecmwST_local_mesh *local_mesh,
                           struct equation_block *eqn_block,
                           int *mpc_index_flag, int *mpc_group2block) {
  int *n_block;
  int block, counter;
  int i;

  /* allocation */
  n_block = (int *)HECMW_calloc(local_mesh->mpc->n_mpc, sizeof(int));
  if (n_block == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  /* count MPC groups in each MPC block */
  for (i = 0; i < local_mesh->mpc->n_mpc; i++) {
    block = mpc_index_flag[i];
    n_block[block]++;
  }

  /* count MPC blocks                                              */
  /* create conversion table from "MPC group ID" to "MPC block ID" */
  for (counter = 0, i = 0; i < local_mesh->mpc->n_mpc; i++) {
    if (n_block[i]) mpc_group2block[i] = counter++;
  }

  /* number of MPC blocks */
  eqn_block->n_mpc_block = counter;
  eqn_block->n_eqn_block = counter;

  /* free */
  HECMW_free(n_block);

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  count EQUATION_BLOCKs                                                     */
/*----------------------------------------------------------------------------*/
static int count_eqn_block(struct hecmwST_local_mesh *local_mesh,
                           struct equation_block *eqn_block,
                           int *mpc_node_flag) {
  int i;

  for (i = 0; i < local_mesh->n_node; i++) {
    if (mpc_node_flag[i] < 0) (eqn_block->n_eqn_block)++;
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  set EQUATION_BLOCK id to node                                             */
/*----------------------------------------------------------------------------*/
static int set_eqn_block_of_node(struct hecmwST_local_mesh *local_mesh,
                                 struct equation_block *eqn_block,
                                 int *mpc_node_flag, int *mpc_index_flag,
                                 int *mpc_group2block) {
  int counter;
  int i;

  for (counter = eqn_block->n_mpc_block, i = 0; i < local_mesh->n_node; i++) {
    if (mpc_node_flag[i] >= 0) {
      mpc_node_flag[i] = mpc_group2block[mpc_index_flag[mpc_node_flag[i]]];
    } else {
      mpc_node_flag[i] = counter++;
    }
  }

  HECMW_assert(counter == eqn_block->n_eqn_block);

  return 0;
}

/*----------------------------------------------------------------------------*/
/*   create EQUATION_BLOCK                                                    */
/*----------------------------------------------------------------------------*/
static int create_eqn_block_index(struct hecmwST_local_mesh *local_mesh,
                                  struct equation_block *eqn_block,
                                  int *mpc_node_flag, int *mpc_group2block) {
  int *n_block;
  int i;

  /* allocation */
  n_block = (int *)HECMW_calloc(eqn_block->n_eqn_block, sizeof(int));
  if (n_block == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  eqn_block->eqn_block_index =
      (int *)HECMW_calloc(eqn_block->n_eqn_block + 1, sizeof(int));
  if (eqn_block->eqn_block_index == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  /* count nodes in each EQUATION_BLOCK */
  for (i = 0; i < local_mesh->n_node; i++) {
    (n_block[mpc_node_flag[i]])++;
  }

  /* create index for EQUATION_BLOCK */
  for (i = 0; i < eqn_block->n_eqn_block; i++) {
    eqn_block->eqn_block_index[i + 1] =
        eqn_block->eqn_block_index[i] + n_block[i];
  }

  /* free */
  HECMW_free(n_block);

  HECMW_assert(eqn_block->eqn_block_index[eqn_block->n_eqn_block] ==
               local_mesh->n_node);

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  create conversion table of node                                           */
/*----------------------------------------------------------------------------*/
static int create_old2new_node(struct hecmwST_local_mesh *local_mesh,
                               struct equation_block *eqn_block,
                               int *mpc_node_flag, int *node_old2new,
                               int *node_new2old) {
  int *n_block;
  int new_id;
  int i;

  n_block = (int *)HECMW_calloc(eqn_block->n_eqn_block, sizeof(int));
  if (n_block == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  for (i = 0; i < local_mesh->n_node; i++) {
    new_id = eqn_block->eqn_block_index[mpc_node_flag[i]] +
             n_block[mpc_node_flag[i]];
    node_old2new[i]      = new_id + 1;
    node_new2old[new_id] = i + 1;
    (n_block[mpc_node_flag[i]])++;
  }

  HECMW_free(n_block);

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  reconstruct node group information                                        */
/*----------------------------------------------------------------------------*/
static int reconstruct_node_grp(struct hecmwST_local_mesh *local_mesh,
                                struct equation_block *eqn_block) {
  struct hecmwST_node_grp *grp = local_mesh->node_group;
  int i;

  /* number of node groups */
  (grp->n_grp)++;

  /* index for node group */
  grp->grp_index =
      (int *)HECMW_realloc(grp->grp_index, sizeof(int) * (grp->n_grp + 1));
  if (grp->grp_index == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  grp->grp_index[grp->n_grp] =
      grp->grp_index[grp->n_grp - 1] + eqn_block->n_eqn_block;

  /* name of node group */
  grp->grp_name =
      (char **)HECMW_realloc(grp->grp_name, sizeof(char *) * grp->n_grp);
  if (grp->grp_name == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }
  grp->grp_name[grp->n_grp - 1] =
      (char *)HECMW_malloc(sizeof(char *) * (HECMW_NAME_LEN + 1));
  if (grp->grp_name[grp->n_grp - 1] == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  strcpy(grp->grp_name[grp->n_grp - 1], HECMW_COMMON_EQUATION_BLOCK_NAME);

  /* member of node group */
  grp->grp_item = (int *)HECMW_realloc(
      grp->grp_item, sizeof(int) * grp->grp_index[grp->n_grp]);
  if (grp->grp_item == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  for (i = 0; i < eqn_block->n_eqn_block; i++) {
    grp->grp_item[grp->grp_index[grp->n_grp - 1] + i] =
        eqn_block->eqn_block_index[i + 1];
  }

  return 0;
}

/*============================================================================*/
/*  reorder node according to MPC group                                       */
/*============================================================================*/
extern int HECMW_reorder_node_mpc(struct hecmwST_local_mesh *local_mesh) {
  int *mpc_node_flag, *mpc_index_flag, *mpc_group2block;
  int *node_old2new, *node_new2old;
  int i;
  struct equation_block *eqn_block;

  if (local_mesh->mpc->n_mpc == 0) return 0;

  /* allocation */
  mpc_node_flag = (int *)HECMW_calloc(local_mesh->n_node, sizeof(int));
  if (mpc_node_flag == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }
  for (i = 0; i < local_mesh->n_node; i++) {
    mpc_node_flag[i] = -1;
  }

  mpc_index_flag = (int *)HECMW_malloc(sizeof(int) * local_mesh->mpc->n_mpc);
  if (mpc_index_flag == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }
  for (i = 0; i < local_mesh->mpc->n_mpc; i++) {
    mpc_index_flag[i] = i;
  }

  eqn_block =
      (struct equation_block *)HECMW_malloc(sizeof(struct equation_block));
  if (eqn_block == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  /* group together */
  if (set_mpc_block(local_mesh, mpc_node_flag, mpc_index_flag)) {
    return -1;
  }

  /* count MPC blocks */
  mpc_group2block = (int *)HECMW_malloc(sizeof(int) * local_mesh->mpc->n_mpc);
  if (mpc_group2block == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }
  for (i = 0; i < local_mesh->mpc->n_mpc; i++) mpc_group2block[i] = -1;

  if (count_mpc_block(local_mesh, eqn_block, mpc_index_flag, mpc_group2block)) {
    return -1;
  }

  /* count EQUATION_BLOCKs */
  if (count_eqn_block(local_mesh, eqn_block, mpc_node_flag)) {
    return -1;
  }

  /* set EQUATION_BLOCK to node */
  if (set_eqn_block_of_node(local_mesh, eqn_block, mpc_node_flag,
                            mpc_index_flag, mpc_group2block)) {
    return -1;
  }

  /* create EQUATION_BLOCK */
  if (create_eqn_block_index(local_mesh, eqn_block, mpc_node_flag,
                             mpc_group2block)) {
    return -1;
  }

  HECMW_free(mpc_index_flag);
  HECMW_free(mpc_group2block);

  /* create conversion table between old and new node id */
  node_old2new = (int *)HECMW_malloc(sizeof(int) * local_mesh->n_node);
  if (node_old2new == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  node_new2old = (int *)HECMW_malloc(sizeof(int) * local_mesh->n_node);
  if (node_new2old == NULL) {
    HECMW_set_error(HECMW_COMMON_E_ALLOCATION, "");
    return -1;
  }

  if (create_old2new_node(local_mesh, eqn_block, mpc_node_flag, node_old2new,
                          node_new2old)) {
    return -1;
  }

  if (old2new_node_info(local_mesh, node_new2old, node_old2new)) {
    return -1;
  }

  HECMW_free(mpc_node_flag);
  HECMW_free(node_old2new);
  HECMW_free(node_new2old);

  /* reconstruct node group */
  if (reconstruct_node_grp(local_mesh, eqn_block)) {
    return -1;
  }

  HECMW_free(eqn_block->eqn_block_index);
  HECMW_free(eqn_block);

  return 0;
}

/*                                                                            */
/*  reorder node & element                                                    */
/*                                                                            */

extern int HECMW_reorder(struct hecmwST_local_mesh *local_mesh) {
  HECMW_assert(local_mesh);

  /* reorder element according to finite element type */
  if (HECMW_reorder_elem_type(local_mesh)) return -1;

  /* reorder node according to MPC group */
  /* commented out by K.Goto; begin */
  /* if(HECMW_reorder_node_mpc( local_mesh )) return -1; */
  /* commented out by K.Goto; end */

  /* reorder node according to node dof */
  if (HECMW_reorder_node_dof(local_mesh)) return -1;

  return 0;
}
