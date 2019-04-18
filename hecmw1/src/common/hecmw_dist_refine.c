/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "hecmw_struct.h"
#include "hecmw_util.h"
#include "hecmw_common_define.h"
#include "hecmw_etype.h"
#include "hecmw_bit_array.h"
#include "hecmw_set_int.h"
#include "hecmw_varray_int.h"
#include "hecmw_dist_alloc.h"
#include "hecmw_dist_free.h"
#include "hecmw_dist_refine.h"

#ifdef HECMW_WITH_REFINER

#include "rcapRefiner.h"

/*============================================================================*/
/*                                                                            */
/*  convert elem_node_item                                                    */
/*                                                                            */
/*============================================================================*/
static const int *get_enode_h2r(int etype, int *ierror) {
  /* tetra, prism, hexa: no need to convert */
  static const int enode_pyr_h2r[13] = {4,  0,  1, 2, 3, 9, 10,
                                        11, 12, 5, 6, 7, 8};
  static const int enode_rod2_h2r[3] = {0, 2, 1};
  const int *enode_h2r;

  *ierror = 0;

  if (HECMW_is_etype_link(etype)) return NULL;
  if (HECMW_is_etype_truss(etype)) return NULL;

  switch (etype) {
    case HECMW_ETYPE_ROD1:
    /* case HECMW_ETYPE_ROD31: */
    case HECMW_ETYPE_TRI1:
    case HECMW_ETYPE_TRI2:
    case HECMW_ETYPE_QUA1:
    case HECMW_ETYPE_QUA2:
    case HECMW_ETYPE_TET1:
    case HECMW_ETYPE_TET1_4:
    case HECMW_ETYPE_TET2:
    case HECMW_ETYPE_PRI1:
    case HECMW_ETYPE_PRI2:
    case HECMW_ETYPE_HEX1:
    case HECMW_ETYPE_HEX1_4:
    case HECMW_ETYPE_HEX2:
    case HECMW_ETYPE_BEM1:
    case HECMW_ETYPE_BEM3:
    case HECMW_ETYPE_SHT1:
    case HECMW_ETYPE_SHT2:
    case HECMW_ETYPE_SHT6:
    case HECMW_ETYPE_SHQ1:
    case HECMW_ETYPE_SHQ2:
    case HECMW_ETYPE_SHQ8:
    case HECMW_ETYPE_PTT1:
    case HECMW_ETYPE_PTT2:
    case HECMW_ETYPE_PTQ1:
    case HECMW_ETYPE_PTQ2:
      enode_h2r = NULL;
      break;
    case HECMW_ETYPE_PYR1:
    case HECMW_ETYPE_PYR2:
      enode_h2r = enode_pyr_h2r;
      break;
    case HECMW_ETYPE_ROD2:
    case HECMW_ETYPE_BEM2:
      enode_h2r = enode_rod2_h2r;
      break;
    default:
      HECMW_log(HECMW_LOG_ERROR,
                "Element type %d not supported for rerinement.\n", etype);
      *ierror   = 1;
      enode_h2r = NULL;
      break;
  }
  return enode_h2r;
}

static const int *get_enode_r2h(int etype, int *ierror) {
  /* tetra, prism, hexa: no need to convert */
  static const int enode_pyr_r2h[13] = {1,  2,  3, 4, 0, 9, 10,
                                        11, 12, 5, 6, 7, 8};
  static const int enode_rod2_r2h[3] = {0, 2, 1};
  const int *enode_r2h;

  *ierror = 0;

  if (HECMW_is_etype_link(etype)) return NULL;
  if (HECMW_is_etype_truss(etype)) return NULL;

  switch (etype) {
    case HECMW_ETYPE_ROD1:
    case HECMW_ETYPE_ROD31:
    case HECMW_ETYPE_TRI1:
    case HECMW_ETYPE_TRI2:
    case HECMW_ETYPE_QUA1:
    case HECMW_ETYPE_QUA2:
    case HECMW_ETYPE_TET1:
    case HECMW_ETYPE_TET1_4:
    case HECMW_ETYPE_TET2:
    case HECMW_ETYPE_PRI1:
    case HECMW_ETYPE_PRI2:
    case HECMW_ETYPE_HEX1:
    case HECMW_ETYPE_HEX1_4:
    case HECMW_ETYPE_HEX2:
    case HECMW_ETYPE_BEM1:
    case HECMW_ETYPE_BEM3:
    case HECMW_ETYPE_SHT1:
    case HECMW_ETYPE_SHT2:
    case HECMW_ETYPE_SHT6:
    case HECMW_ETYPE_SHQ1:
    case HECMW_ETYPE_SHQ2:
    case HECMW_ETYPE_SHQ8:
    case HECMW_ETYPE_PTT1:
    case HECMW_ETYPE_PTT2:
    case HECMW_ETYPE_PTQ1:
    case HECMW_ETYPE_PTQ2:
      enode_r2h = NULL;
      break;
    case HECMW_ETYPE_PYR1:
    case HECMW_ETYPE_PYR2:
      enode_r2h = enode_pyr_r2h;
      break;
    case HECMW_ETYPE_ROD2:
    case HECMW_ETYPE_BEM2:
      enode_r2h = enode_rod2_r2h;
      break;
    default:
      HECMW_log(HECMW_LOG_ERROR,
                "Element type %d not supported for refinement.\n", etype);
      *ierror   = 1;
      enode_r2h = NULL;
  }
  return enode_r2h;
}

static int elem_node_item_hecmw2rcap(struct hecmwST_local_mesh *mesh) {
  int i, j;
  int etype, ierror;
  const int *enode_h2r;
  int *enode;
  int nn, istart, iend;
  int temp_enode[HECMW_MAX_NODE_MAX];

  for (i = 0; i < mesh->n_elem_type; i++) {
    etype     = mesh->elem_type_item[i];
    enode_h2r = get_enode_h2r(etype, &ierror);
    if (ierror) return HECMW_ERROR;
    if (enode_h2r == NULL) continue;
    nn     = HECMW_get_max_node(etype);
    istart = mesh->elem_type_index[i];
    iend   = mesh->elem_type_index[i + 1];
    for (j = istart; j < iend; j++) {
      enode = mesh->elem_node_item + mesh->elem_node_index[i];
      for (j = 0; j < nn; j++) temp_enode[j] = enode[enode_h2r[j]];
      for (j = 0; j < nn; j++) enode[j] = temp_enode[j];
    }
  }
  return HECMW_SUCCESS;
}

static int elem_node_item_rcap2hecmw(struct hecmwST_local_mesh *mesh) {
  int i, j;
  int etype, ierror;
  const int *enode_r2h;
  int *enode;
  int nn, istart, iend;
  int temp_enode[HECMW_MAX_NODE_MAX];

  for (i = 0; i < mesh->n_elem_type; i++) {
    etype     = mesh->elem_type_item[i];
    enode_r2h = get_enode_r2h(etype, &ierror);
    if (ierror) return HECMW_ERROR;
    if (enode_r2h == NULL) continue;
    nn     = HECMW_get_max_node(etype);
    istart = mesh->elem_type_index[i];
    iend   = mesh->elem_type_index[i + 1];
    for (j = istart; j < iend; j++) {
      enode = mesh->elem_node_item + mesh->elem_node_index[i];
      for (j = 0; j < nn; j++) temp_enode[j] = enode[enode_r2h[j]];
      for (j = 0; j < nn; j++) enode[j] = temp_enode[j];
    }
  }
  return HECMW_SUCCESS;
}

/*============================================================================*/
/*                                                                            */
/*  convert surface IDs                                                       */
/*                                                                            */
/*============================================================================*/
static const int *get_sid_h2r(int etype, int *ierror) {
  /* HECMW to RCAP */
  /*
  static const int sid_tri_h2r[4] = {-1, 0, 1, 2};
  static const int sid_qua_h2r[5] = {-1, 3, 1, 0, 2};
  static const int sid_tet_h2r[5] = {-1, 0, 1, 2, 3};
  static const int sid_pri_h2r[6] = {-1, 3, 4, 2, 0, 1};
  static const int sid_hex_h2r[7] = {-1, 5, 3, 2, 4, 0, 1};
  static const int sid_pyr_h2r[6] = {-1, 3, 1, 0, 2, 4};
  */
  /* FSTR to RCAP*/
  static const int sid_tri_h2r[4] = {-1, 2, 0, 1};
  static const int sid_qua_h2r[5] = {-1, 0, 1, 2, 3};
  static const int sid_tet_h2r[5] = {-1, 3, 2, 0, 1};
  static const int sid_pri_h2r[6] = {-1, 0, 1, 2, 3, 4};
  static const int sid_hex_h2r[7] = {-1, 0, 1, 2, 3, 4, 5};
  static const int sid_pyr_h2r[6] = {-1, 3, 1,
                                     0,  2, 4}; /* TEMPORARY: same as HECMW */

  const int *sid_h2r;
  *ierror = 0;
  switch (etype) {
    case HECMW_ETYPE_TRI1:
    case HECMW_ETYPE_TRI2:
      sid_h2r = sid_tri_h2r;
      break;
    case HECMW_ETYPE_QUA1:
    case HECMW_ETYPE_QUA2:
      sid_h2r = sid_qua_h2r;
      break;
    case HECMW_ETYPE_TET1:
    case HECMW_ETYPE_TET1_4:
    case HECMW_ETYPE_TET2:
      sid_h2r = sid_tet_h2r;
      break;
    case HECMW_ETYPE_PRI1:
    case HECMW_ETYPE_PRI2:
      sid_h2r = sid_pri_h2r;
      break;
    case HECMW_ETYPE_HEX1:
    case HECMW_ETYPE_HEX1_4:
    case HECMW_ETYPE_HEX2:
      sid_h2r = sid_hex_h2r;
      break;
    case HECMW_ETYPE_PYR1:
    case HECMW_ETYPE_PYR2:
      sid_h2r = sid_pyr_h2r;
      break;
    case HECMW_ETYPE_SHT1:
    case HECMW_ETYPE_SHT2:
    case HECMW_ETYPE_SHT6:
    case HECMW_ETYPE_SHQ1:
    case HECMW_ETYPE_SHQ2:
    case HECMW_ETYPE_SHQ8:
    case HECMW_ETYPE_PTT1:
    case HECMW_ETYPE_PTT2:
    case HECMW_ETYPE_PTQ1:
    case HECMW_ETYPE_PTQ2:
      sid_h2r = NULL;
      break;
    default:
      HECMW_log(HECMW_LOG_ERROR,
                "Element type %d not supported for refinement.\n", etype);
      *ierror = 1;
      sid_h2r = NULL;
  }
  return sid_h2r;
}

static const int *get_sid_r2h(int etype, int *ierror) {
  /* RCAP to HECMW */
  /*
  static const int sid_tri_r2h[3] = {1, 2, 3};
  static const int sid_qua_r2h[4] = {3, 2, 4, 1};
  static const int sid_tet_r2h[4] = {1, 2, 3, 4};
  static const int sid_pri_r2h[5] = {4, 5, 3, 1, 2};
  static const int sid_hex_r2h[6] = {1, 2, 3, 4, 5, 6};
  static const int sid_pyr_r2h[5] = {3, 2, 4, 1, 5};
  */
  /* RCAP to FSTR */
  static const int sid_tri_r2h[3] = {2, 3, 1};
  static const int sid_qua_r2h[4] = {1, 2, 3, 4};
  static const int sid_tet_r2h[4] = {3, 4, 2, 1};
  static const int sid_pri_r2h[5] = {1, 2, 3, 4, 5};
  static const int sid_hex_r2h[6] = {1, 2, 3, 4, 5, 6};
  static const int sid_pyr_r2h[5] = {3, 2, 4, 1,
                                     5}; /* TEMPORARY: same as HECMW */

  const int *sid_r2h;
  *ierror = 0;
  switch (etype) {
    case HECMW_ETYPE_TRI1:
    case HECMW_ETYPE_TRI2:
      sid_r2h = sid_tri_r2h;
      break;
    case HECMW_ETYPE_QUA1:
    case HECMW_ETYPE_QUA2:
      sid_r2h = sid_qua_r2h;
      break;
    case HECMW_ETYPE_TET1:
    case HECMW_ETYPE_TET1_4:
    case HECMW_ETYPE_TET2:
      sid_r2h = sid_tet_r2h;
      break;
    case HECMW_ETYPE_PRI1:
    case HECMW_ETYPE_PRI2:
      sid_r2h = sid_pri_r2h;
      break;
    case HECMW_ETYPE_HEX1:
    case HECMW_ETYPE_HEX1_4:
    case HECMW_ETYPE_HEX2:
      sid_r2h = sid_hex_r2h;
      break;
    case HECMW_ETYPE_PYR1:
    case HECMW_ETYPE_PYR2:
      sid_r2h = sid_pyr_r2h;
      break;
    case HECMW_ETYPE_SHT1:
    case HECMW_ETYPE_SHT2:
    case HECMW_ETYPE_SHT6:
    case HECMW_ETYPE_SHQ1:
    case HECMW_ETYPE_SHQ2:
    case HECMW_ETYPE_SHQ8:
    case HECMW_ETYPE_PTT1:
    case HECMW_ETYPE_PTT2:
    case HECMW_ETYPE_PTQ1:
    case HECMW_ETYPE_PTQ2:
      sid_r2h = NULL;
      break;
    default:
      HECMW_log(HECMW_LOG_ERROR,
                "Element type %d not supported for refinement.\n", etype);
      *ierror = 1;
      sid_r2h = NULL;
  }
  return sid_r2h;
}

static int surf_ID_hecmw2rcap(struct hecmwST_local_mesh *mesh) {
  struct hecmwST_surf_grp *grp = mesh->surf_group;
  int i;

  if (grp->n_grp == 0) return HECMW_SUCCESS;

  for (i = 0; i < grp->grp_index[grp->n_grp]; i++) {
    int eid    = grp->grp_item[2 * i];
    int *sid_p = &(grp->grp_item[2 * i + 1]);
    const int *sid_h2r;
    int ierror;

    HECMW_assert(0 < eid && eid <= mesh->n_elem_gross);

    sid_h2r = get_sid_h2r(mesh->elem_type[eid - 1], &ierror);
    if (ierror) return HECMW_ERROR;
    if (sid_h2r) *sid_p = sid_h2r[*sid_p];
  }
  return HECMW_SUCCESS;
}

static int surf_ID_rcap2hecmw(struct hecmwST_local_mesh *mesh) {
  struct hecmwST_surf_grp *grp = mesh->surf_group;
  int i;

  if (grp->n_grp == 0) return HECMW_SUCCESS;

  for (i = 0; i < grp->grp_index[grp->n_grp]; i++) {
    int eid    = grp->grp_item[2 * i];
    int *sid_p = &(grp->grp_item[2 * i + 1]);
    const int *sid_r2h;
    int ierror;

    HECMW_assert(0 < eid && eid <= mesh->n_elem_gross);

    sid_r2h = get_sid_r2h(mesh->elem_type[eid - 1], &ierror);
    if (ierror) return HECMW_ERROR;
    if (sid_r2h) *sid_p = sid_r2h[*sid_p];
  }
  return HECMW_SUCCESS;
}

/*============================================================================*/
/*                                                                            */
/* for models with multiple element types                                     */
/*                                                                            */
/*============================================================================*/
#define NDIV_SEG 2
#define NDIV_SURF 4
#define NDIV_BODY 8
#define NDIV_OTHER 1
#define NDIV_MAX 8

static int get_elem_ndiv(int etype, int *ierror) {
  int ndiv;
  *ierror = 0;
  if (HECMW_is_etype_truss(etype)) {
    ndiv = NDIV_OTHER;
  } else if (HECMW_is_etype_rod(etype) || HECMW_is_etype_beam(etype)) {
    ndiv = NDIV_SEG;
  } else if (HECMW_is_etype_surface(etype) || HECMW_is_etype_shell(etype)
             || HECMW_is_etype_patch(etype)) {
    ndiv = NDIV_SURF;
  } else if (HECMW_is_etype_solid(etype)) {
    ndiv = NDIV_BODY;
  } else if (HECMW_is_etype_link(etype)) {
    ndiv = NDIV_OTHER;
  } else {
    HECMW_log(HECMW_LOG_ERROR,
              "Element type %d not supported for refinement.\n", etype);
    *ierror = 1;
    ndiv    = NDIV_OTHER;
  }
  return ndiv;
}

/*============================================================================*/
/*                                                                            */
/*  prepare, clear and terminate REVOCAP_Refiner                              */
/*                                                                            */
/*============================================================================*/
static void register_node(struct hecmwST_local_mesh *mesh) {
  HECMW_log(HECMW_LOG_DEBUG, "Original Node Count = %d\n", mesh->n_node_gross);

  rcapSetNode64(mesh->n_node_gross, mesh->node, mesh->global_node_ID, NULL);

  HECMW_assert(rcapGetNodeCount() == mesh->n_node_gross);
}

static void register_node_groups(struct hecmwST_node_grp *grp) {
  int i;
  char rcap_name[80];

  HECMW_assert(strcmp(grp->grp_name[0], "ALL") == 0);

  for (i = 1; i < grp->n_grp; i++) {
    int start  = grp->grp_index[i];
    int num    = grp->grp_index[i + 1] - start;
    int *array = grp->grp_item + start;
    sprintf(rcap_name, "NG_%s", grp->grp_name[i]);
    rcapAppendNodeGroup(rcap_name, num, array);
  }
}

static int register_surf_groups(struct hecmwST_local_mesh *mesh) {
  struct hecmwST_surf_grp *grp = mesh->surf_group;
  int i, j;
  char rcap_name[80];

  struct hecmw_varray_int other;

  for (i = 0; i < grp->n_grp; i++) {
    int start  = grp->grp_index[i];
    int num    = grp->grp_index[i + 1] - start;
    int *array = grp->grp_item + start * 2;

    if (HECMW_varray_int_init(&other) != HECMW_SUCCESS) return HECMW_ERROR;

    for (j = 0; j < num; j++) {
      int elem  = array[j * 2];
      int surf  = array[j * 2 + 1];
      int etype = mesh->elem_type[elem - 1];

      /* ignore shell/patch surface */
      if (HECMW_is_etype_shell(etype) || HECMW_is_etype_patch(etype)) continue;

      if (HECMW_varray_int_append(&other, elem) != HECMW_SUCCESS)
        return HECMW_ERROR;
      if (HECMW_varray_int_append(&other, surf) != HECMW_SUCCESS)
        return HECMW_ERROR;
    }

    sprintf(rcap_name, "SG_%s", grp->grp_name[i]);
    num   = HECMW_varray_int_nval(&other) / 2;
    array = HECMW_varray_int_get_v(&other);
    rcapAppendFaceGroup(rcap_name, num, array);

    HECMW_varray_int_finalize(&other);
  }
  return HECMW_SUCCESS;
}

static int prepare_refiner(struct hecmwST_local_mesh *mesh,
                           const char *cad_filename,
                           const char *part_filename) {
  HECMW_log(HECMW_LOG_DEBUG, "Started preparing refiner...\n");

  if (elem_node_item_hecmw2rcap(mesh) != HECMW_SUCCESS) return HECMW_ERROR;
  if (surf_ID_hecmw2rcap(mesh) != HECMW_SUCCESS) return HECMW_ERROR;

  rcapInitRefiner(1, 1);
  if (cad_filename != NULL) {
    rcapSetCADFilename(cad_filename);
  }
  if (part_filename != NULL) {
    rcapSetPartitionFilename(part_filename);
  }
  register_node(mesh);

  /* node groups are refined by REVOCAP_Refiner */
  register_node_groups(mesh->node_group);

  /* element groups are refined by myself */

  /* surface groups are refined by REVOCAP_Refiner except for shell/patch surface */
  if (register_surf_groups(mesh) != HECMW_SUCCESS) return HECMW_ERROR;

  HECMW_log(HECMW_LOG_DEBUG, "Finished preparing refiner.\n");
  return HECMW_SUCCESS;
}

static void clear_refiner(void) { rcapClearRefiner(); }

static int terminate_refiner(struct hecmwST_local_mesh *mesh) {
  rcapTermRefiner();
  if (surf_ID_rcap2hecmw(mesh) != HECMW_SUCCESS) return HECMW_ERROR;
  if (elem_node_item_rcap2hecmw(mesh) != HECMW_SUCCESS) return HECMW_ERROR;
  return HECMW_SUCCESS;
}

/*============================================================================*/
/*                                                                            */
/*  call REVOCAP_Refiner                                                      */
/*                                                                            */
/*============================================================================*/
static int elem_type_hecmw2rcap(int etype) {
  switch (etype) {
    case HECMW_ETYPE_ROD1:
      return RCAP_SEGMENT;
    case HECMW_ETYPE_ROD2:
      return RCAP_SEGMENT2;
    case HECMW_ETYPE_TRI1:
      return RCAP_TRIANGLE;
    case HECMW_ETYPE_TRI2:
      return RCAP_TRIANGLE2;
    case HECMW_ETYPE_QUA1:
      return RCAP_QUAD;
    case HECMW_ETYPE_QUA2:
      return RCAP_QUAD2;
    case HECMW_ETYPE_ROD31:
      return RCAP_SEGMENT;
    case HECMW_ETYPE_TET1:
      return RCAP_TETRAHEDRON;
    case HECMW_ETYPE_TET1_4:
      return RCAP_TETRAHEDRON;
    case HECMW_ETYPE_TET2:
      return RCAP_TETRAHEDRON2;
    case HECMW_ETYPE_PRI1:
      return RCAP_WEDGE;
    case HECMW_ETYPE_PRI2:
      return RCAP_WEDGE2;
    case HECMW_ETYPE_PYR1:
      return RCAP_PYRAMID;
    case HECMW_ETYPE_PYR2:
      return RCAP_PYRAMID2;
    case HECMW_ETYPE_HEX1:
      return RCAP_HEXAHEDRON;
    case HECMW_ETYPE_HEX1_4:
      return RCAP_HEXAHEDRON;
    case HECMW_ETYPE_HEX2:
      return RCAP_HEXAHEDRON2;
    case HECMW_ETYPE_BEM1:
      return RCAP_SEGMENT;
    case HECMW_ETYPE_BEM2:
      return RCAP_SEGMENT2;
    case HECMW_ETYPE_BEM3:
      return RCAP_SEGMENT;
    case HECMW_ETYPE_SHT1:
      return RCAP_TRIANGLE;
    case HECMW_ETYPE_SHT2:
      return RCAP_TRIANGLE2;
    case HECMW_ETYPE_SHT6:
      return RCAP_TRIANGLE;
    case HECMW_ETYPE_SHQ1:
      return RCAP_QUAD;
    case HECMW_ETYPE_SHQ2:
      return RCAP_QUAD2;
    case HECMW_ETYPE_SHQ8:
      return RCAP_QUAD;
    case HECMW_ETYPE_PTT1:
      return RCAP_TRIANGLE;
    case HECMW_ETYPE_PTT2:
      return RCAP_TRIANGLE2;
    case HECMW_ETYPE_PTQ1:
      return RCAP_QUAD;
    case HECMW_ETYPE_PTQ2:
      return RCAP_QUAD2;
  }
  return RCAP_UNKNOWNTYPE;
}

static int reorder_enode_33struct(int etype, int n_elem,
                                  int *elem_node_item_ref) {
  int nn, nn_2, ndiv, nn_tot, ierror;
  int j, k, l;
  int tmp[32];
  int *enode1, *enode2;
  nn   = HECMW_get_max_node(etype);
  nn_2 = nn / 2;
  ndiv = get_elem_ndiv(etype, &ierror);
  if (ierror) return HECMW_ERROR;
  nn_tot = ndiv * nn;
  enode1 = elem_node_item_ref;
  enode2 = enode1 + nn_2 * ndiv;
  for (j = 0; j < n_elem; j++) {
    int *tmp_p = tmp;
    int *en1   = enode1;
    for (k = 0; k < ndiv; k++) {
      for (l = 0; l < nn_2; l++) {
        tmp_p[l]        = enode1[l];
        tmp_p[nn_2 + l] = enode2[l];
      }
      tmp_p += nn;
      enode1 += nn_2;
      enode2 += nn_2;
    }
    for (k = 0; k < nn_tot; k++) {
      en1[k] = tmp[k];
    }
    enode1 += nn_2 * ndiv;
    enode2 += nn_2 * ndiv;
  }
  return HECMW_SUCCESS;
}

static int refine_element(struct hecmwST_local_mesh *mesh,
                          struct hecmwST_local_mesh *ref_mesh) {
  int i, j;
  int n_elem_ref_tot;
  size_t n_enode_ref_tot;
  int *elem_node_index_ref;
  int *elem_node_item_ref;

  HECMW_log(HECMW_LOG_DEBUG, "Original Element Count = %d\n", mesh->n_elem_gross);

  /*
   * count num of elems after refinement
   */
  n_elem_ref_tot  = 0;
  n_enode_ref_tot = 0;
  for (i = 0; i < mesh->n_elem_type; i++) {
    int etype, istart, iend, n_elem, etype_rcap, nn;
    size_t n_elem_ref;
    int *elem_node_item;
    int ndiv, ierror;
    etype = mesh->elem_type_item[i];
    ndiv  = get_elem_ndiv(etype, &ierror);
    if (ierror) return HECMW_ERROR;
    istart = mesh->elem_type_index[i];
    iend   = mesh->elem_type_index[i + 1];
    n_elem = iend - istart;
    if (HECMW_is_etype_link(etype) || HECMW_is_etype_truss(etype)) {
      n_elem_ref = n_elem;
    } else {
      int is33 = HECMW_is_etype_33struct(etype);
      if (is33) n_elem *= 2;
      etype_rcap = elem_type_hecmw2rcap(etype);
      HECMW_assert(etype_rcap != RCAP_UNKNOWNTYPE);
      elem_node_item = mesh->elem_node_item + mesh->elem_node_index[istart];
      n_elem_ref = rcapRefineElement(n_elem, etype_rcap, elem_node_item, NULL);
      if (is33) {
        n_elem_ref /= 2;
        n_elem /= 2;
      }
    }
    HECMW_assert(n_elem_ref == n_elem * ndiv);
    n_elem_ref_tot += n_elem_ref;
    nn = HECMW_get_max_node(etype);
    n_enode_ref_tot += n_elem_ref * nn;
  }
  ref_mesh->n_elem_gross = n_elem_ref_tot;

  HECMW_log(HECMW_LOG_DEBUG, "Refined Element Count will be %d\n", ref_mesh->n_elem_gross);

  /*
   * allocate memory
   */
  ref_mesh->elem_node_index =
      (int *)HECMW_malloc(sizeof(int) * (n_elem_ref_tot + 1));
  if (ref_mesh->elem_node_index == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  ref_mesh->elem_node_item = (int *)HECMW_calloc(n_enode_ref_tot, sizeof(int));
  if (ref_mesh->elem_node_item == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }

  HECMW_assert(rcapGetNodeCount() == mesh->n_node_gross);

  /*
   * perform actual refinement
   */
  ref_mesh->elem_node_index[0] = 0;
  elem_node_index_ref          = ref_mesh->elem_node_index;
  elem_node_item_ref           = ref_mesh->elem_node_item;
  n_elem_ref_tot               = 0;
  for (i = 0; i < mesh->n_elem_type; i++) {
    int etype, istart, iend, n_elem, etype_rcap, n_elem_ref, nn, jstart;
    int *elem_node_item;
    etype  = mesh->elem_type_item[i];
    istart = mesh->elem_type_index[i];
    iend   = mesh->elem_type_index[i + 1];
    n_elem = iend - istart;
    nn     = HECMW_get_max_node(etype);
    if (HECMW_is_etype_link(etype) || HECMW_is_etype_truss(etype)) {
      n_elem_ref = n_elem;
      for (j = 0; j < n_elem * nn; j++) {
        /* elem_node_item_ref[j] = elem_node_item[j]; */
        jstart = mesh->elem_node_index[istart];
        /* printf("aa %d %d\n",jstart,j); */
        elem_node_item_ref[j] = mesh->elem_node_item[jstart + j];
      }
    } else {
      int is33 = HECMW_is_etype_33struct(etype);
      if (is33) n_elem *= 2;
      etype_rcap     = elem_type_hecmw2rcap(etype);
      elem_node_item = mesh->elem_node_item + mesh->elem_node_index[istart];
      n_elem_ref     = rcapRefineElement(n_elem, etype_rcap, elem_node_item,
                                     elem_node_item_ref);
      if (is33) {
        n_elem /= 2;
        n_elem_ref /= 2;
        reorder_enode_33struct(etype, n_elem, elem_node_item_ref);
      }
    }
    for (j = 0; j < n_elem_ref; j++) {
      elem_node_index_ref[j + 1] = elem_node_index_ref[j] + nn;
      if (elem_node_index_ref[j + 1] < elem_node_index_ref[j]) {
        HECMW_log(HECMW_LOG_ERROR,
                  "Integer overflow detected while creating elem_node_index\n");
        return HECMW_ERROR;
      }
    }
    elem_node_index_ref += n_elem_ref;
    elem_node_item_ref += nn * n_elem_ref;
    n_elem_ref_tot += n_elem_ref;
  }

  HECMW_log(HECMW_LOG_DEBUG, "Refined Element Count = %d\n", n_elem_ref_tot);

  HECMW_assert(n_elem_ref_tot == ref_mesh->n_elem_gross);

  rcapCommit();
  return HECMW_SUCCESS;
}

static int refine_node(struct hecmwST_local_mesh *mesh,
                       struct hecmwST_local_mesh *ref_mesh) {
  /* n_node_gross */
  ref_mesh->n_node_gross = rcapGetNodeCount();

  HECMW_log(HECMW_LOG_DEBUG, "Refined Node Count = %d\n", ref_mesh->n_node_gross);

  HECMW_assert(ref_mesh->n_node_gross > mesh->n_node_gross);

  /* node */
  ref_mesh->node =
      (double *)HECMW_malloc(sizeof(double) * ref_mesh->n_node_gross * 3);
  if (ref_mesh->node == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  rcapGetNodeSeq64(ref_mesh->n_node_gross, 1, ref_mesh->node);

  return HECMW_SUCCESS;
}

static int refine_node_group_info(struct hecmwST_node_grp *grp,
                                  struct hecmwST_node_grp *ref_grp,
                                  int ref_n_node_gross) {
  int i;
  char rcap_name[80];

  ref_grp->n_grp = grp->n_grp;

  /* grp_name (COPY) */
  ref_grp->grp_name = (char **)HECMW_malloc(sizeof(char *) * grp->n_grp);
  if (ref_grp->grp_name == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  /* ALL */
  ref_grp->grp_name[0] = HECMW_strdup("ALL");
  if (ref_grp->grp_name[0] == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  /* other groups */
  for (i = 1; i < grp->n_grp; i++) {
    ref_grp->grp_name[i] = HECMW_strdup(grp->grp_name[i]);
    if (ref_grp->grp_name[i] == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
  }

  /* grp_index */
  ref_grp->grp_index = (int *)HECMW_malloc(sizeof(int) * (grp->n_grp + 1));
  if (ref_grp->grp_index == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  ref_grp->grp_index[0] = 0;
  /* ALL */
  ref_grp->grp_index[1] = ref_n_node_gross;
  /* other groups */
  for (i = 1; i < grp->n_grp; i++) {
    sprintf(rcap_name, "NG_%s", grp->grp_name[i]);
    ref_grp->grp_index[i + 1] =
        ref_grp->grp_index[i] + rcapGetNodeGroupCount(rcap_name);
  }

  /* grp_item */
  if (ref_grp->grp_index[grp->n_grp] > 0) {
    ref_grp->grp_item =
        (int *)HECMW_calloc(ref_grp->grp_index[grp->n_grp], sizeof(int));
    if (ref_grp->grp_item == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    /* ALL */
    for (i = 0; i < ref_n_node_gross; i++) {
      ref_grp->grp_item[i] = i + 1;
    }
    /* other groups */
    for (i = 1; i < grp->n_grp; i++) {
      int start = ref_grp->grp_index[i];
      int num   = ref_grp->grp_index[i + 1] - start;
      if (num > 0) {
        int *array = ref_grp->grp_item + start;
        sprintf(rcap_name, "NG_%s", grp->grp_name[i]);
        rcapGetNodeGroup(rcap_name, num, array);
      }
    }
  } else {
    ref_grp->grp_item = NULL;
  }
  return HECMW_SUCCESS;
}

/*
 * static functions for refining element groups
 */

static int get_refined_element_count(const struct hecmwST_local_mesh *mesh,
                                     int n_elem, int *elem_list) {
  int i, eid, etype, ierror, ndiv;
  int n_elem_ref = 0;
  for (i = 0; i < n_elem; i++) {
    eid = elem_list[i];
    HECMW_assert(0 < eid && eid <= mesh->n_elem_gross);
    etype = mesh->elem_type[eid - 1];
    ndiv  = get_elem_ndiv(etype, &ierror);
    if (ierror) return -1;
    n_elem_ref += ndiv;
  }
  return n_elem_ref;
}

static int get_refined_elements(const struct hecmwST_local_mesh *mesh,
                                struct hecmwST_local_mesh *ref_mesh, int eid,
                                int *eid_ref) {
  int i, etype, ierror, ndiv;
  int etype_idx, eid_type;
  int offset;

  HECMW_assert(0 < eid);
  HECMW_assert(eid <= mesh->n_elem_gross);
  HECMW_assert(eid_ref);

  etype = mesh->elem_type[eid - 1];
  ndiv  = get_elem_ndiv(etype, &ierror);
  if (ierror) return -1;
  etype_idx = -1;
  eid_type  = -1;
  for (i = 0; i < mesh->n_elem_type; i++) {
    if (mesh->elem_type_item[i] == etype) {
      etype_idx = i;
      HECMW_assert(mesh->elem_type_index[i] < eid);
      HECMW_assert(eid <= mesh->elem_type_index[i + 1]);
      eid_type = eid - mesh->elem_type_index[i];
      break;
    }
  }
  HECMW_assert(etype_idx >= 0);
  HECMW_assert(eid_type > 0);

  offset = ref_mesh->elem_type_index[etype_idx] + (eid_type - 1) * ndiv + 1;
  for (i = 0; i < ndiv; i++) {
    eid_ref[i] = offset + i;
  }
  return ndiv;
}

static int get_refined_elem_list(const struct hecmwST_local_mesh *mesh,
                                 struct hecmwST_local_mesh *ref_mesh,
                                 int n_elem, int *elem_list,
                                 int *elem_list_ref) {
  int i, j;
  int eid_ref[NDIV_MAX];
  int *el_ref    = elem_list_ref;
  int n_elem_ref = 0;
  for (i = 0; i < n_elem; i++) {
    int ndiv = get_refined_elements(mesh, ref_mesh, elem_list[i], eid_ref);
    if (ndiv < 0) return -1;
    for (j = 0; j < ndiv; j++) {
      el_ref[j] = eid_ref[j];
    }
    el_ref += ndiv;
    n_elem_ref += ndiv;
  }
  return n_elem_ref;
}

/*
 * refinement of element groups is done by myself
 */

static int refine_elem_group_info(struct hecmwST_local_mesh *mesh,
                                  struct hecmwST_local_mesh *ref_mesh) {
  int i;
  struct hecmwST_elem_grp *grp     = mesh->elem_group;
  struct hecmwST_elem_grp *ref_grp = ref_mesh->elem_group;
  int ref_n_elem_gross             = ref_mesh->n_elem_gross;

  ref_grp->n_grp = grp->n_grp;

  /* grp_name (COPY) */
  ref_grp->grp_name = (char **)HECMW_malloc(sizeof(char *) * grp->n_grp);
  if (ref_grp->grp_name == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  /* ALL */
  HECMW_assert(strcmp(grp->grp_name[0], "ALL") == 0);
  ref_grp->grp_name[0] = HECMW_strdup("ALL");
  if (ref_grp->grp_name[0] == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  /* other groups */
  for (i = 1; i < grp->n_grp; i++) {
    ref_grp->grp_name[i] = HECMW_strdup(grp->grp_name[i]);
    if (ref_grp->grp_name[i] == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
  }

  /* grp_index */
  ref_grp->grp_index = (int *)HECMW_malloc(sizeof(int) * (grp->n_grp + 1));
  if (ref_grp->grp_index == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  ref_grp->grp_index[0] = 0;
  /* ALL */
  ref_grp->grp_index[1] = ref_n_elem_gross;
  /* other groups */
  for (i = 1; i < grp->n_grp; i++) {
    int start, n_elem, n_elem_ref;
    int *elem_list;
    start      = grp->grp_index[i];
    n_elem     = grp->grp_index[i + 1] - start;
    elem_list  = grp->grp_item + start;
    n_elem_ref = get_refined_element_count(mesh, n_elem, elem_list);
    if (n_elem_ref < 0) return HECMW_ERROR;
    ref_grp->grp_index[i + 1] = ref_grp->grp_index[i] + n_elem_ref;
  }

  /* grp_item */
  if (ref_grp->grp_index[grp->n_grp] > 0) {
    ref_grp->grp_item =
        (int *)HECMW_calloc(ref_grp->grp_index[grp->n_grp], sizeof(int));
    if (ref_grp->grp_item == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    /* ALL */
    for (i = 0; i < ref_n_elem_gross; i++) {
      ref_grp->grp_item[i] = i + 1;
    }
    /* other groups */
    for (i = 1; i < grp->n_grp; i++) {
      int start, start_ref, n_elem, n_elem_ref, ret;
      int *elem_list, *elem_list_ref;
      start_ref  = ref_grp->grp_index[i];
      n_elem_ref = ref_grp->grp_index[i + 1] - start_ref;
      HECMW_assert(n_elem_ref >= 0);
      if (n_elem_ref == 0) continue;
      start         = grp->grp_index[i];
      n_elem        = grp->grp_index[i + 1] - start;
      elem_list     = grp->grp_item + start;
      elem_list_ref = ref_grp->grp_item + start_ref;
      ret           = get_refined_elem_list(mesh, ref_mesh, n_elem, elem_list,
                                  elem_list_ref);
      HECMW_assert(ret == n_elem_ref);
    }
  } else {
    ref_grp->grp_item = NULL;
  }
  return HECMW_SUCCESS;
}

/*
 * refinement of surface groups except for shell/patch surface is done
 * by REVOCAP_Refiner
 * refinement of shell/patch surfaces is done by myself
 */

static int refine_surf_group_info(struct hecmwST_local_mesh *mesh,
                                  struct hecmwST_local_mesh *ref_mesh) {
  struct hecmwST_surf_grp *grp     = mesh->surf_group;
  struct hecmwST_surf_grp *ref_grp = ref_mesh->surf_group;
  int i, j;
  char rcap_name[80];

  ref_grp->n_grp = grp->n_grp;

  /* grp_name (COPY) */
  ref_grp->grp_name = (char **)HECMW_malloc(sizeof(char *) * grp->n_grp);
  if (ref_grp->grp_name == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < grp->n_grp; i++) {
    ref_grp->grp_name[i] = HECMW_strdup(grp->grp_name[i]);
    if (ref_grp->grp_name[i] == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
  }

  /* grp_index */
  ref_grp->grp_index = (int *)HECMW_malloc(sizeof(int) * (grp->n_grp + 1));
  if (ref_grp->grp_index == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  ref_grp->grp_index[0] = 0;
  for (i = 0; i < grp->n_grp; i++) {
    int start  = grp->grp_index[i];
    int num    = grp->grp_index[i + 1] - start;
    int *array = grp->grp_item + start * 2;
    int num_sh, num_other;

    /* count surfaces except for shell/patch */
    sprintf(rcap_name, "SG_%s", grp->grp_name[i]);
    num_other = rcapGetFaceGroupCount(rcap_name);

    /* count shell/patch surfaces */
    num_sh = 0;
    for (j = 0; j < num; j++) {
      int elem  = array[j * 2];
      int etype = mesh->elem_type[elem - 1];
      if (!(HECMW_is_etype_shell(etype) || HECMW_is_etype_patch(etype))) continue;
      num_sh += 1;
    }
    num_sh *= NDIV_SURF;

    ref_grp->grp_index[i + 1] = ref_grp->grp_index[i] + num_other + num_sh;
  }
  HECMW_assert(ref_grp->grp_index[grp->n_grp] >= 0);

  /* grp_item */
  if (ref_grp->grp_index[grp->n_grp] == 0) {
    ref_grp->grp_item = NULL;
    return HECMW_SUCCESS;
  }
  ref_grp->grp_item =
      (int *)HECMW_calloc(ref_grp->grp_index[grp->n_grp] * 2, sizeof(int));
  if (ref_grp->grp_item == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < grp->n_grp; i++) {
    struct hecmw_varray_int sh1, sh2;
    int start, num_tot, n_elem, ret;
    int start_ref, num_ref, num_tot_ref, n_elem_ref;
    int *array, *array_ref, *elem_list, *elem_list_ref;

    start_ref   = ref_grp->grp_index[i];
    num_tot_ref = ref_grp->grp_index[i + 1] - start_ref;
    if (num_tot_ref == 0) continue;
    array_ref = ref_grp->grp_item + start_ref * 2;

    /* get surfaces from REVOCAP_Refiner */
    sprintf(rcap_name, "SG_%s", grp->grp_name[i]);
    num_ref = rcapGetFaceGroupCount(rcap_name);
    rcapGetFaceGroup(rcap_name, num_ref, array_ref);

    num_tot_ref -= num_ref;
    if (num_tot_ref == 0) continue;
    array_ref += num_ref * 2;

    /*
     * make refined shell/patch surfaces by myself
     */

    if (HECMW_varray_int_init(&sh1) != HECMW_SUCCESS) return HECMW_ERROR;
    if (HECMW_varray_int_init(&sh2) != HECMW_SUCCESS) return HECMW_ERROR;

    /* collect shell/patch surface */
    start   = grp->grp_index[i];
    num_tot = grp->grp_index[i + 1] - start;
    array   = grp->grp_item + start * 2;
    for (j = 0; j < num_tot; j++) {
      int elem, etype, surf;
      elem  = array[j * 2];
      surf  = array[j * 2 + 1];
      etype = mesh->elem_type[elem - 1];
      if (!(HECMW_is_etype_shell(etype) || HECMW_is_etype_patch(etype))) continue;
      if (surf == 1) {
        if (HECMW_varray_int_append(&sh1, elem) != HECMW_SUCCESS)
          return HECMW_ERROR;
      } else {
        HECMW_assert(surf == 2);
        if (HECMW_varray_int_append(&sh2, elem) != HECMW_SUCCESS)
          return HECMW_ERROR;
      }
    }

    /* 1st surface of shells/patches */
    n_elem        = HECMW_varray_int_nval(&sh1);
    if (n_elem > 0) {
      elem_list     = HECMW_varray_int_get_v(&sh1);
      n_elem_ref    = n_elem * NDIV_SURF;
      elem_list_ref = (int *)HECMW_calloc(n_elem_ref, sizeof(int));
      if (elem_list_ref == NULL) {
        HECMW_set_error(errno, "");
        return HECMW_ERROR;
      }
      ret = get_refined_elem_list(mesh, ref_mesh, n_elem, elem_list, elem_list_ref);
      HECMW_assert(ret == n_elem_ref);
      for (j = 0; j < n_elem_ref; j++) {
        array_ref[j * 2]     = elem_list_ref[j];
        array_ref[j * 2 + 1] = 1;
      }
      HECMW_free(elem_list_ref);

      num_tot_ref -= n_elem_ref;
      array_ref += n_elem_ref * 2;
    }

    /* 2nd surface of shells */
    n_elem        = HECMW_varray_int_nval(&sh2);
    if (n_elem > 0) {
      elem_list     = HECMW_varray_int_get_v(&sh2);
      n_elem_ref    = n_elem * NDIV_SURF;
      elem_list_ref = (int *)HECMW_calloc(n_elem_ref, sizeof(int));
      if (elem_list_ref == NULL) {
        HECMW_set_error(errno, "");
        return HECMW_ERROR;
      }
      ret = get_refined_elem_list(mesh, ref_mesh, n_elem, elem_list, elem_list_ref);
      HECMW_assert(ret == n_elem_ref);
      for (j = 0; j < n_elem_ref; j++) {
        array_ref[j * 2]     = elem_list_ref[j];
        array_ref[j * 2 + 1] = 2;
      }
      HECMW_free(elem_list_ref);

      num_tot_ref -= n_elem_ref;
    }

    HECMW_varray_int_finalize(&sh1);
    HECMW_varray_int_finalize(&sh2);

    HECMW_assert(num_tot_ref == 0);
  }
  return HECMW_SUCCESS;
}

static int call_refiner(struct hecmwST_local_mesh *mesh,
                        struct hecmwST_local_mesh *ref_mesh) {
  HECMW_log(HECMW_LOG_DEBUG, "Started calling refiner...\n");

  if (refine_element(mesh, ref_mesh) != HECMW_SUCCESS) {
    return HECMW_ERROR;
  }
  if (refine_node(mesh, ref_mesh) != HECMW_SUCCESS) {
    return HECMW_ERROR;
  }
  if (refine_node_group_info(mesh->node_group, ref_mesh->node_group,
                             ref_mesh->n_node_gross) != HECMW_SUCCESS) {
    return HECMW_ERROR;
  }
  if (refine_elem_group_info(mesh, ref_mesh) != HECMW_SUCCESS) {
    return HECMW_ERROR;
  }
  if (refine_surf_group_info(mesh, ref_mesh) != HECMW_SUCCESS) {
    return HECMW_ERROR;
  }
  HECMW_log(HECMW_LOG_DEBUG, "Finished calling refiner.\n");
  return HECMW_SUCCESS;
}

/*============================================================================*/
/*                                                                            */
/*  rebuild information                                                       */
/*                                                                            */
/*============================================================================*/
static int rebuild_elem_ID(const struct hecmwST_local_mesh *mesh,
                           struct hecmwST_local_mesh *ref_mesh) {
  int i, j, k;
  int count;
  int *elem_ID_org;
  int *elem_ID_ref;

  ref_mesh->elem_ID =
      (int *)HECMW_malloc(sizeof(int) * ref_mesh->n_elem_gross * 2);
  if (ref_mesh->elem_ID == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  count       = 0;
  elem_ID_org = mesh->elem_ID;
  elem_ID_ref = ref_mesh->elem_ID;
  for (i = 0; i < mesh->n_elem_type; i++) {
    int etype, istart, iend, n_elem, ierror, ndiv;
    etype  = mesh->elem_type_item[i];
    istart = mesh->elem_type_index[i];
    iend   = mesh->elem_type_index[i + 1];
    n_elem = iend - istart;
    ndiv   = get_elem_ndiv(etype, &ierror);
    if (ierror) return HECMW_ERROR;

    for (j = 0; j < n_elem; j++) {
      int rank = elem_ID_org[1];
      for (k = 0; k < ndiv; k++) {
        count++;
        elem_ID_ref[0] = -count; /* TEMPORARY */
        elem_ID_ref[1] = rank;   /* to be corrected later */
        elem_ID_ref += 2;
      }
      elem_ID_org += 2;
    }
  }
  return HECMW_SUCCESS;
}

static int rebuild_global_elem_ID(const struct hecmwST_local_mesh *mesh,
                                  struct hecmwST_local_mesh *ref_mesh) {
  int i, j, k;
  int min_ID;
  int *global_elem_ID_org;
  int *global_elem_ID_ref;

  ref_mesh->global_elem_ID =
      (int *)HECMW_malloc(sizeof(int) * ref_mesh->n_elem_gross);
  if (ref_mesh->global_elem_ID == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  min_ID = 0;
  global_elem_ID_org = mesh->global_elem_ID;
  for (i = 0; i < mesh->n_elem_gross; i++) {
    if (global_elem_ID_org[i] < min_ID)
      min_ID = global_elem_ID_org[i];
  }
  global_elem_ID_ref = ref_mesh->global_elem_ID;
  for (i = 0; i < mesh->n_elem_type; i++) {
    int etype, istart, iend, n_elem, ierror, ndiv;
    etype  = mesh->elem_type_item[i];
    istart = mesh->elem_type_index[i];
    iend   = mesh->elem_type_index[i + 1];
    n_elem = iend - istart;
    ndiv   = get_elem_ndiv(etype, &ierror);
    if (ierror) return HECMW_ERROR;

    for (j = 0; j < n_elem; j++) {
      global_elem_ID_ref[0] = global_elem_ID_org[0];
      for (k = 1; k < ndiv; k++) {
        min_ID--;
        global_elem_ID_ref[k] = min_ID;
      }
      global_elem_ID_org += 1;
      global_elem_ID_ref += ndiv;
    }
  }
  return HECMW_SUCCESS;
}

static int rebuild_elem_type(const struct hecmwST_local_mesh *mesh,
                             struct hecmwST_local_mesh *ref_mesh) {
  int i, j, k;
  int *elem_type_ref;

  ref_mesh->elem_type =
      (int *)HECMW_malloc(sizeof(int) * ref_mesh->n_elem_gross);
  if (ref_mesh->elem_type == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  elem_type_ref = ref_mesh->elem_type;
  for (i = 0; i < mesh->n_elem_type; i++) {
    int etype, istart, iend, n_elem, ierror, ndiv;
    etype  = mesh->elem_type_item[i];
    istart = mesh->elem_type_index[i];
    iend   = mesh->elem_type_index[i + 1];
    n_elem = iend - istart;
    ndiv   = get_elem_ndiv(etype, &ierror);
    if (ierror) return HECMW_ERROR;

    for (j = 0; j < n_elem; j++) {
      for (k = 0; k < ndiv; k++) {
        elem_type_ref[k] = etype;
      }
      elem_type_ref += ndiv;
    }
  }
  return HECMW_SUCCESS;
}

static int rebuild_section_ID(const struct hecmwST_local_mesh *mesh,
                              struct hecmwST_local_mesh *ref_mesh) {
  int i, j, k;
  int *section_ID_org;
  int *section_ID_ref;
  ref_mesh->section_ID =
      (int *)HECMW_malloc(sizeof(int) * ref_mesh->n_elem_gross);
  if (ref_mesh->section_ID == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  section_ID_org = mesh->section_ID;
  section_ID_ref = ref_mesh->section_ID;
  for (i = 0; i < mesh->n_elem_type; i++) {
    int etype, istart, iend, n_elem, ierror, ndiv;
    etype  = mesh->elem_type_item[i];
    istart = mesh->elem_type_index[i];
    iend   = mesh->elem_type_index[i + 1];
    n_elem = iend - istart;
    ndiv   = get_elem_ndiv(etype, &ierror);
    if (ierror) return HECMW_ERROR;

    for (j = 0; j < n_elem; j++) {
      int sect_id = section_ID_org[0];
      for (k = 0; k < ndiv; k++) {
        section_ID_ref[k] = sect_id;
      }
      section_ID_org += 1;
      section_ID_ref += ndiv;
    }
  }
  return HECMW_SUCCESS;
}

static int rebuild_elem_mat_ID_index(const struct hecmwST_local_mesh *mesh,
                                     struct hecmwST_local_mesh *ref_mesh) {
  int i, j, k;
  int *elem_mat_ID_index_org;
  int *elem_mat_ID_index_ref;

  ref_mesh->elem_mat_ID_index =
      (int *)HECMW_malloc(sizeof(int) * (ref_mesh->n_elem_gross + 1));
  if (ref_mesh->elem_mat_ID_index == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  ref_mesh->elem_mat_ID_index[0] = 0;
  elem_mat_ID_index_org          = mesh->elem_mat_ID_index;
  elem_mat_ID_index_ref          = ref_mesh->elem_mat_ID_index;
  for (i = 0; i < mesh->n_elem_type; i++) {
    int etype, istart, iend, n_elem, ierror, ndiv;
    etype  = mesh->elem_type_item[i];
    istart = mesh->elem_type_index[i];
    iend   = mesh->elem_type_index[i + 1];
    n_elem = iend - istart;
    ndiv   = get_elem_ndiv(etype, &ierror);
    if (ierror) return HECMW_ERROR;

    for (j = 0; j < n_elem; j++) {
      int n = elem_mat_ID_index_org[1] - elem_mat_ID_index_org[0];
      for (k = 0; k < ndiv; k++) {
        elem_mat_ID_index_ref[k + 1] = elem_mat_ID_index_ref[k] + n;
      }
      elem_mat_ID_index_org += 1;
      elem_mat_ID_index_ref += ndiv;
    }
  }
  return HECMW_SUCCESS;
}

static int rebuild_elem_mat_ID_item(const struct hecmwST_local_mesh *mesh,
                                    struct hecmwST_local_mesh *ref_mesh) {
  int i, j, k, l;
  int *elem_mat_ID_index_org;
  int *elem_mat_ID_item_org;
  int *elem_mat_ID_item_ref;

  ref_mesh->elem_mat_ID_item = (int *)HECMW_malloc(
      sizeof(int) * ref_mesh->elem_mat_ID_index[ref_mesh->n_elem_gross]);
  if (ref_mesh->elem_mat_ID_item == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  elem_mat_ID_index_org = mesh->elem_mat_ID_index;
  elem_mat_ID_item_org  = mesh->elem_mat_ID_item;
  elem_mat_ID_item_ref  = ref_mesh->elem_mat_ID_item;
  for (i = 0; i < mesh->n_elem_type; i++) {
    int etype, istart, iend, n_elem, ierror, ndiv;
    etype  = mesh->elem_type_item[i];
    istart = mesh->elem_type_index[i];
    iend   = mesh->elem_type_index[i + 1];
    n_elem = iend - istart;
    ndiv   = get_elem_ndiv(etype, &ierror);
    if (ierror) return HECMW_ERROR;

    for (j = 0; j < n_elem; j++) {
      int n = elem_mat_ID_index_org[1] - elem_mat_ID_index_org[0];
      for (k = 0; k < ndiv; k++) {
        for (l = 0; l < n; l++) {
          elem_mat_ID_item_ref[l] = elem_mat_ID_item_org[l];
        }
        elem_mat_ID_item_ref += n;
      }
      elem_mat_ID_index_org += 1;
      elem_mat_ID_item_org += n;
    }
  }
  return HECMW_SUCCESS;
}

static int rebuild_elem_info(const struct hecmwST_local_mesh *mesh,
                             struct hecmwST_local_mesh *ref_mesh) {
  HECMW_log(HECMW_LOG_DEBUG, "Started rebuilding element info...\n");

  if (rebuild_elem_ID(mesh, ref_mesh) != HECMW_SUCCESS) return HECMW_ERROR;
  if (rebuild_global_elem_ID(mesh, ref_mesh) != HECMW_SUCCESS)
    return HECMW_ERROR;
  if (rebuild_elem_type(mesh, ref_mesh) != HECMW_SUCCESS) return HECMW_ERROR;
  if (rebuild_section_ID(mesh, ref_mesh) != HECMW_SUCCESS) return HECMW_ERROR;
  if (rebuild_elem_mat_ID_index(mesh, ref_mesh) != HECMW_SUCCESS)
    return HECMW_ERROR;
  if (rebuild_elem_mat_ID_item(mesh, ref_mesh) != HECMW_SUCCESS)
    return HECMW_ERROR;

  ref_mesh->n_elem_mat_ID = ref_mesh->elem_mat_ID_index[ref_mesh->n_elem_gross];

  HECMW_log(HECMW_LOG_DEBUG, "Finished rebuilding element info.\n");
  return HECMW_SUCCESS;
}

static int rebuild_node_info(const struct hecmwST_local_mesh *mesh,
                             struct hecmwST_local_mesh *ref_mesh) {
  int i;
  int count, min_ID;

  HECMW_log(HECMW_LOG_DEBUG, "Started rebuilding node info...\n");

  /* allocate node_ID */
  ref_mesh->node_ID =
      (int *)HECMW_malloc(sizeof(int) * ref_mesh->n_node_gross * 2);
  if (ref_mesh->node_ID == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }

  /* node_ID */
  for (i = 0; i < mesh->n_node_gross; i++) {
    ref_mesh->node_ID[2 * i]     = mesh->node_ID[2 * i];
    ref_mesh->node_ID[2 * i + 1] = mesh->node_ID[2 * i + 1];
  }
  count = 0;
  for (i = mesh->n_node_gross; i < ref_mesh->n_node_gross; i++) {
    count++;
    ref_mesh->node_ID[2 * i]     = -count;        /* TEMPORARY */
    ref_mesh->node_ID[2 * i + 1] = mesh->my_rank; /* to be corrected later */
  }

  /* global_node_ID */
  ref_mesh->global_node_ID =
      (int *)HECMW_malloc(sizeof(int) * ref_mesh->n_node_gross);
  if (ref_mesh->global_node_ID == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  min_ID = 0;
  for (i = 0; i < mesh->n_node_gross; i++) {
    ref_mesh->global_node_ID[i] = mesh->global_node_ID[i];
    if (mesh->global_node_ID[i] < min_ID)
      min_ID = mesh->global_node_ID[i];
  }
  for (i = mesh->n_node_gross; i < ref_mesh->n_node_gross; i++) {
    min_ID--;
    ref_mesh->global_node_ID[i] = min_ID;
  }

  if (mesh->hecmw_flag_initcon && mesh->n_node_gross > 0) {
    /* node_init_val_index */
    ref_mesh->node_init_val_index =
        (int *)HECMW_malloc(sizeof(int) * (ref_mesh->n_node_gross + 1));
    if (ref_mesh->node_init_val_index == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    for (i = 0; i < mesh->n_node_gross + 1; i++) {
      ref_mesh->node_init_val_index[i] = mesh->node_init_val_index[i];
    }
    for (i = mesh->n_node_gross + 1; i < ref_mesh->n_node_gross + 1; i++) {
      ref_mesh->node_init_val_index[i] =
          mesh->node_init_val_index[mesh->n_node_gross];
    }

    /* node_init_val_item */
    if (ref_mesh->node_init_val_index[ref_mesh->n_node_gross] == 0) {
      ref_mesh->node_init_val_item = NULL;
    } else {
      ref_mesh->node_init_val_item = (double *)HECMW_malloc(
          sizeof(double) * mesh->node_init_val_index[mesh->n_node_gross]);
      if (ref_mesh->node_init_val_item == NULL) {
        HECMW_set_error(errno, "");
        return HECMW_ERROR;
      }
      for (i = 0; i < mesh->node_init_val_index[mesh->n_node_gross]; i++) {
        ref_mesh->node_init_val_item[i] = mesh->node_init_val_item[i];
      }
    }
  }

  HECMW_log(HECMW_LOG_DEBUG, "Finished rebuilding node info.\n");
  return HECMW_SUCCESS;
}

/*
 * static functions called by rebuild_comm_tables
 */

static int get_refined_shared(const struct hecmwST_local_mesh *mesh,
                              struct hecmwST_local_mesh *ref_mesh,
                              struct hecmw_varray_int *shared) {
  int i;

  for (i = 0; i < mesh->n_neighbor_pe; i++) {
    int start, n_elem, n_elem_ref, ret;
    int *elem_list, *elem_list_ref;
    start      = mesh->shared_index[i];
    n_elem     = mesh->shared_index[i + 1] - start;
    elem_list  = mesh->shared_item + start;
    n_elem_ref = get_refined_element_count(mesh, n_elem, elem_list);
    if (n_elem_ref < 0) return HECMW_ERROR;
    HECMW_varray_int_resize(shared + i, n_elem_ref);
    elem_list_ref = HECMW_varray_int_get_v(shared + i);
    ret =
        get_refined_elem_list(mesh, ref_mesh, n_elem, elem_list, elem_list_ref);
    HECMW_assert(ret == n_elem_ref);
  }
  return HECMW_SUCCESS;
}

static int elem_type_rcap2hecmw(int etype) {
  switch (etype) {
    case RCAP_SEGMENT:
      return HECMW_ETYPE_ROD1;
    case RCAP_SEGMENT2:
      return HECMW_ETYPE_ROD2;
    case RCAP_TRIANGLE:
      return HECMW_ETYPE_TRI1;
    case RCAP_TRIANGLE2:
      return HECMW_ETYPE_TRI2;
    case RCAP_QUAD:
      return HECMW_ETYPE_QUA1;
    case RCAP_QUAD2:
      return HECMW_ETYPE_QUA2;
    case RCAP_TETRAHEDRON:
      return HECMW_ETYPE_TET1;
    case RCAP_TETRAHEDRON2:
      return HECMW_ETYPE_TET2;
    case RCAP_WEDGE:
      return HECMW_ETYPE_PRI1;
    case RCAP_WEDGE2:
      return HECMW_ETYPE_PRI2;
    case RCAP_PYRAMID:
      return HECMW_ETYPE_PYR1;
    case RCAP_PYRAMID2:
      return HECMW_ETYPE_PYR2;
    case RCAP_HEXAHEDRON:
      return HECMW_ETYPE_HEX1;
    case RCAP_HEXAHEDRON2:
      return HECMW_ETYPE_HEX2;
  }
  return HECMW_ERROR;
}

static int get_rcap_elem_max_node(int etype_rcap) {
  int etype_hecmw = elem_type_rcap2hecmw(etype_rcap);
  return HECMW_get_max_node(etype_hecmw);
}

static int determine_node_rank_of_purely_external_elems(
  const struct hecmwST_local_mesh *mesh,
  struct hecmwST_local_mesh *ref_mesh) {
  int eid, rank, nn, k, nid;
  int *elem_nodes;

  for (eid = 1; eid <= ref_mesh->n_elem_gross; eid++) {
    rank = ref_mesh->elem_ID[2 * (eid - 1) + 1];

    if (rank >= 0) continue; /* purely external elems have rank < 0 */

    nn = ref_mesh->elem_node_index[eid] - ref_mesh->elem_node_index[eid - 1];
    elem_nodes = ref_mesh->elem_node_item + ref_mesh->elem_node_index[eid - 1];
    for (k = 0; k < nn; k++) {
      nid = elem_nodes[k];
      if (nid <= mesh->n_node_gross) continue; /* skip old nodes */

      /* new nodes that have not been assigned new rank have default rank = my_rank */
      if (ref_mesh->node_ID[2 * (nid - 1) + 1] == ref_mesh->my_rank) {
        /* assign out-of-range rank = n_subdomain as purely external number = -rank-1 */
        ref_mesh->node_ID[2 * (nid - 1) + 1] = - ref_mesh->n_subdomain - 1;
      }
    }
  }
  return HECMW_SUCCESS;
}

static int determine_node_rank(const struct hecmwST_local_mesh *mesh,
                               const struct hecmw_varray_int *shared,
                               struct hecmwST_local_mesh *ref_mesh) {
  int i, j, k;
  struct hecmw_set_int boundary_nodes;
  int bnode;

  HECMW_log(HECMW_LOG_DEBUG, "Started determine_node_rank...\n");

  HECMW_set_int_init(&boundary_nodes);

  /* Collect nodes of the shared elements, except for the old nodes */
  for (i = 0; i < mesh->n_neighbor_pe; i++) {
    for (j = 0; j < HECMW_varray_int_nval(shared + i); j++) {
      int eid = HECMW_varray_int_get(shared + i, j);
      int nn =
          ref_mesh->elem_node_index[eid] - ref_mesh->elem_node_index[eid - 1];
      int *elem_nodes =
          ref_mesh->elem_node_item + ref_mesh->elem_node_index[eid - 1];

      for (k = 0; k < nn; k++) {
        if (elem_nodes[k] <= mesh->n_node_gross) continue;

        HECMW_set_int_add(&boundary_nodes, elem_nodes[k]);
      }
    }
  }
  HECMW_set_int_check_dup(&boundary_nodes);

  /* for each node: determine which domain it belongs to; save in node_ID */
  HECMW_set_int_iter_init(&boundary_nodes);
  while (HECMW_set_int_iter_next(&boundary_nodes, &bnode)) {
    int original[HECMW_MAX_NODE_MAX];
    int orig_type_rcap = rcapGetOriginal(bnode, original);
    int min_rank       = mesh->n_subdomain;
    int max_rank       = -1;
    int nn;

    HECMW_assert(orig_type_rcap != RCAP_UNKNOWNTYPE);

    nn = get_rcap_elem_max_node(orig_type_rcap);

    for (k = 0; k < nn; k++) {
      int rank;

      HECMW_assert(1 <= original[k] && original[k] < bnode);

      rank = ref_mesh->node_ID[2 * (original[k] - 1) + 1];

      /* for purely external nodes, rank (0, 1, ...) is stored as -rank-1 (-1,
       * -2, ...) */
      if (rank < 0) rank = -rank - 1;

      if (rank < min_rank) min_rank = rank;
      if (max_rank < rank) max_rank = rank;
    }
    HECMW_assert(0 <= min_rank && min_rank < mesh->n_subdomain);
    HECMW_assert(0 <= max_rank && max_rank < mesh->n_subdomain);

    /* min-rule for the first step, max-rule for the following steps */
    if (mesh->n_refine == 0) {
      ref_mesh->node_ID[2 * (bnode - 1) + 1] = min_rank;
    } else {
      ref_mesh->node_ID[2 * (bnode - 1) + 1] = max_rank;
    }
  }

  HECMW_set_int_finalize(&boundary_nodes);

  /* new nodes in purely external elems */
  if (determine_node_rank_of_purely_external_elems(mesh, ref_mesh) != HECMW_SUCCESS)
    return HECMW_ERROR;

  HECMW_log(HECMW_LOG_DEBUG, "Finished determine_node_rank.\n");
  return HECMW_SUCCESS;
}

static int *conv_index_ucd2hec(int etype)
{
  static int conv_index_ucd2hec_rod1[] = {0, 1};
  static int conv_index_ucd2hec_rod2[] = {0, -1, 2};
  static int conv_index_ucd2hec_tri1[] = {0, 1, 2};
  static int conv_index_ucd2hec_tri2[] = {0, 1, 2, -1, -1, -1};
  static int conv_index_ucd2hec_qua1[] = {0, 1, 2, 3};
  static int conv_index_ucd2hec_qua2[] = {0, 1, 2, 3, -1, -1, -1, -1};
  static int conv_index_ucd2hec_tet1[] = {0, 3, 2, 1};
  static int conv_index_ucd2hec_tet2[] = {0, 3, 2, 1, -1, -1, -1, -1, -1, -1};
  static int conv_index_ucd2hec_pri1[] = {3, 4, 5, 0, 1, 2};
  static int conv_index_ucd2hec_pri2[] = {3,  4,  5,  0,  1,  2,  -1, -1,
                                          -1, -1, -1, -1, -1, -1, -1};
  static int conv_index_ucd2hec_hex1[] = {4, 5, 6, 7, 0, 1, 2, 3};
  static int conv_index_ucd2hec_hex2[] = {4,  5,  6,  7,  0,  1,  2,  3,  -1, -1,
                                          -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
  static int conv_index_ucd2hec_pyr1[] = {4, 0, 1, 2, 3};
  static int conv_index_ucd2hec_pyr2[] = {4,  0,  1,  2,  3,  -1, -1,
                                          -1, -1, -1, -1, -1, -1};
  static int conv_index_ucd2hec_mst1[] = {0, 1, 2, 3};
  static int conv_index_ucd2hec_mst2[] = {0, 1, 2, 3, -1, -1, -1, -1};
  static int conv_index_ucd2hec_msq1[] = {0, 1, 2, 3, 4};
  static int conv_index_ucd2hec_msq2[] = {0, 1, 2, 3, 4, -1, -1, -1, -1};
  static int conv_index_ucd2hec_jtt1[] = {3, 4, 5, 0, 1, 2};
  static int conv_index_ucd2hec_jtt2[] = {3,  4,  5,  0,  1,  2,
                                          -1, -1, -1, -1, -1, -1};
  static int conv_index_ucd2hec_jtq1[] = {4, 5, 6, 7, 0, 1, 2, 3};
  static int conv_index_ucd2hec_jtq2[] = {4,  5,  6,  7,  0,  1,  2,  3,
                                          -1, -1, -1, -1, -1, -1, -1, -1};
#if 0
  static int conv_index_ucd2hec_bem1[] = {0, 1};
  static int conv_index_ucd2hec_bem2[] = {0, -1, 1};
#endif
  static int conv_index_ucd2hec_sht1[] = {0, 1, 2};
  static int conv_index_ucd2hec_sht2[] = {0, 1, 2, -1, -1, -1};
  static int conv_index_ucd2hec_shq1[] = {0, 1, 2, 3};
  static int conv_index_ucd2hec_shq2[] = {0, 1, 2, 3, -1, -1, -1, -1};
  static int conv_index_ucd2hec_ln[] = {0, 1};
  switch (etype) {
  case HECMW_ETYPE_ROD1:
    return conv_index_ucd2hec_rod1;
  case HECMW_ETYPE_ROD2:
    return conv_index_ucd2hec_rod2;
  case HECMW_ETYPE_TRI1:
    return conv_index_ucd2hec_tri1;
  case HECMW_ETYPE_TRI2:
    return conv_index_ucd2hec_tri2;
  case HECMW_ETYPE_QUA1:
    return conv_index_ucd2hec_qua1;
  case HECMW_ETYPE_QUA2:
    return conv_index_ucd2hec_qua2;
  case HECMW_ETYPE_TET1:
    return conv_index_ucd2hec_tet1;
  case HECMW_ETYPE_TET2:
    return conv_index_ucd2hec_tet2;
  case HECMW_ETYPE_PRI1:
    return conv_index_ucd2hec_pri1;
  case HECMW_ETYPE_PRI2:
    return conv_index_ucd2hec_pri2;
  case HECMW_ETYPE_HEX1:
    return conv_index_ucd2hec_hex1;
  case HECMW_ETYPE_HEX2:
    return conv_index_ucd2hec_hex2;
  case HECMW_ETYPE_PYR1:
    return conv_index_ucd2hec_pyr1;
  case HECMW_ETYPE_PYR2:
    return conv_index_ucd2hec_pyr2;
  case HECMW_ETYPE_MST1:
    return conv_index_ucd2hec_mst1;
  case HECMW_ETYPE_MST2:
    return conv_index_ucd2hec_mst2;
  case HECMW_ETYPE_MSQ1:
    return conv_index_ucd2hec_msq1;
  case HECMW_ETYPE_MSQ2:
    return conv_index_ucd2hec_msq2;
  case HECMW_ETYPE_JTT1:
    return conv_index_ucd2hec_jtt1;
  case HECMW_ETYPE_JTT2:
    return conv_index_ucd2hec_jtt2;
  case HECMW_ETYPE_JTQ1:
    return conv_index_ucd2hec_jtq1;
  case HECMW_ETYPE_JTQ2:
    return conv_index_ucd2hec_jtq2;
  case HECMW_ETYPE_SHT1:
    return conv_index_ucd2hec_sht1;
  case HECMW_ETYPE_SHT2:
    return conv_index_ucd2hec_sht2;
  case HECMW_ETYPE_SHQ1:
    return conv_index_ucd2hec_shq1;
  case HECMW_ETYPE_SHQ2:
    return conv_index_ucd2hec_shq2;
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
    return conv_index_ucd2hec_ln;
  default:
    return NULL;
  }
}

static int write_refined_shared_in_ucd(const struct hecmwST_local_mesh *mesh,
                                       const struct hecmwST_local_mesh *ref_mesh,
                                       const struct hecmw_varray_int *shared,
                                       const struct hecmw_varray_int *shared_nodes,
                                       int i)
{
  char fname[64];
  FILE *fp;
  int irank = mesh->neighbor_pe[i];
  int j, k, l, nid;
  int n_shared_elem = HECMW_varray_int_nval(shared + i);
  int n_shared_node = HECMW_varray_int_nval(shared_nodes + i);
  sprintf(fname, "shared_%d_%d.inp", mesh->my_rank, irank);
  HECMW_log(HECMW_LOG_DEBUG, "writing refined shared elements to %s\n", fname);
  fp = fopen(fname, "w");
  if (fp == NULL) {
    HECMW_log(HECMW_LOG_ERROR, "failed to open file %s\n", fname);
  } else {
    fprintf(fp, "%d %d 1 1 0\n", n_shared_node, n_shared_elem);
    for (j = 0; j < n_shared_node; j++) {
      nid = HECMW_varray_int_get(shared_nodes + i, j);
      fprintf(fp, "%d %e %e %e\n", j+1, ref_mesh->node[3*(nid-1)], ref_mesh->node[3*(nid-1)+1], ref_mesh->node[3*(nid-1)+2]);
    }
    for (j = 0; j < n_shared_elem; j++) {
      int eid = HECMW_varray_int_get(shared + i, j);
      int etype = ref_mesh->elem_type[eid-1];
      int *elem_nodes = ref_mesh->elem_node_item + ref_mesh->elem_node_index[eid-1];
      int nn = ref_mesh->elem_node_index[eid] - ref_mesh->elem_node_index[eid-1];
      int *index_ucd = conv_index_ucd2hec(etype);
      HECMW_assert(index_ucd);
      fprintf(fp, "%d 0 %s", j+1, HECMW_get_ucd_label(etype));
      for (k = 0; k < nn; k++) {
        if (index_ucd[k] < 0) continue;
        nid = elem_nodes[index_ucd[k]];
        for (l = 0; l < n_shared_node; l++) {
          if (HECMW_varray_int_get(shared_nodes + i, l) == nid) break;
        }
        fprintf(fp, " %d", l+1);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
}

static int check_node_rank(const struct hecmwST_local_mesh *mesh,
                           const struct hecmw_varray_int *shared,
                           const struct hecmwST_local_mesh *ref_mesh)
{
  int *n_shared_nodes;
  struct hecmw_varray_int *shared_nodes, *shared_ranks;
  HECMW_Request *requests;
  HECMW_Status *statuses;
  HECMW_Comm comm;
  int *send_buf;
  int recv_buf[4];
  int n_error = 0, n_error_g;
  int i, j, k, irank, tag, nid, rank;

  HECMW_log(HECMW_LOG_DEBUG, "Started check_node_rank...\n");

  /* alloc memory */
  n_shared_nodes = (int *)HECMW_malloc(sizeof(int) * mesh->n_neighbor_pe);
  if (n_shared_nodes == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }

  shared_nodes = (struct hecmw_varray_int *)HECMW_malloc(
    sizeof(struct hecmw_varray_int) * mesh->n_neighbor_pe);
  if (shared_nodes == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }

  shared_ranks = (struct hecmw_varray_int *)HECMW_malloc(
    sizeof(struct hecmw_varray_int) * mesh->n_neighbor_pe);
  if (shared_ranks == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }

  send_buf = (int *)HECMW_malloc(sizeof(int) * 4 * mesh->n_neighbor_pe);
  if (send_buf == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }

  comm = HECMW_comm_get_comm();

  requests = (HECMW_Request *)HECMW_malloc(sizeof(HECMW_Request) * mesh->n_neighbor_pe);
  if (requests == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }

  statuses = (HECMW_Status *)HECMW_malloc(sizeof(HECMW_Status) * mesh->n_neighbor_pe);
  if (statuses == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }

  /*
   * check number of refined shared nodes
   */

  /* Collect nodes of the shared elements, and send size */
  for (i = 0; i < mesh->n_neighbor_pe; i++) {
    /* numver of shared elements in original mesh */
    send_buf[4*i] = mesh->shared_index[i+1] - mesh->shared_index[i];

    /* number of refined shared elements */
    send_buf[4*i+1] = HECMW_varray_int_nval(shared + i);

    HECMW_varray_int_init(shared_nodes + i);
    for (j = 0; j < HECMW_varray_int_nval(shared + i); j++) {
      int eid = HECMW_varray_int_get(shared + i, j);
      int nn = ref_mesh->elem_node_index[eid] - ref_mesh->elem_node_index[eid - 1];
      int *elem_nodes = ref_mesh->elem_node_item + ref_mesh->elem_node_index[eid - 1];

      for (k = 0; k < nn; k++) {
#if 0  /* keep old nodes for plotting */
        if (elem_nodes[k] <= mesh->n_node_gross) continue; /* skip old nodes */
#endif

        HECMW_varray_int_append(shared_nodes + i, elem_nodes[k]);
      }
    }
    /* number of shared nodes before removing duplication */
    send_buf[4*i+2] = HECMW_varray_int_nval(shared_nodes + i);

#if 1
    HECMW_varray_int_rmdup(shared_nodes + i);
#endif
    n_shared_nodes[i] = HECMW_varray_int_nval(shared_nodes + i);

    /* number of shared nodes after removing duplication */
    send_buf[4*i+3] = n_shared_nodes[i];

    /* send sizes */
    irank = mesh->neighbor_pe[i];
    tag = 9901;
    HECMW_Isend(send_buf + 4*i, 4, HECMW_INT, irank, tag, comm, requests + i);
  }

  /* recieve size and check */
  for (i = 0; i < mesh->n_neighbor_pe; i++) {
    irank = mesh->neighbor_pe[i];
    tag = 9901;
    HECMW_Recv(recv_buf, 4, HECMW_INT, irank, tag, comm, statuses + i);

    if (recv_buf[0] != send_buf[4*i]) {
      HECMW_log(HECMW_LOG_ERROR,
                "number of shared elements %d not match with %d in rank %d\n",
                send_buf[4*i], recv_buf[0], irank);
      n_error++;
    }
    if (recv_buf[1] != send_buf[4*i+1]) {
      HECMW_log(HECMW_LOG_ERROR,
                "number of refined shared elements %d not match with %d in rank %d\n",
                send_buf[4*i+1], recv_buf[1], irank);
    }
    if (recv_buf[2] != send_buf[4*i+2]) {
      HECMW_log(HECMW_LOG_ERROR,
                "number of new shared nodes before rmdup %d not match with %d in rank %d\n",
                send_buf[4*i+2], recv_buf[2], irank);
      n_error++;
    }
    if (recv_buf[3] != send_buf[4*i+3]) {
      HECMW_log(HECMW_LOG_ERROR,
                "number of new shared nodes %d not match with %d in rank %d\n",
                send_buf[4*i+3], recv_buf[3], irank);
      n_error++;
      write_refined_shared_in_ucd(mesh, ref_mesh, shared, shared_nodes, i);
    }
  }
  HECMW_Waitall(mesh->n_neighbor_pe, requests, statuses);

  HECMW_Allreduce(&n_error, &n_error_g, 1, HECMW_INT, HECMW_SUM, comm);
  if (n_error_g > 0) return HECMW_ERROR;

  /*
   * check rank of new shared nodes
   */

  n_error = 0;

  /* send rank of new shared nodes */
  for (i = 0; i < mesh->n_neighbor_pe; i++) {
    HECMW_varray_int_init(shared_ranks + i);
    /* HECMW_varray_int_resize(shared_ranks + i, n_shared_nodes[i]); */
    for (j = 0; j < n_shared_nodes[i]; j++) {
      nid = HECMW_varray_int_get(shared_nodes + i, j);
      rank = ref_mesh->node_ID[2 * (nid - 1) + 1];
      /* HECMW_varray_int_set(shared_ranks + i, j, rank); */
      HECMW_varray_int_append(shared_ranks + i, rank);
    }

    irank = mesh->neighbor_pe[i];
    tag = 9902;
    int *send_ranks = HECMW_varray_int_get_v(shared_ranks + i);
    HECMW_Isend(send_ranks, n_shared_nodes[i], HECMW_INT, irank, tag, comm, requests + i);
  }

  /* recieve rank of new shared nodes */
  for (i = 0; i < mesh->n_neighbor_pe; i++) {
    int *recv_ranks = (int *)HECMW_malloc(sizeof(int) * n_shared_nodes[i]);
    if (recv_ranks == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }

    irank = mesh->neighbor_pe[i];
    tag = 9902;
    HECMW_Recv(recv_ranks, n_shared_nodes[i], HECMW_INT, irank, tag, comm, statuses + i);

    for (j = 0; j < n_shared_nodes[i]; j++) {
      if (recv_ranks[j] != HECMW_varray_int_get(shared_ranks + i, j)) {
        HECMW_log(HECMW_LOG_ERROR,
                  "rank of new shared node %d: %d not match with %d in rank %d\n",
                  HECMW_varray_int_get(shared_nodes + i, j),
                  HECMW_varray_int_get(shared_ranks + i, j), recv_ranks[j], irank);
        n_error++;
      }
    }
    HECMW_free(recv_ranks);
  }
  HECMW_Waitall(mesh->n_neighbor_pe, requests, statuses);

  HECMW_Allreduce(&n_error, &n_error_g, 1, HECMW_INT, HECMW_SUM, comm);
  if (n_error_g > 0) return HECMW_ERROR;

  /* free memory */
  for (i = 0; i < mesh->n_neighbor_pe; i++) {
    HECMW_varray_int_finalize(shared_nodes + i);
    HECMW_varray_int_finalize(shared_ranks + i);
  }
  HECMW_free(shared_nodes);
  HECMW_free(shared_ranks);
  HECMW_free(n_shared_nodes);
  HECMW_free(requests);
  HECMW_free(statuses);

  HECMW_log(HECMW_LOG_DEBUG, "Finished check_node_rank.\n");

  return HECMW_SUCCESS;
}

static int count_internal(int my_rank, int n, int *ID) {
  int i, count = 0;
  for (i = 0; i < n; i++) {
    if (ID[2 * i + 1] == my_rank) {
      ID[2 * i] = ++count;
    }
  }
  return count;
}

static int count_purely_external(int n, int *ID) {
  int i, count = 0;
  for (i = 0; i < n; i++) {
    if (ID[2 * i + 1] < 0) {
      count++;
    }
  }
  return count;
}

static int *new_internal_list(int my_rank, int n, int n_internal, int *ID) {
  int i, k;
  int *internal_list;

  internal_list = (int *)HECMW_malloc(sizeof(int) * n_internal);
  if (internal_list == NULL) {
    HECMW_set_error(errno, "");
    return NULL;
  }

  for (i = 0, k = 0; i < n; i++) {
    if (ID[2 * i + 1] == my_rank) {
      internal_list[k] = i + 1;
      k++;
    }
  }
  return internal_list;
}

static int create_index_item(int n, const struct hecmw_varray_int *arrays,
                             int **index, int **item) {
  int i, j, k;

  /* index */
  *index = (int *)HECMW_malloc(sizeof(int) * (n + 1));
  if (*index == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  (*index)[0] = 0;
  for (i = 0; i < n; i++) {
    (*index)[i + 1] = (*index)[i] + HECMW_varray_int_nval(arrays + i);
    if ((*index)[i + 1] < (*index)[i]) {
      HECMW_log(HECMW_LOG_ERROR,
                "Integer overflow detected while creating index array\n");
      return HECMW_ERROR;
    }
  }

  /* item */
  if ((*index)[n] == 0) {
    (*item) = NULL;
    return HECMW_SUCCESS;
  }
  *item = (int *)HECMW_malloc(sizeof(int) * (*index)[n]);
  if (*item == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  k = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < HECMW_varray_int_nval(arrays + i); j++) {
      (*item)[k] = HECMW_varray_int_get(arrays + i, j);
      k++;
    }
  }
  return HECMW_SUCCESS;
}

static int new_shared_export_import(const struct hecmwST_local_mesh *mesh,
                                    const struct hecmw_varray_int *shared,
                                    struct hecmwST_local_mesh *ref_mesh) {
  struct hecmw_varray_int *new_shared, *new_export, *new_import;
  int i;

  /* prepare new shared, import and export lists */
  new_shared = (struct hecmw_varray_int *)HECMW_malloc(
      sizeof(struct hecmw_varray_int) * mesh->n_neighbor_pe);
  if (new_shared == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  new_export = (struct hecmw_varray_int *)HECMW_malloc(
      sizeof(struct hecmw_varray_int) * mesh->n_neighbor_pe);
  if (new_export == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  new_import = (struct hecmw_varray_int *)HECMW_malloc(
      sizeof(struct hecmw_varray_int) * mesh->n_neighbor_pe);
  if (new_import == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }

  /* for each refined shared list */
  for (i = 0; i < mesh->n_neighbor_pe; i++) {
    int j;
    int rank_i = mesh->neighbor_pe[i];

    HECMW_varray_int_init(new_shared + i);
    HECMW_varray_int_init(new_export + i);
    HECMW_varray_int_init(new_import + i);

    /* for each element in the list */
    for (j = 0; j < HECMW_varray_int_nval(shared + i); j++) {
      struct hecmw_varray_int exp_cand, imp_cand;

      int eid = HECMW_varray_int_get(shared + i, j);
      int nn =
          ref_mesh->elem_node_index[eid] - ref_mesh->elem_node_index[eid - 1];
      int *elem_nodes =
          ref_mesh->elem_node_item + ref_mesh->elem_node_index[eid - 1];
      int etype = ref_mesh->elem_type[eid - 1];

      int min_rank        = mesh->n_subdomain;
      int myrank_included = 0;
      int rank_i_included = 0;
      int k;

      HECMW_varray_int_init(&exp_cand);
      HECMW_varray_int_init(&imp_cand);

      /* for each node in the element */
      for (k = 0; k < nn; k++) {
        int rank = ref_mesh->node_ID[2 * (elem_nodes[k] - 1) + 1];

        /* for purely external nodes, rank (0, 1, ...) is stored as -rank-1 (-1,
         * -2, ...) */
        if (rank < 0) rank = -rank - 1;

        if (rank < min_rank) min_rank = rank;

        if (rank == mesh->my_rank) {
          myrank_included = 1;
          HECMW_varray_int_append(&exp_cand, elem_nodes[k]);
        } else if (rank == rank_i) {
          rank_i_included = 1;
          HECMW_varray_int_append(&imp_cand, elem_nodes[k]);
        }
      }
      HECMW_assert(0 <= min_rank && min_rank < mesh->n_subdomain);

      ref_mesh->elem_ID[2 * (eid - 1) + 1] = min_rank;

      if (HECMW_is_etype_patch(etype)) {
        if (myrank_included || rank_i_included) HECMW_varray_int_append(new_shared + i, eid);
        if (myrank_included) HECMW_varray_int_cat(new_export + i, &exp_cand);
        if (rank_i_included) HECMW_varray_int_cat(new_import + i, &imp_cand);
      } else {
        if (myrank_included && rank_i_included) { /* the element is shared */
          HECMW_varray_int_append(new_shared + i, eid);
          HECMW_varray_int_cat(new_export + i, &exp_cand);
          HECMW_varray_int_cat(new_import + i, &imp_cand);
        } else if (myrank_included == 0) {
          /* for purely external elements, rank (0, 1, ...) is stored as -rank-1
           * (-1, -2, ...) */
          ref_mesh->elem_ID[2 * (eid - 1) + 1] = -min_rank - 1;
        }
      }

      HECMW_varray_int_finalize(&exp_cand);
      HECMW_varray_int_finalize(&imp_cand);
    }
    /* remove duplication */
    HECMW_varray_int_rmdup(new_export + i);
    HECMW_varray_int_rmdup(new_import + i);
  }

  /* new shared */
  if (create_index_item(mesh->n_neighbor_pe, new_shared, &(ref_mesh->shared_index),
                        &(ref_mesh->shared_item)) != HECMW_SUCCESS) {
    HECMW_log(HECMW_LOG_ERROR, "Create shared_index and shared_item failed\n");
    return HECMW_ERROR;
  }
  HECMW_log(HECMW_LOG_DEBUG, "Total number of shared elements = %d\n",
            ref_mesh->shared_index[mesh->n_neighbor_pe]);

  /* new export */
  if (create_index_item(mesh->n_neighbor_pe, new_export, &(ref_mesh->export_index),
                        &(ref_mesh->export_item)) != HECMW_SUCCESS) {
    HECMW_log(HECMW_LOG_ERROR, "Create export_index and export_item failed\n");
    return HECMW_ERROR;
  }
  HECMW_log(HECMW_LOG_DEBUG, "Total number of export nodes = %d\n",
            ref_mesh->export_index[mesh->n_neighbor_pe]);

  /* new import */
  if (create_index_item(mesh->n_neighbor_pe, new_import, &(ref_mesh->import_index),
                        &(ref_mesh->import_item)) != HECMW_SUCCESS) {
    HECMW_log(HECMW_LOG_ERROR, "Create import_index and import_item failed\n");
    return HECMW_ERROR;
  }
  HECMW_log(HECMW_LOG_DEBUG, "Total number of import nodes = %d\n",
            ref_mesh->import_index[mesh->n_neighbor_pe]);

  /* deallocate new shared, import and export lists */
  for (i = 0; i < mesh->n_neighbor_pe; i++) {
    HECMW_varray_int_finalize(new_shared + i);
    HECMW_varray_int_finalize(new_export + i);
    HECMW_varray_int_finalize(new_import + i);
  }
  HECMW_free(new_shared);
  HECMW_free(new_export);
  HECMW_free(new_import);

  return HECMW_SUCCESS;
}

static int rebuild_comm_tables_serial(const struct hecmwST_local_mesh *mesh,
                                      struct hecmwST_local_mesh *ref_mesh) {
  int i;

  HECMW_assert(mesh->hecmw_flag_parttype == HECMW_FLAG_PARTTYPE_UNKNOWN ||
               mesh->hecmw_flag_parttype == HECMW_FLAG_PARTTYPE_NODEBASED);

  /* nn_internal, n_node */
  ref_mesh->nn_internal = count_internal(mesh->my_rank, ref_mesh->n_node_gross, ref_mesh->node_ID);
  ref_mesh->n_node = ref_mesh->n_node_gross;
  ref_mesh->nn_middle = ref_mesh->n_node;

  HECMW_log(HECMW_LOG_DEBUG, "nn_internal = %d, n_node = %d, nn_middle = %d\n",
            ref_mesh->nn_internal, ref_mesh->n_node, ref_mesh->nn_middle);

  /* node_internal_list */
  ref_mesh->node_internal_list =
      (int *)HECMW_malloc(sizeof(int) * ref_mesh->nn_internal);
  if (ref_mesh->node_internal_list == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < ref_mesh->nn_internal; i++) {
    ref_mesh->node_internal_list[i] = i + 1;
  }

  /* ne_internal, n_elem */
  ref_mesh->ne_internal = count_internal(mesh->my_rank, ref_mesh->n_elem_gross, ref_mesh->elem_ID);
  ref_mesh->n_elem = ref_mesh->n_elem_gross;

  HECMW_log(HECMW_LOG_DEBUG, "ne_internal = %d, n_elem = %d\n",
            ref_mesh->ne_internal, ref_mesh->n_elem);

  /* elem_internal_list */
  ref_mesh->elem_internal_list =
      (int *)HECMW_malloc(sizeof(int) * ref_mesh->ne_internal);
  if (ref_mesh->elem_internal_list == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < ref_mesh->ne_internal; i++) {
    ref_mesh->elem_internal_list[i] = i + 1;
  }
  return HECMW_SUCCESS;
}

static int mark_purely_external_nodes(int my_rank,
                                      struct hecmwST_local_mesh *ref_mesh) {
  struct hecmw_bit_array mark;
  int i;

  if (HECMW_bit_array_init(&mark, ref_mesh->n_node_gross) != HECMW_SUCCESS) {
    return HECMW_ERROR;
  }
  /* mark all nodes in the new import lists */
  for (i = 0; i < ref_mesh->import_index[ref_mesh->n_neighbor_pe]; i++) {
    HECMW_bit_array_set(&mark, ref_mesh->import_item[i] - 1);
  }
  /* scan for all external, unmarked nodes and set their rank to -rank-1 */
  for (i = 0; i < ref_mesh->n_node_gross; i++) {
    int rank = ref_mesh->node_ID[2 * i + 1];
    if (rank == my_rank || rank < 0) continue;
    if (!HECMW_bit_array_get(&mark, i)) {
      /* for purely external nodes, rank (0, 1, ...) is stored as -rank-1 (-1,
       * -2, ...) */
      ref_mesh->node_ID[2 * i + 1] = -rank - 1;
    }
  }
  HECMW_bit_array_finalize(&mark);
  return HECMW_SUCCESS;
}

static int rebuild_comm_tables(const struct hecmwST_local_mesh *mesh,
                               struct hecmwST_local_mesh *ref_mesh) {
  int i;
  struct hecmw_varray_int *shared;

  if (mesh->n_neighbor_pe == 0) { /* SERIAL */
    return rebuild_comm_tables_serial(mesh, ref_mesh);
  }

  /* PARALLEL */
  HECMW_log(HECMW_LOG_DEBUG, "Started rebuilding communication tables...\n");

  /* check part-type */
  if (mesh->hecmw_flag_parttype != HECMW_FLAG_PARTTYPE_NODEBASED) {
    HECMW_log(HECMW_LOG_ERROR,
              "Partitioning-type must be node-based for refinement.\n");
    return HECMW_ERROR;
  }

  /* Get refined shared element list */
  shared = (struct hecmw_varray_int *)HECMW_malloc(
      sizeof(struct hecmw_varray_int) * mesh->n_neighbor_pe);
  if (shared == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < mesh->n_neighbor_pe; i++) {
    HECMW_varray_int_init(shared + i);
  }
  if (get_refined_shared(mesh, ref_mesh, shared) != HECMW_SUCCESS) {
    return HECMW_ERROR;
  }

  /* Determine rank of new nodes */
  if (determine_node_rank(mesh, shared, ref_mesh) != HECMW_SUCCESS) {
    return HECMW_ERROR;
  }

  /* Check rank of new nodes with neighbor processes */
  if (check_node_rank(mesh, shared, ref_mesh) != HECMW_SUCCESS) {
    return HECMW_ERROR;
  }

  /* nn_internal */
  ref_mesh->nn_internal =
      count_internal(mesh->my_rank, ref_mesh->n_node_gross, ref_mesh->node_ID);

  HECMW_assert(mesh->node_internal_list ==
               NULL); /* because PARTTYPE is NODEBASED. */

  /* new shared, import and export lists */
  if (new_shared_export_import(mesh, shared, ref_mesh) != HECMW_SUCCESS) {
    return HECMW_ERROR;
  }

  /* deallocate refined shared list */
  for (i = 0; i < mesh->n_neighbor_pe; i++) {
    HECMW_varray_int_finalize(shared + i);
  }
  HECMW_free(shared);

  /* ne_internal */
  ref_mesh->ne_internal =
      count_internal(mesh->my_rank, ref_mesh->n_elem_gross, ref_mesh->elem_ID);

  /* n_elem */
  ref_mesh->n_elem =
      ref_mesh->n_elem_gross -
      count_purely_external(ref_mesh->n_elem_gross, ref_mesh->elem_ID);

  /* elem_internal_list */
  if (ref_mesh->ne_internal > 0) {
    ref_mesh->elem_internal_list =
        new_internal_list(mesh->my_rank, ref_mesh->n_elem_gross,
                          ref_mesh->ne_internal, ref_mesh->elem_ID);
    if (ref_mesh->elem_internal_list == NULL) return HECMW_ERROR;
  }

  /* mark purely external nodes */
  if (mark_purely_external_nodes(mesh->my_rank, ref_mesh) != HECMW_SUCCESS) {
    return HECMW_ERROR;
  }

  /* n_node */
  ref_mesh->n_node =
      ref_mesh->n_node_gross -
      count_purely_external(ref_mesh->n_node_gross, ref_mesh->node_ID);
  ref_mesh->nn_middle = ref_mesh->n_node;  /* not set properly (since it's never used) */

  HECMW_log(HECMW_LOG_DEBUG, "nn_internal = %d, n_node = %d nn_middle = %d, \n",
            ref_mesh->nn_internal, ref_mesh->n_node, ref_mesh->nn_middle);
  HECMW_log(HECMW_LOG_DEBUG, "ne_internal = %d, nelem = %d\n",
            ref_mesh->ne_internal, ref_mesh->n_elem);

  HECMW_log(HECMW_LOG_DEBUG, "Finished rebuilding communication tables.\n");
  return HECMW_SUCCESS;
}

static int check_comm_table_len(const struct hecmwST_local_mesh *ref_mesh) {
  int len_tot;
  int i, j, irank, js, je, len, tag, nsend;
  HECMW_Request *requests;
  HECMW_Status *statuses;
  HECMW_Comm comm;
  int *send_buf;
  int n_error = 0, n_error_g;

  if (ref_mesh->n_neighbor_pe == 0) return HECMW_SUCCESS;

  HECMW_log(HECMW_LOG_DEBUG, "Started checking communication tables...\n");

  send_buf = (int *)HECMW_malloc(sizeof(int) * ref_mesh->n_neighbor_pe);
  if (send_buf == 0) {
    HECMW_set_error(errno, "");
    HECMW_abort(comm);
  }

  comm = HECMW_comm_get_comm();

  requests = (HECMW_Request *)HECMW_malloc(sizeof(HECMW_Request) *
                                           ref_mesh->n_neighbor_pe);
  if (requests == NULL) {
    HECMW_set_error(errno, "");
    HECMW_abort(comm);
  }
  statuses = (HECMW_Status *)HECMW_malloc(sizeof(HECMW_Status) *
                                          ref_mesh->n_neighbor_pe);
  if (statuses == NULL) {
    HECMW_set_error(errno, "");
    HECMW_abort(comm);
  }

  /* export and import */

  for (i = 0; i < ref_mesh->n_neighbor_pe; i++) {
    irank     = ref_mesh->neighbor_pe[i];
    js        = ref_mesh->export_index[i];
    je        = ref_mesh->export_index[i + 1];
    send_buf[i] = je - js;
    /* isend number of export nodes */
    tag = 0;
    HECMW_Isend(send_buf + i, 1, HECMW_INT, irank, tag, comm, requests + i);
  }
  for (i = 0; i < ref_mesh->n_neighbor_pe; i++) {
    irank     = ref_mesh->neighbor_pe[i];
    js        = ref_mesh->import_index[i];
    je        = ref_mesh->import_index[i + 1];
    /* recv number of import nodes */
    tag = 0;
    HECMW_Recv(&len, 1, HECMW_INT, irank, tag, comm, statuses + i);
    if (len != je - js) {
      HECMW_log(HECMW_LOG_ERROR,
                "inconsistent length of import (%d) with export on rank %d (%d)\n",
                je-js, irank, len);
      n_error++;
    }
  }
  HECMW_Waitall(ref_mesh->n_neighbor_pe, requests, statuses);

  /* shared */

  nsend = 0;
  for (i = 0; i < ref_mesh->n_neighbor_pe; i++) {
    irank   = ref_mesh->neighbor_pe[i];
    if (irank < ref_mesh->my_rank) continue;
    js      = ref_mesh->shared_index[i];
    je      = ref_mesh->shared_index[i + 1];
    send_buf[i] = je - js;
    /* isend number of shared elems */
    tag = 1;
    HECMW_Isend(send_buf + i, 1, HECMW_INT, irank, tag, comm, requests + nsend);
    nsend++;
  }
  for (i = 0; i < ref_mesh->n_neighbor_pe; i++) {
    irank   = ref_mesh->neighbor_pe[i];
    if (ref_mesh->my_rank < irank) continue;
    js      = ref_mesh->shared_index[i];
    je      = ref_mesh->shared_index[i + 1];
    /* recv number of shared elems */
    tag = 1;
    HECMW_Recv(&len, 1, HECMW_INT, irank, tag, comm, statuses);
    if (len != je - js) {
      HECMW_log(HECMW_LOG_ERROR,
                "inconsistent length of shared (%d) with rank %d (%d)\n",
                je-js, irank, len);
      n_error++;
    }
  }
  HECMW_Waitall(nsend, requests, statuses);

  HECMW_free(send_buf);
  HECMW_free(requests);
  HECMW_free(statuses);

  HECMW_Allreduce(&n_error, &n_error_g, 1, HECMW_INT, HECMW_SUM, comm);
  if (n_error_g > 0) return HECMW_ERROR;

  HECMW_log(HECMW_LOG_DEBUG, "Finished checking communication tables.\n");
  return HECMW_SUCCESS;
}

static int rebuild_ID_external(struct hecmwST_local_mesh *ref_mesh) {
  int len_tot;
  int *sendbuf, *recvbuf, *srbuf;
  int i, j, irank, js, je, len, nid, cnt, tag, nsend, rank, lid;
  int *item_p, *sendbuf_p, *recvbuf_p, *srbuf_p;
  HECMW_Request *requests;
  HECMW_Status *statuses;
  HECMW_Comm comm;

  if (ref_mesh->n_neighbor_pe == 0) return HECMW_SUCCESS;

  HECMW_log(HECMW_LOG_DEBUG, "Started rebuilding external IDs...\n");

  comm = HECMW_comm_get_comm();

  requests = (HECMW_Request *)HECMW_malloc(sizeof(HECMW_Request) *
                                           ref_mesh->n_neighbor_pe);
  if (requests == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  statuses = (HECMW_Status *)HECMW_malloc(sizeof(HECMW_Status) *
                                          ref_mesh->n_neighbor_pe);
  if (statuses == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }

  /* node_ID (for external nodes) */

  /* send local IDs of export nodes */
  len_tot = ref_mesh->export_index[ref_mesh->n_neighbor_pe];
  sendbuf = (int *)HECMW_malloc(sizeof(int) * len_tot);
  if (sendbuf == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < ref_mesh->n_neighbor_pe; i++) {
    irank     = ref_mesh->neighbor_pe[i];
    js        = ref_mesh->export_index[i];
    je        = ref_mesh->export_index[i + 1];
    len       = je - js;
    item_p    = ref_mesh->export_item + js;
    sendbuf_p = sendbuf + js;
    for (j = 0; j < len; j++) {
      nid = item_p[j] - 1;
      lid = ref_mesh->node_ID[2 * nid];
      HECMW_assert(0 < lid);
      HECMW_assert(lid <= ref_mesh->nn_internal);
      sendbuf_p[j] = lid;
    }
    /* isend local id of export nodes */
    tag = 0;
    HECMW_Isend(sendbuf_p, len, HECMW_INT, irank, tag, comm, requests + i);
  }

  /* recv local IDs of import nodes */
  len_tot = ref_mesh->import_index[ref_mesh->n_neighbor_pe];
  recvbuf = (int *)HECMW_malloc(sizeof(int) * len_tot);
  if (recvbuf == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < ref_mesh->n_neighbor_pe; i++) {
    irank     = ref_mesh->neighbor_pe[i];
    js        = ref_mesh->import_index[i];
    je        = ref_mesh->import_index[i + 1];
    len       = je - js;
    item_p    = ref_mesh->import_item + js;
    recvbuf_p = recvbuf + js;
    /* recv local id of import nodes */
    tag = 0;
    HECMW_Recv(recvbuf_p, len, HECMW_INT, irank, tag, comm, statuses + i);
    /* set node_ID[2*j] */
    for (j = 0; j < len; j++) {
      nid = item_p[j] - 1;
      lid = recvbuf_p[j];
      HECMW_assert(0 < lid);
      ref_mesh->node_ID[2 * nid] = lid;
    }
  }

  HECMW_Waitall(ref_mesh->n_neighbor_pe, requests, statuses);
  HECMW_free(sendbuf);
  HECMW_free(recvbuf);

  /* elem_ID (for external elements) */

  len_tot = ref_mesh->shared_index[ref_mesh->n_neighbor_pe];
  srbuf   = (int *)HECMW_malloc(sizeof(int) * len_tot);
  if (srbuf == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  nsend = 0;
  for (i = 0; i < ref_mesh->n_neighbor_pe; i++) {
    irank   = ref_mesh->neighbor_pe[i];
    if (irank < ref_mesh->my_rank) continue;
    js      = ref_mesh->shared_index[i];
    je      = ref_mesh->shared_index[i + 1];
    len     = je - js;
    srbuf_p = srbuf + js;
    item_p  = ref_mesh->shared_item + js;
    for (j = 0; j < len; j++) {
      nid  = item_p[j] - 1;
      rank = ref_mesh->elem_ID[2 * nid + 1];
      lid  = ref_mesh->elem_ID[2 * nid];
      if (rank == ref_mesh->my_rank && lid > 0) {
        HECMW_assert(lid <= ref_mesh->ne_internal);
        srbuf_p[j] = lid;
      } else {
        srbuf_p[j] = -1;
      }
    }
    /* isend list of local_ID of those elems */
    tag = 1;
    HECMW_Isend(srbuf_p, len, HECMW_INT, irank, tag, comm, requests + nsend);
    nsend++;
  }
  for (i = 0; i < ref_mesh->n_neighbor_pe; i++) {
    irank   = ref_mesh->neighbor_pe[i];
    if (ref_mesh->my_rank < irank) continue;
    js      = ref_mesh->shared_index[i];
    je      = ref_mesh->shared_index[i + 1];
    len     = je - js;
    srbuf_p = srbuf + js;
    item_p  = ref_mesh->shared_item + js;
    /* recv list of local_ID of those elems */
    tag = 1;
    HECMW_Recv(srbuf_p, len, HECMW_INT, irank, tag, comm, statuses);
    /* set elem_ID[2*j] */
    for (j = 0; j < len; j++) {
      lid = srbuf_p[j];
      if (lid < 0) continue;
      nid  = item_p[j] - 1;
      rank = ref_mesh->elem_ID[2 * nid + 1];
      HECMW_assert(rank == irank || rank == -irank - 1);
      if (rank < 0) continue;
      ref_mesh->elem_ID[2 * nid] = lid;
    }
  }
  HECMW_Waitall(nsend, requests, statuses);
  HECMW_free(srbuf);

  HECMW_free(requests);
  HECMW_free(statuses);

  HECMW_log(HECMW_LOG_DEBUG, "Finished rebuilding external IDs.\n");
  return HECMW_SUCCESS;
}

static int rebuild_refine_origin(const struct hecmwST_local_mesh *mesh,
                                 struct hecmwST_local_mesh *ref_mesh) {
  const struct hecmwST_refine_origin *reforg;
  struct hecmwST_refine_origin *ref_reforg;
  struct hecmw_varray_int ro_item;
  const int *ro_item_v;
  int n_refine, offset, cnt, cnt2, i, j;

  reforg   = mesh->refine_origin;
  n_refine = mesh->n_refine + 1;

  ref_reforg = (struct hecmwST_refine_origin *)HECMW_malloc(
      sizeof(struct hecmwST_refine_origin));
  if (ref_reforg == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }

  ref_reforg->index = (int *)HECMW_malloc(sizeof(int) * (n_refine + 1));
  if (ref_reforg->index == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  ref_reforg->index[0] = 0;
  for (i = 1; i < n_refine; i++) {
    ref_reforg->index[i] = reforg->index[i];
  }
  ref_reforg->index[n_refine] =
      ref_reforg->index[n_refine - 1] + ref_mesh->n_node;

  ref_reforg->item_index =
      (int *)HECMW_malloc(sizeof(int) * (ref_reforg->index[n_refine] + 1));
  if (ref_reforg->item_index == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }

  /* copy original-node infomation up to previous refinement */
  ref_reforg->item_index[0] = 0;
  for (i = 1; i <= ref_reforg->index[n_refine - 1]; i++) {
    ref_reforg->item_index[i] = reforg->item_index[i];
  }
  offset = ref_reforg->index[n_refine - 1];
  cnt    = ref_reforg->item_index[offset];
  /* fprintf(stderr, "offset=%d, cnt=%d\n", offset, cnt); */

  /* get original nodes of newly generated nodes at current refinement */
  HECMW_varray_int_init(&ro_item);
  for (i = 0; i < ref_mesh->n_node; i++) {
    int iold;
    iold = ref_mesh->node_new2old ? ref_mesh->node_new2old[i] : i;
    if (iold < mesh->n_node_gross) {
      int org;
      org = mesh->node_old2new ? mesh->node_old2new[iold] + 1 : iold + 1;
      HECMW_varray_int_append(&ro_item, org);
      cnt++;
    } else {
      int original[HECMW_MAX_NODE_MAX];
      int orig_type_rcap = rcapGetOriginal(iold + 1, original);
      int nn;
      HECMW_assert(orig_type_rcap != RCAP_UNKNOWNTYPE);
      nn = get_rcap_elem_max_node(orig_type_rcap);
      for (j = 0; j < nn; j++) {
        int org                     = original[j];
        if (mesh->node_old2new) org = mesh->node_old2new[org - 1] + 1;
        HECMW_varray_int_append(&ro_item, org);
      }
      cnt += nn;
    }
    ref_reforg->item_index[offset + i + 1] = cnt;
  }

  ref_reforg->item_item = (int *)HECMW_malloc(sizeof(int) * cnt);
  if (ref_reforg->item_item == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }

  /* copy original nodes generated until previous refinements */
  cnt = ref_reforg->item_index[ref_reforg->index[n_refine - 1]];
  for (i = 0; i < cnt; i++) {
    ref_reforg->item_item[i] = reforg->item_item[i];
  }
  /* copy original nodes generated at current refinement */
  cnt2      = HECMW_varray_int_nval(&ro_item);
  ro_item_v = HECMW_varray_int_get_cv(&ro_item);
  for (i = 0; i < cnt2; i++) {
    ref_reforg->item_item[cnt + i] = ro_item_v[i];
  }
  HECMW_varray_int_finalize(&ro_item);
#if 0
	fprintf(stderr, "refine_origin->index:\n");
	for( i=0; i <= n_refine; i++ ) {
		fprintf(stderr, " %d", ref_reforg->index[i]);
		if( i % 10 == 9 ) fprintf(stderr, "\n");
	}
	if( i % 10 != 0 ) fprintf(stderr, "\n");
	fprintf(stderr, "refine_origin->item_index:\n");
	for( i=0; i <= ref_reforg->index[n_refine]; i++) {
		fprintf(stderr, " %d", ref_reforg->item_index[i]);
		if( i % 10 == 9 ) fprintf(stderr, "\n");
	}
	if( i % 10 != 0 ) fprintf(stderr, "\n");
	fprintf(stderr, "refine_origin->item_item:\n");
	for( i=0; i < ref_reforg->item_index[ref_reforg->index[n_refine]]; i++ ) {
		fprintf(stderr, " %d", ref_reforg->item_item[i]);
		if( i % 10 == 9 ) fprintf(stderr, "\n");
	}
	if( i % 10 != 0 ) fprintf(stderr, "\n");
#endif
  ref_mesh->refine_origin = ref_reforg;
  return HECMW_SUCCESS;
}

static int rebuild_n_node_refine_hist(const struct hecmwST_local_mesh *mesh,
                                      struct hecmwST_local_mesh *ref_mesh) {
  int i;

  ref_mesh->n_node_refine_hist =
    (int *)HECMW_malloc((mesh->n_refine + 1) * sizeof(int));
  if (ref_mesh->n_node_refine_hist == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }

  for (i = 0; i < mesh->n_refine; i++) {
    ref_mesh->n_node_refine_hist[i] = mesh->n_node_refine_hist[i];
  }
  ref_mesh->n_node_refine_hist[mesh->n_refine] = ref_mesh->n_node_gross;

  return HECMW_SUCCESS;
}

static int rebuild_refine_info(const struct hecmwST_local_mesh *mesh,
                               struct hecmwST_local_mesh *ref_mesh) {
  ref_mesh->n_refine = mesh->n_refine + 1;

  if (rebuild_refine_origin(mesh, ref_mesh) != HECMW_SUCCESS)
    return HECMW_ERROR;
  if (rebuild_n_node_refine_hist(mesh, ref_mesh) != HECMW_SUCCESS)
    return HECMW_ERROR;

  return HECMW_SUCCESS;
}

static int renumber_nodes_generate_tables(struct hecmwST_local_mesh *mesh);
static int renumber_elements_generate_tables(struct hecmwST_local_mesh *mesh);

static int rebuild_info(const struct hecmwST_local_mesh *mesh,
                        struct hecmwST_local_mesh *ref_mesh) {
  HECMW_log(HECMW_LOG_DEBUG, "Started rebuilding info...\n");

  if (rebuild_elem_info(mesh, ref_mesh) != HECMW_SUCCESS) return HECMW_ERROR;
  if (rebuild_node_info(mesh, ref_mesh) != HECMW_SUCCESS) return HECMW_ERROR;
  if (rebuild_comm_tables(mesh, ref_mesh) != HECMW_SUCCESS) return HECMW_ERROR;
  if (check_comm_table_len(ref_mesh) != HECMW_SUCCESS) {
    HECMW_log(HECMW_LOG_ERROR, "Check communication table failed. Checking original mesh...\n");
    if (check_comm_table_len(mesh) != HECMW_SUCCESS) {
      HECMW_log(HECMW_LOG_ERROR, "Original mesh has error\n");
    } else {
      HECMW_log(HECMW_LOG_ERROR, "Original mesh is OK\n");
    }
    return HECMW_ERROR;
  }
  if (rebuild_ID_external(ref_mesh) != HECMW_SUCCESS) return HECMW_ERROR;

  if (renumber_nodes_generate_tables(ref_mesh) != HECMW_SUCCESS)
    return HECMW_ERROR;
  if (renumber_elements_generate_tables(ref_mesh) != HECMW_SUCCESS)
    return HECMW_ERROR;

  if (rebuild_refine_info(mesh, ref_mesh) != HECMW_SUCCESS) return HECMW_ERROR;

  HECMW_log(HECMW_LOG_DEBUG, "Finished rebuilding info.\n");
  return HECMW_SUCCESS;
}

/*============================================================================*/
/*                                                                            */
/*  copy unchanging infomation                                                */
/*                                                                            */
/*============================================================================*/
static int copy_global_info(const struct hecmwST_local_mesh *mesh,
                            struct hecmwST_local_mesh *ref_mesh) {
  int i;

  ref_mesh->hecmw_flag_adapt       = mesh->hecmw_flag_adapt;
  ref_mesh->hecmw_flag_initcon     = mesh->hecmw_flag_initcon;
  ref_mesh->hecmw_flag_parttype    = mesh->hecmw_flag_parttype;
  ref_mesh->hecmw_flag_partdepth   = mesh->hecmw_flag_partdepth;
  ref_mesh->hecmw_flag_version     = mesh->hecmw_flag_version;
  ref_mesh->hecmw_flag_partcontact = mesh->hecmw_flag_partcontact;
  strcpy(ref_mesh->gridfile, mesh->gridfile);
  ref_mesh->hecmw_n_file       = mesh->hecmw_n_file;

  /* files */
  if (mesh->hecmw_n_file > 0) {
    ref_mesh->files = (char **)HECMW_calloc(mesh->hecmw_n_file, sizeof(char *));
    if (ref_mesh->files == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    for (i = 0; i < mesh->hecmw_n_file; i++) {
      ref_mesh->files[i] = HECMW_strdup(mesh->files[i]);
      if (ref_mesh->files[i] == NULL) {
        HECMW_set_error(errno, "");
        return HECMW_ERROR;
      }
    }
  } else {
    ref_mesh->files = NULL;
  }

  strcpy(ref_mesh->header, mesh->header);
  ref_mesh->zero_temp = mesh->zero_temp;

  return HECMW_SUCCESS;
}

static int copy_node_ndof(const struct hecmwST_local_mesh *mesh,
                          struct hecmwST_local_mesh *ref_mesh) {
  int i;

  ref_mesh->n_dof     = mesh->n_dof;
  ref_mesh->n_dof_grp = mesh->n_dof_grp;
  ref_mesh->n_dof_tot = mesh->n_dof_tot;

  if (mesh->n_dof_grp > 0) {
    /* node_dof_index */
    ref_mesh->node_dof_index =
        (int *)HECMW_malloc(sizeof(int) * (mesh->n_dof_grp + 1));
    if (ref_mesh->node_dof_index == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    for (i = 0; i < mesh->n_dof_grp + 1; i++) {
      ref_mesh->node_dof_index[i] = mesh->node_dof_index[i];
    }

    /* node_dof_item */
    ref_mesh->node_dof_item =
        (int *)HECMW_malloc(sizeof(int) * mesh->n_dof_grp);
    if (ref_mesh->node_dof_item == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    for (i = 0; i < mesh->n_dof_grp; i++) {
      ref_mesh->node_dof_item[i] = mesh->node_dof_item[i];
    }
  }
  return HECMW_SUCCESS;
}

static int copy_elem_etype(const struct hecmwST_local_mesh *mesh,
                           struct hecmwST_local_mesh *ref_mesh) {
  int i;

  ref_mesh->n_elem_type = mesh->n_elem_type;

  if (mesh->n_elem_type > 0) {
    /* elem_type_index */
    ref_mesh->elem_type_index =
        (int *)HECMW_malloc(sizeof(int) * (mesh->n_elem_type + 1));
    if (ref_mesh->elem_type_index == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }

    /* elem_type_item */
    ref_mesh->elem_type_item =
        (int *)HECMW_malloc(sizeof(int) * mesh->n_elem_type);
    if (ref_mesh->elem_type_item == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }

    ref_mesh->elem_type_index[0] = 0;
    for (i = 0; i < mesh->n_elem_type; i++) {
      int etype, ierror, ndiv, nelem;
      etype = mesh->elem_type_item[i];
      ndiv  = get_elem_ndiv(etype, &ierror);
      if (ierror) return HECMW_ERROR;
      nelem = mesh->elem_type_index[i + 1] - mesh->elem_type_index[i];
      ref_mesh->elem_type_index[i + 1] =
          ref_mesh->elem_type_index[i] + nelem * ndiv;
      ref_mesh->elem_type_item[i] = etype;
    }
  }
  return HECMW_SUCCESS;
}

static int copy_comm_info(const struct hecmwST_local_mesh *mesh,
                          struct hecmwST_local_mesh *ref_mesh) {
  int i;

  ref_mesh->zero          = mesh->zero;
  ref_mesh->HECMW_COMM    = mesh->HECMW_COMM;
  ref_mesh->PETOT         = mesh->PETOT;
  ref_mesh->PEsmpTOT      = mesh->PEsmpTOT;
  ref_mesh->my_rank       = mesh->my_rank;
  ref_mesh->errnof        = mesh->errnof;
  ref_mesh->n_subdomain   = mesh->n_subdomain;
  ref_mesh->n_neighbor_pe = mesh->n_neighbor_pe;

  if (mesh->n_neighbor_pe == 0) {
    ref_mesh->neighbor_pe = NULL;
    ref_mesh->import_item = NULL;
    ref_mesh->export_item = NULL;
    ref_mesh->shared_item = NULL;

    ref_mesh->import_index = (int *)HECMW_malloc(sizeof(int));
    if (ref_mesh->import_index == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    ref_mesh->import_index[0] = 0;

    ref_mesh->export_index = (int *)HECMW_malloc(sizeof(int));
    if (ref_mesh->export_index == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    ref_mesh->export_index[0] = 0;

    ref_mesh->shared_index = (int *)HECMW_malloc(sizeof(int));
    if (ref_mesh->shared_index == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    ref_mesh->shared_index[0] = 0;

    return HECMW_SUCCESS;
  }

  /* neighbor_pe */
  ref_mesh->neighbor_pe =
      (int *)HECMW_malloc(sizeof(int) * mesh->n_neighbor_pe);
  if (ref_mesh->neighbor_pe == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < mesh->n_neighbor_pe; i++) {
    ref_mesh->neighbor_pe[i] = mesh->neighbor_pe[i];
  }

  return HECMW_SUCCESS;
}

static int copy_adapt_info(const struct hecmwST_local_mesh *mesh,
                           struct hecmwST_local_mesh *ref_mesh) {
  if (mesh->hecmw_flag_adapt != 0) {
    HECMW_log(HECMW_LOG_ERROR, "Refinement of adaptation not supported\n");
    return HECMW_ERROR;
  }
  ref_mesh->coarse_grid_level       = 0;
  ref_mesh->n_adapt                 = 0;
  ref_mesh->when_i_was_refined_node = NULL;
  ref_mesh->when_i_was_refined_elem = NULL;
  ref_mesh->adapt_parent_type       = NULL;
  ref_mesh->adapt_type              = NULL;
  ref_mesh->adapt_level             = NULL;
  ref_mesh->adapt_parent            = NULL;
  ref_mesh->adapt_children_index    = NULL;
  ref_mesh->adapt_children_item     = NULL;
  return HECMW_SUCCESS;
}

static int copy_section_info(const struct hecmwST_section *sect,
                             struct hecmwST_section *ref_sect) {
  int i;

  ref_sect->n_sect = sect->n_sect;

  if (sect->n_sect == 0) {
    ref_sect->sect_type         = NULL;
    ref_sect->sect_opt          = NULL;
    ref_sect->sect_mat_ID_index = NULL;
    ref_sect->sect_mat_ID_item  = NULL;
    ref_sect->sect_I_index      = NULL;
    ref_sect->sect_I_item       = NULL;
    ref_sect->sect_R_index      = NULL;
    ref_sect->sect_R_item       = NULL;
    return HECMW_SUCCESS;
  }

  /* sect_type */
  ref_sect->sect_type = (int *)HECMW_malloc(sizeof(int) * sect->n_sect);
  if (ref_sect->sect_type == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < sect->n_sect; i++) {
    ref_sect->sect_type[i] = sect->sect_type[i];
  }

  /* sect_opt */
  ref_sect->sect_opt = (int *)HECMW_malloc(sizeof(int) * sect->n_sect);
  if (ref_sect->sect_opt == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < sect->n_sect; i++) {
    ref_sect->sect_opt[i] = sect->sect_opt[i];
  }

  /* sect_mat_ID_index */
  ref_sect->sect_mat_ID_index =
      (int *)HECMW_malloc(sizeof(int) * (sect->n_sect + 1));
  if (ref_sect->sect_mat_ID_index == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < sect->n_sect + 1; i++) {
    ref_sect->sect_mat_ID_index[i] = sect->sect_mat_ID_index[i];
  }

  /* sect_mat_ID_item */
  if (sect->sect_mat_ID_index[sect->n_sect] > 0) {
    ref_sect->sect_mat_ID_item = (int *)HECMW_malloc(
        sizeof(int) * sect->sect_mat_ID_index[sect->n_sect]);
    if (ref_sect->sect_mat_ID_item == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    for (i = 0; i < sect->sect_mat_ID_index[sect->n_sect]; i++) {
      ref_sect->sect_mat_ID_item[i] = sect->sect_mat_ID_item[i];
    }
  }

  /* sect_I_index */
  ref_sect->sect_I_index =
      (int *)HECMW_malloc(sizeof(int) * (sect->n_sect + 1));
  if (ref_sect->sect_I_index == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < sect->n_sect + 1; i++) {
    ref_sect->sect_I_index[i] = sect->sect_I_index[i];
  }

  /* sect_I_item */
  if (sect->sect_I_index[sect->n_sect] > 0) {
    ref_sect->sect_I_item =
        (int *)HECMW_malloc(sizeof(int) * sect->sect_I_index[sect->n_sect]);
    if (ref_sect->sect_I_item == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    for (i = 0; i < sect->sect_I_index[sect->n_sect]; i++) {
      ref_sect->sect_I_item[i] = sect->sect_I_item[i];
    }
  }

  /* sect_R_index */
  ref_sect->sect_R_index =
      (int *)HECMW_malloc(sizeof(int) * (sect->n_sect + 1));
  if (ref_sect->sect_R_index == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < sect->n_sect + 1; i++) {
    ref_sect->sect_R_index[i] = sect->sect_R_index[i];
  }

  /* sect_R_item */
  if (sect->sect_R_index[sect->n_sect] > 0) {
    ref_sect->sect_R_item = (double *)HECMW_malloc(
        sizeof(double) * sect->sect_R_index[sect->n_sect]);
    if (ref_sect->sect_R_item == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    for (i = 0; i < sect->sect_R_index[sect->n_sect]; i++) {
      ref_sect->sect_R_item[i] = sect->sect_R_item[i];
    }
  }

  return HECMW_SUCCESS;
}

static int copy_material_info(const struct hecmwST_material *mat,
                              struct hecmwST_material *ref_mat) {
  int i;

  ref_mat->n_mat = mat->n_mat;

  if (mat->n_mat == 0) {
    ref_mat->n_mat_item        = 0;
    ref_mat->n_mat_subitem     = 0;
    ref_mat->n_mat_table       = 0;
    ref_mat->mat_name          = NULL;
    ref_mat->mat_item_index    = NULL;
    ref_mat->mat_subitem_index = NULL;
    ref_mat->mat_table_index   = NULL;
    ref_mat->mat_val           = NULL;
    ref_mat->mat_temp          = NULL;
    return HECMW_SUCCESS;
  }

  ref_mat->n_mat_item    = mat->n_mat_item;
  ref_mat->n_mat_subitem = mat->n_mat_subitem;
  ref_mat->n_mat_table   = mat->n_mat_table;

  /* mat_name */
  ref_mat->mat_name = (char **)HECMW_malloc(sizeof(char *) * mat->n_mat);
  if (ref_mat->mat_name == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < mat->n_mat; i++) {
    ref_mat->mat_name[i] = HECMW_strdup(mat->mat_name[i]);
    if (ref_mat->mat_name[i] == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
  }

  /* mat_item_index */
  ref_mat->mat_item_index = (int *)HECMW_malloc(sizeof(int) * (mat->n_mat + 1));
  if (ref_mat->mat_item_index == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < mat->n_mat + 1; i++) {
    ref_mat->mat_item_index[i] = mat->mat_item_index[i];
  }

  /* mat_subitem_index */
  ref_mat->mat_subitem_index =
      (int *)HECMW_malloc(sizeof(int) * (mat->n_mat_item + 1));
  if (ref_mat->mat_subitem_index == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < mat->n_mat_item + 1; i++) {
    ref_mat->mat_subitem_index[i] = mat->mat_subitem_index[i];
  }

  /* mat_table_index */
  ref_mat->mat_table_index =
      (int *)HECMW_malloc(sizeof(int) * (mat->n_mat_subitem + 1));
  if (ref_mat->mat_table_index == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < mat->n_mat_subitem + 1; i++) {
    ref_mat->mat_table_index[i] = mat->mat_table_index[i];
  }

  /* mat_val */
  ref_mat->mat_val = (double *)HECMW_malloc(sizeof(double) * mat->n_mat_table);
  if (ref_mat->mat_val == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < mat->n_mat_table; i++) {
    ref_mat->mat_val[i] = mat->mat_val[i];
  }

  /* mat_temp */
  ref_mat->mat_temp = (double *)HECMW_malloc(sizeof(double) * mat->n_mat_table);
  if (ref_mat->mat_temp == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < mat->n_mat_table; i++) {
    ref_mat->mat_temp[i] = mat->mat_temp[i];
  }

  return HECMW_SUCCESS;
}

static int copy_mpc_info(const struct hecmwST_mpc *mpc,
                         struct hecmwST_mpc *ref_mpc) {
  int i;

  ref_mpc->n_mpc = mpc->n_mpc;

  if (mpc->n_mpc == 0) {
    ref_mpc->mpc_item  = NULL;
    ref_mpc->mpc_dof   = NULL;
    ref_mpc->mpc_val   = NULL;
    ref_mpc->mpc_const = NULL;
    ref_mpc->mpc_index = HECMW_malloc(sizeof(int));
    if (ref_mpc->mpc_index == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    ref_mpc->mpc_index[0] = 0;
    return HECMW_SUCCESS;
  }

#if 0
  HECMW_log(HECMW_LOG_ERROR, "refinement of MPC information is not supported\n");
  return HECMW_ERROR;

#else
  HECMW_log(HECMW_LOG_WARN, "MPC information is not refined\n");

  /* mpc_index */
  ref_mpc->mpc_index = (int *) HECMW_malloc( sizeof(int) * (mpc->n_mpc+1) );
  if( ref_mpc->mpc_index == NULL ) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for ( i = 0; i < mpc->n_mpc + 1; i++ ) {
    ref_mpc->mpc_index[i] = mpc->mpc_index[i];
  }

  /* mpc_item */
  ref_mpc->mpc_item = (int *) HECMW_malloc( sizeof(int) * mpc->mpc_index[mpc->n_mpc] );
  if( ref_mpc->mpc_item == NULL ) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for ( i = 0; i < mpc->mpc_index[mpc->n_mpc]; i++ ) {
    ref_mpc->mpc_item[i] = mpc->mpc_item[i];
  }

  /* mpc_dof */
  ref_mpc->mpc_dof = (int *) HECMW_malloc( sizeof(int) * mpc->mpc_index[mpc->n_mpc] );
  if( ref_mpc->mpc_dof == NULL ) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for ( i = 0; i < mpc->mpc_index[mpc->n_mpc]; i++ ) {
    ref_mpc->mpc_dof[i] = mpc->mpc_dof[i];
  }

  /* mpc_val */
  ref_mpc->mpc_val = (double *) HECMW_malloc( sizeof(double) * mpc->mpc_index[mpc->n_mpc] );
  if( ref_mpc->mpc_val == NULL ) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for ( i = 0; i < mpc->mpc_index[mpc->n_mpc]; i++ ) {
    ref_mpc->mpc_val[i] = mpc->mpc_val[i];
  }

  /* mpc_const */
  ref_mpc->mpc_const = (double *) HECMW_malloc( sizeof(double) * mpc->mpc_index[mpc->n_mpc] );
  if( ref_mpc->mpc_const == NULL ) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for ( i = 0; i < mpc->mpc_index[mpc->n_mpc]; i++ ) {
    ref_mpc->mpc_const[i] = mpc->mpc_const[i];
  }

  return HECMW_SUCCESS;
#endif
}

static int copy_amp_info(const struct hecmwST_amplitude *amp,
                         struct hecmwST_amplitude *ref_amp) {
  int i;

  ref_amp->n_amp = amp->n_amp;

  if (amp->n_amp == 0) {
    ref_amp->amp_name            = NULL;
    ref_amp->amp_type_definition = NULL;
    ref_amp->amp_type_time       = NULL;
    ref_amp->amp_type_value      = NULL;
    ref_amp->amp_val             = NULL;
    ref_amp->amp_table           = NULL;
    ref_amp->amp_index           = (int *)HECMW_malloc(sizeof(int));
    if (ref_amp->amp_index == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    ref_amp->amp_index[0] = 0;
    return HECMW_SUCCESS;
  }

  /* amp_name */
  ref_amp->amp_name = (char **)HECMW_malloc(sizeof(char *) * amp->n_amp);
  if (ref_amp->amp_name == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < amp->n_amp; i++) {
    ref_amp->amp_name[i] = HECMW_strdup(amp->amp_name[i]);
    if (ref_amp->amp_name[i] == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
  }

  /* amp_type_definition */
  ref_amp->amp_type_definition = (int *)HECMW_malloc(sizeof(int) * amp->n_amp);
  if (ref_amp->amp_type_definition == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < amp->n_amp; i++) {
    ref_amp->amp_type_definition[i] = amp->amp_type_definition[i];
  }

  /* amp_type_time */
  ref_amp->amp_type_time = (int *)HECMW_malloc(sizeof(int) * amp->n_amp);
  if (ref_amp->amp_type_time == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < amp->n_amp; i++) {
    ref_amp->amp_type_time[i] = amp->amp_type_time[i];
  }

  /* amp_type_value */
  ref_amp->amp_type_value = (int *)HECMW_malloc(sizeof(int) * amp->n_amp);
  if (ref_amp->amp_type_value == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < amp->n_amp; i++) {
    ref_amp->amp_type_value[i] = amp->amp_type_value[i];
  }

  /* amp_index */
  ref_amp->amp_index = (int *)HECMW_malloc(sizeof(int) * (amp->n_amp + 1));
  if (ref_amp->amp_index == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < amp->n_amp + 1; i++) {
    ref_amp->amp_index[i] = amp->amp_index[i];
  }

  /* amp_val */
  ref_amp->amp_val =
      (double *)HECMW_malloc(sizeof(double) * amp->amp_index[amp->n_amp]);
  if (ref_amp->amp_val == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < amp->amp_index[amp->n_amp]; i++) {
    ref_amp->amp_val[i] = amp->amp_val[i];
  }

  /* amp_table */
  ref_amp->amp_table =
      (double *)HECMW_malloc(sizeof(double) * amp->amp_index[amp->n_amp]);
  if (ref_amp->amp_table == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < amp->amp_index[amp->n_amp]; i++) {
    ref_amp->amp_table[i] = amp->amp_table[i];
  }

  return HECMW_SUCCESS;
}

static int copy_contact_pair(const struct hecmwST_contact_pair *cp,
                             struct hecmwST_contact_pair *ref_cp) {
  int i;

  ref_cp->n_pair = cp->n_pair;

  /* name */
  ref_cp->name = (char **)HECMW_malloc(sizeof(char *) * cp->n_pair);
  if (ref_cp->name == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < cp->n_pair; i++) {
    ref_cp->name[i] = HECMW_strdup(cp->name[i]);
    if (ref_cp->name[i] == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
  }

  /* type */
  ref_cp->type = (int *)HECMW_malloc(sizeof(int) * cp->n_pair);
  if (ref_cp->type == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < cp->n_pair; i++) {
    ref_cp->type[i] = cp->type[i];
  }

  /* slave_grp_id */
  ref_cp->slave_grp_id = (int *)HECMW_malloc(sizeof(int) * cp->n_pair);
  if (ref_cp->slave_grp_id == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < cp->n_pair; i++) {
    ref_cp->slave_grp_id[i] = cp->slave_grp_id[i];
  }

  /* master_grp_id */
  ref_cp->master_grp_id = (int *)HECMW_malloc(sizeof(int) * cp->n_pair);
  if (ref_cp->master_grp_id == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  for (i = 0; i < cp->n_pair; i++) {
    ref_cp->master_grp_id[i] = cp->master_grp_id[i];
  }

  return HECMW_SUCCESS;
}

static int copy_unchanging_info(const struct hecmwST_local_mesh *mesh,
                                struct hecmwST_local_mesh *ref_mesh) {
  HECMW_log(HECMW_LOG_DEBUG, "Started copying unchanging info...\n");

  if (copy_global_info(mesh, ref_mesh) != HECMW_SUCCESS) return HECMW_ERROR;
  if (copy_node_ndof(mesh, ref_mesh) != HECMW_SUCCESS) return HECMW_ERROR;
  if (copy_elem_etype(mesh, ref_mesh) != HECMW_SUCCESS) return HECMW_ERROR;
  if (copy_comm_info(mesh, ref_mesh) != HECMW_SUCCESS) return HECMW_ERROR;
  if (copy_adapt_info(mesh, ref_mesh) != HECMW_SUCCESS) return HECMW_ERROR;
  if (copy_section_info(mesh->section, ref_mesh->section) != HECMW_SUCCESS)
    return HECMW_ERROR;
  if (copy_material_info(mesh->material, ref_mesh->material) != HECMW_SUCCESS)
    return HECMW_ERROR;
  if (copy_mpc_info(mesh->mpc, ref_mesh->mpc) != HECMW_SUCCESS)
    return HECMW_ERROR;
  if (copy_amp_info(mesh->amp, ref_mesh->amp) != HECMW_SUCCESS)
    return HECMW_ERROR;
  if (copy_contact_pair(mesh->contact_pair, ref_mesh->contact_pair) !=
      HECMW_SUCCESS)
    return HECMW_ERROR;

  HECMW_log(HECMW_LOG_DEBUG, "Finished copying unchanging info.\n");
  return HECMW_SUCCESS;
}

/*============================================================================*/
/*                                                                            */
/*  renumber nodes                                                            */
/*                                                                            */
/*============================================================================*/
static int renumber_nodes_generate_tables(struct hecmwST_local_mesh *mesh) {
  int i;
  int count_in, count_ex, count_pex;
  int *old2new, *new2old;

  if (mesh->n_node_gross == mesh->nn_internal) return HECMW_SUCCESS;

  HECMW_log(HECMW_LOG_DEBUG, "Started generating renumbering tables...\n");

  old2new = (int *)HECMW_malloc(sizeof(int) * mesh->n_node_gross);
  if (old2new == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  new2old = (int *)HECMW_malloc(sizeof(int) * mesh->n_node_gross);
  if (new2old == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  count_in  = 0;
  count_ex  = mesh->nn_internal;
  count_pex = mesh->n_node;
  for (i = 0; i < mesh->n_node_gross; i++) {
    int rank = mesh->node_ID[i * 2 + 1];
    if (rank == mesh->my_rank) {
      old2new[i]        = count_in;
      new2old[count_in] = i;
      count_in++;
    } else if (rank >= 0) {
      old2new[i]        = count_ex;
      new2old[count_ex] = i;
      count_ex++;
    } else {
      old2new[i]         = count_pex;
      new2old[count_pex] = i;
      count_pex++;
    }
  }
  mesh->node_old2new = old2new;
  mesh->node_new2old = new2old;

  HECMW_log(HECMW_LOG_DEBUG, "Finished generating renumbering tables.\n");
  return HECMW_SUCCESS;
}

static int reorder_nodes(struct hecmwST_local_mesh *mesh, int *old2new,
                         int *new2old) {
  int i;

#ifndef NDEBUG
  for (i = 0; i < mesh->n_node_gross; i++) {
    HECMW_assert(new2old[old2new[i]] == i);
  }
#endif

  /*
   * Reorder using new2old
   */
  /* node_ID */
  {
    int *new_node_ID =
        (int *)HECMW_malloc(sizeof(int) * mesh->n_node_gross * 2);
    if (new_node_ID == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    for (i = 0; i < mesh->n_node_gross; i++) {
      int old                = new2old[i];
      new_node_ID[2 * i]     = mesh->node_ID[2 * old];
      new_node_ID[2 * i + 1] = mesh->node_ID[2 * old + 1];
    }
    HECMW_free(mesh->node_ID);
    mesh->node_ID = new_node_ID;
  }
  /* global_node_ID */
  {
    int *new_global_node_ID =
        (int *)HECMW_malloc(sizeof(int) * mesh->n_node_gross);
    if (new_global_node_ID == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    for (i = 0; i < mesh->n_node_gross; i++) {
      int old               = new2old[i];
      new_global_node_ID[i] = mesh->global_node_ID[old];
    }
    HECMW_free(mesh->global_node_ID);
    mesh->global_node_ID = new_global_node_ID;
  }
  /* node */
  {
    double *new_node =
        (double *)HECMW_malloc(sizeof(double) * mesh->n_node_gross * 3);
    if (new_node == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    for (i = 0; i < mesh->n_node_gross; i++) {
      int old = new2old[i];
      int j;
      for (j = 0; j < 3; j++) new_node[3 * i + j] = mesh->node[3 * old + j];
    }
    HECMW_free(mesh->node);
    mesh->node = new_node;
  }
  /* node_init_val_index, node_init_val_item */
  if (mesh->hecmw_flag_initcon && mesh->n_node_gross > 0) {
    /* node_init_val_index */
    int *new_node_init_val_index =
        (int *)HECMW_malloc(sizeof(int) * (mesh->n_node_gross + 1));
    if (new_node_init_val_index == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    new_node_init_val_index[0] = 0;
    for (i = 0; i < mesh->n_node_gross; i++) {
      int old = new2old[i];
      int ninit =
          mesh->node_init_val_index[old + 1] - mesh->node_init_val_index[old];
      new_node_init_val_index[i + 1] = new_node_init_val_index[i] + ninit;
    }

    /* node_init_val_item */
    if (mesh->node_init_val_index[mesh->n_node_gross] != 0) {
      double *new_node_init_val_item = (double *)HECMW_malloc(
          sizeof(double) * mesh->node_init_val_index[mesh->n_node_gross]);
      if (new_node_init_val_item == NULL) {
        HECMW_set_error(errno, "");
        return HECMW_ERROR;
      }
      for (i = 0; i < mesh->n_node_gross; i++) {
        int old     = new2old[i];
        int idx_new = new_node_init_val_index[i];
        int idx_old = mesh->node_init_val_index[old];
        int ninit   = mesh->node_init_val_index[old + 1] - idx_old;
        int j;
        for (j = 0; j < ninit; j++) {
          new_node_init_val_item[idx_new + j] =
              mesh->node_init_val_item[idx_old + j];
        }
      }
      HECMW_free(mesh->node_init_val_item);
      mesh->node_init_val_item = new_node_init_val_item;
    }
    HECMW_free(mesh->node_init_val_index);
    mesh->node_init_val_index = new_node_init_val_index;
  }

  /*
   * Update using old2new
   */
  /* elem_node_item */
  for (i = 0; i < mesh->elem_node_index[mesh->n_elem_gross]; i++) {
    int old                 = mesh->elem_node_item[i];
    int new                 = old2new[old - 1] + 1;
    mesh->elem_node_item[i] = new;
  }
  /* import_item */
  for (i = 0; i < mesh->import_index[mesh->n_neighbor_pe]; i++) {
    int old              = mesh->import_item[i];
    int new              = old2new[old - 1] + 1;
    mesh->import_item[i] = new;
  }
  /* export_item */
  for (i = 0; i < mesh->export_index[mesh->n_neighbor_pe]; i++) {
    int old              = mesh->export_item[i];
    int new              = old2new[old - 1] + 1;
    mesh->export_item[i] = new;
  }
  /* mpc->mpc_item */
  for (i = 0; i < mesh->mpc->mpc_index[mesh->mpc->n_mpc]; i++) {
    int old                = mesh->mpc->mpc_item[i];
    int new                = old2new[old - 1] + 1;
    mesh->mpc->mpc_item[i] = new;
  }
  /* node_group->grp_item */
  for (i = 0; i < mesh->node_group->grp_index[mesh->node_group->n_grp]; i++) {
    int old                       = mesh->node_group->grp_item[i];
    int new                       = old2new[old - 1] + 1;
    mesh->node_group->grp_item[i] = new;
  }
  return HECMW_SUCCESS;
}

static int delete_external_items(int n_elem, int n_grp, int *index, int *item,
                                 int n_memb) {
  int i, j, k;
  int n_del = 0;
  int start = 0;
  for (i = 0; i < n_grp; i++) {
    int end = index[i + 1];
    for (j = start; j < end; j++) {
      if (item[j * n_memb] > n_elem) {
        n_del++;
      } else {
        for (k                           = 0; k < n_memb; k++)
          item[(j - n_del) * n_memb + k] = item[j * n_memb + k];
      }
    }
    /* start for next group */
    start        = end;
    index[i + 1] = end - n_del;
  }
  return HECMW_SUCCESS;
}

static int renumber_nodes(struct hecmwST_local_mesh *mesh) {
  if (mesh->n_node_gross == mesh->nn_internal) return HECMW_SUCCESS;

  HECMW_log(HECMW_LOG_DEBUG, "Started renumbering nodes...\n");

  if (reorder_nodes(mesh, mesh->node_old2new, mesh->node_new2old) !=
      HECMW_SUCCESS) {
    return HECMW_ERROR;
  }

  delete_external_items(mesh->n_node, mesh->node_group->n_grp,
                        mesh->node_group->grp_index, mesh->node_group->grp_item,
                        1);

  HECMW_log(HECMW_LOG_DEBUG, "Finished renumbering nodes.\n");
  return HECMW_SUCCESS;
}

static int renumber_back_nodes(struct hecmwST_local_mesh *mesh) {
  if (mesh->n_node_gross == mesh->nn_internal) return HECMW_SUCCESS;

  HECMW_log(HECMW_LOG_DEBUG, "Started renumbering back nodes...\n");

  if (reorder_nodes(mesh, mesh->node_new2old, mesh->node_old2new) !=
      HECMW_SUCCESS) {
    return HECMW_ERROR;
  }

  HECMW_log(HECMW_LOG_DEBUG, "Finished renumbering back nodes.\n");
  return HECMW_SUCCESS;
}

/*============================================================================*/
/*                                                                            */
/*  renumber elements                                                         */
/*                                                                            */
/*============================================================================*/
static int renumber_elements_generate_tables(struct hecmwST_local_mesh *mesh) {
  int i;
  int count_rel, count_pex;
  int *old2new, *new2old;

  if (mesh->n_node_gross == mesh->nn_internal) return HECMW_SUCCESS;

  HECMW_log(HECMW_LOG_DEBUG, "Started generating renumbering tables...\n");

  old2new = (int *)HECMW_malloc(sizeof(int) * mesh->n_elem_gross);
  if (old2new == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  new2old = (int *)HECMW_malloc(sizeof(int) * mesh->n_elem_gross);
  if (new2old == NULL) {
    HECMW_set_error(errno, "");
    return HECMW_ERROR;
  }
  count_rel = 0;
  count_pex = mesh->n_elem;
  for (i = 0; i < mesh->n_elem_gross; i++) {
    int rank = mesh->elem_ID[i * 2 + 1];
    if (rank >= 0) {
      old2new[i]         = count_rel;
      new2old[count_rel] = i;
      count_rel++;
    } else {
      old2new[i]         = count_pex;
      new2old[count_pex] = i;
      count_pex++;
    }
  }
  mesh->elem_old2new = old2new;
  mesh->elem_new2old = new2old;

  HECMW_log(HECMW_LOG_DEBUG, "Finished generating renumbering tables.\n");
  return HECMW_SUCCESS;
}

static int reorder_elems(struct hecmwST_local_mesh *mesh, int *old2new,
                         int *new2old) {
  int i;

#ifndef NDEBUG
  for (i = 0; i < mesh->n_elem_gross; i++) {
    HECMW_assert(new2old[old2new[i]] == i);
  }
#endif

  /*
   * Reorder using new2old
   */
  /* elem_ID */
  {
    int *new_elem_ID =
        (int *)HECMW_malloc(sizeof(int) * mesh->n_elem_gross * 2);
    if (new_elem_ID == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    for (i = 0; i < mesh->n_elem_gross; i++) {
      int old                = new2old[i];
      new_elem_ID[2 * i]     = mesh->elem_ID[2 * old];
      new_elem_ID[2 * i + 1] = mesh->elem_ID[2 * old + 1];
    }
    HECMW_free(mesh->elem_ID);
    mesh->elem_ID = new_elem_ID;
  }
  /* global_elem_ID */
  {
    int *new_global_elem_ID =
        (int *)HECMW_malloc(sizeof(int) * mesh->n_elem_gross);
    if (new_global_elem_ID == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    for (i = 0; i < mesh->n_elem_gross; i++) {
      int old               = new2old[i];
      new_global_elem_ID[i] = mesh->global_elem_ID[old];
    }
    HECMW_free(mesh->global_elem_ID);
    mesh->global_elem_ID = new_global_elem_ID;
  }
  /* elem_type */
  {
    int *new_elem_type = (int *)HECMW_malloc(sizeof(int) * mesh->n_elem_gross);
    if (new_elem_type == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    for (i = 0; i < mesh->n_elem_gross; i++) {
      int old          = new2old[i];
      new_elem_type[i] = mesh->elem_type[old];
    }
    HECMW_free(mesh->elem_type);
    mesh->elem_type = new_elem_type;
  }
  /* elem_type_index, elem_type_item */
  {
    for (i = 0; i < mesh->n_elem_type; i++) {
      int j;
      for (j = mesh->elem_type_index[i]; j < mesh->n_elem; j++) {
        if (mesh->elem_type[j] != mesh->elem_type_item[i]) {
          break;
        }
      }
      mesh->elem_type_index[i + 1] = j;
    }
  }
  /* elem_node_index, elem_node_item */
  {
    int *new_elem_node_index =
        (int *)HECMW_malloc(sizeof(int) * (mesh->n_elem_gross + 1));
    int *new_elem_node_item = (int *)HECMW_malloc(
        sizeof(int) * mesh->elem_node_index[mesh->n_elem_gross]);
    if (new_elem_node_index == NULL || new_elem_node_item == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    new_elem_node_index[0] = 0;
    for (i = 0; i < mesh->n_elem_gross; i++) {
      int old = new2old[i];
      int nn  = mesh->elem_node_index[old + 1] - mesh->elem_node_index[old];
      new_elem_node_index[i + 1] = new_elem_node_index[i] + nn;
    }
    for (i = 0; i < mesh->n_elem_gross; i++) {
      int old       = new2old[i];
      int start_old = mesh->elem_node_index[old];
      int start_new = new_elem_node_index[i];
      int nn        = new_elem_node_index[i + 1] - start_new;
      int j;
      for (j                              = 0; j < nn; j++)
        new_elem_node_item[start_new + j] = mesh->elem_node_item[start_old + j];
    }
    HECMW_free(mesh->elem_node_index);
    HECMW_free(mesh->elem_node_item);
    mesh->elem_node_index = new_elem_node_index;
    mesh->elem_node_item  = new_elem_node_item;
  }
  /* section_ID */
  {
    int *new_section_ID = (int *)HECMW_malloc(sizeof(int) * mesh->n_elem_gross);
    if (new_section_ID == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    for (i = 0; i < mesh->n_elem_gross; i++) {
      int old           = new2old[i];
      new_section_ID[i] = mesh->section_ID[old];
    }
    HECMW_free(mesh->section_ID);
    mesh->section_ID = new_section_ID;
  }
  /* elem_mat_ID_index, elem_mat_ID_item */
  {
    int *new_emID_index, *new_emID_item;
    new_emID_index =
        (int *)HECMW_malloc(sizeof(int) * (mesh->n_elem_gross + 1));
    if (new_emID_index == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    new_emID_item = (int *)HECMW_malloc(
        sizeof(int) * mesh->elem_mat_ID_index[mesh->n_elem_gross]);
    if (new_emID_item == NULL) {
      HECMW_set_error(errno, "");
      return HECMW_ERROR;
    }
    new_emID_index[0] = 0;
    for (i = 0; i < mesh->n_elem_gross; i++) {
      int old, js_old, je_old, len, js_new, j;
      old                   = new2old[i];
      js_old                = mesh->elem_mat_ID_index[old];
      je_old                = mesh->elem_mat_ID_index[old + 1];
      len                   = je_old - js_old;
      js_new                = new_emID_index[i];
      new_emID_index[i + 1] = js_new + len;
      for (j = 0; j < len; j++) {
        new_emID_item[js_new + j] = mesh->elem_mat_ID_item[js_old + j];
      }
    }
    HECMW_free(mesh->elem_mat_ID_index);
    HECMW_free(mesh->elem_mat_ID_item);
    mesh->elem_mat_ID_index = new_emID_index;
    mesh->elem_mat_ID_item  = new_emID_item;
  }

  /*
   * Update using old2new
   */
  /* elem_internal_list */
  for (i = 0; i < mesh->ne_internal; i++) {
    int old                     = mesh->elem_internal_list[i];
    int new                     = old2new[old - 1] + 1;
    mesh->elem_internal_list[i] = new;
  }
  /* shared_item */
  for (i = 0; i < mesh->shared_index[mesh->n_neighbor_pe]; i++) {
    int old              = mesh->shared_item[i];
    int new              = old2new[old - 1] + 1;
    mesh->shared_item[i] = new;
  }
  /* elem_groups->grp_item */
  for (i = 0; i < mesh->elem_group->grp_index[mesh->elem_group->n_grp]; i++) {
    int old                       = mesh->elem_group->grp_item[i];
    int new                       = old2new[old - 1] + 1;
    mesh->elem_group->grp_item[i] = new;
  }
  /* surf_groups->grp_item */
  for (i = 0; i < mesh->surf_group->grp_index[mesh->surf_group->n_grp]; i++) {
    int old                           = mesh->surf_group->grp_item[2 * i];
    int new                           = old2new[old - 1] + 1;
    mesh->surf_group->grp_item[2 * i] = new;
  }
  return HECMW_SUCCESS;
}

static int renumber_elements(struct hecmwST_local_mesh *mesh) {
  if (mesh->n_elem_gross == mesh->n_elem) return HECMW_SUCCESS;

  HECMW_log(HECMW_LOG_DEBUG, "Started renumbering elements...\n");

  if (reorder_elems(mesh, mesh->elem_old2new, mesh->elem_new2old) !=
      HECMW_SUCCESS) {
    return HECMW_ERROR;
  }

  /* delete_external_items( mesh->n_elem, mesh->n_neighbor_pe,
   * mesh->shared_index, mesh->shared_item, 1 ); */
  delete_external_items(mesh->n_elem, mesh->elem_group->n_grp,
                        mesh->elem_group->grp_index, mesh->elem_group->grp_item,
                        1);
  delete_external_items(mesh->n_elem, mesh->surf_group->n_grp,
                        mesh->surf_group->grp_index, mesh->surf_group->grp_item,
                        2);

  HECMW_log(HECMW_LOG_DEBUG, "Finished renumbering elements.\n");
  return HECMW_SUCCESS;
}

static int renumber_back_elements(struct hecmwST_local_mesh *mesh) {
  if (mesh->n_elem_gross == mesh->n_elem) return HECMW_SUCCESS;

  HECMW_log(HECMW_LOG_DEBUG, "Started renumbering back elements...\n");

  if (reorder_elems(mesh, mesh->elem_new2old, mesh->elem_old2new) !=
      HECMW_SUCCESS) {
    return HECMW_ERROR;
  }

  HECMW_log(HECMW_LOG_DEBUG, "Finished renumbering back elements.\n");
  return HECMW_SUCCESS;
}

/*============================================================================*/
/*                                                                            */
/*  refine distributed mesh ONCE                                              */
/*                                                                            */
/*============================================================================*/
static int refine_dist_mesh(struct hecmwST_local_mesh *mesh,
                            struct hecmwST_local_mesh *ref_mesh) {
  HECMW_log(HECMW_LOG_DEBUG, "Started refining mesh...\n");

  if (copy_unchanging_info(mesh, ref_mesh) != HECMW_SUCCESS) return HECMW_ERROR;
  if (call_refiner(mesh, ref_mesh) != HECMW_SUCCESS) return HECMW_ERROR;
  if (rebuild_info(mesh, ref_mesh) != HECMW_SUCCESS) return HECMW_ERROR;

  HECMW_log(HECMW_LOG_DEBUG, "Finished refining mesh.\n");
  return HECMW_SUCCESS;
}

/*============================================================================*/
/*                                                                            */
/*  refine HEC-MW distributed mesh data                                       */
/*                                                                            */
/*============================================================================*/
int HECMW_dist_refine(struct hecmwST_local_mesh **mesh, int refine,
                      const char *cad_filename, const char *part_filename) {
  int error_flag = HECMW_SUCCESS;
  int i;

  if (refine <= 0) return HECMW_SUCCESS;

  HECMW_log(HECMW_LOG_DEBUG, "Refinement requested; starting...\n");

  if ((*mesh)->n_refine > 0) {
    if (renumber_back_elements(*mesh) != HECMW_SUCCESS) return HECMW_ERROR;
    if (renumber_back_nodes(*mesh) != HECMW_SUCCESS) return HECMW_ERROR;
  }

  if (prepare_refiner(*mesh, cad_filename, part_filename) != HECMW_SUCCESS)
    return HECMW_ERROR;

  for (i = 0; i < refine; i++) {
    struct hecmwST_local_mesh *ref_mesh = HECMW_dist_alloc();
    if (ref_mesh == NULL) {
      error_flag = HECMW_ERROR;
      break;
    }

    if (i > 0) {
      clear_refiner();
    }

    HECMW_log(HECMW_LOG_DEBUG, "Refining(%d)...\n", i + 1);

    if (refine_dist_mesh(*mesh, ref_mesh) != HECMW_SUCCESS) {
      HECMW_log(HECMW_LOG_ERROR, "Refinement failed\n");
      error_flag = HECMW_ERROR;
      break;
    }

    HECMW_dist_free(*mesh);
    *mesh = ref_mesh;
  }

  if (terminate_refiner(*mesh) != HECMW_SUCCESS) return HECMW_ERROR;

  if (error_flag != HECMW_SUCCESS) return HECMW_ERROR;

  if (renumber_nodes(*mesh) != HECMW_SUCCESS) return HECMW_ERROR;
  if (renumber_elements(*mesh) != HECMW_SUCCESS) return HECMW_ERROR;

  HECMW_log(HECMW_LOG_DEBUG, "Refinement finished.\n");

  return HECMW_SUCCESS;
}

#else  /* HECMW_WITH_REFINER */

int HECMW_dist_refine(struct hecmwST_local_mesh **mesh, int refine,
                      const char *cad_filename, const char *part_filename) {
  if (refine > 0) {
    HECMW_log(HECMW_LOG_WARN, "Refiner not enabled; ignoring...\n");
  }
  return HECMW_SUCCESS;
}
#endif /* HECMW_WITH_REFINER */
