/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

/**
 * @brief Utility for reading and processing results computed in parallel
 */

#include "fstr_rmerge_util.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "hecmw_util.h"
#include "hecmw_dist_free.h"
#include "hecmw_io_dist.h"
#include "hecmw_io_get_mesh.h"
#include "hecmw_etype.h"

static FILE* Log_FP;

/**
 * @brief Set file pointer for log output
 */
void fstr_set_log_fp(FILE *log_fp) {
  Log_FP = log_fp;
}

/**
 * @brief Log output
 */
void fstr_out_log(const char* fmt, ...) {
  va_list arg;
  va_start(arg, fmt);
  vfprintf(Log_FP, fmt, arg);
  va_end(arg);
}

/**
 * @brief Read all distributed meshes
 */

static int get_dist_fname(char* name_ID, char* fheader, int* fg_single,
                          int* refine, int nrank, int irank) {
  struct hecmw_ctrl_meshfiles* files;

  files = HECMW_ctrl_get_meshfiles_header_sub(name_ID, nrank, irank);
  if (!files) return -1;

  if (files->n_mesh == 1) {
    strcpy(fheader, files->meshfiles[0].filename);
    if (files->meshfiles[0].type == HECMW_CTRL_FTYPE_HECMW_DIST) {
      *fg_single = 0;
    } else {
      *fg_single = 1;
    }
    *refine = files->meshfiles[0].refine;
    fstr_out_log("refine number is %d\n", *refine);
  } else {
    HECMW_ctrl_free_meshfiles(files);
    return -1;
  }

  HECMW_ctrl_free_meshfiles(files);
  return 0;
}

static int get_area_n(char* fheader) {
  FILE* fp;
  char buff[HECMW_FILENAME_LEN + 1];
  int area = 0;

  while (1) {
    sprintf(buff, "%s.%d", fheader, area);
    fstr_out_log("try open : %s  ... ", buff);
    fp = fopen(buff, "r");
    if (!fp) {
      fstr_out_log("fail\n");
      fstr_out_log("area number is %d\n", area);
      return area;
    } else {
      fstr_out_log("success\n");
      fclose(fp);
    }
    area++;
  }
}

/**
 * @brief Read all distributed meshes
 */

struct hecmwST_local_mesh** fstr_get_all_local_mesh(char* name_ID,
                                                    int nrank,
                                                    int* area_number,
                                                    int* refine) {
  int i;
  int area_n;
  int fg_single;
  char fheader[HECMW_HEADER_LEN + 1];
  char fname[HECMW_FILENAME_LEN + 1];
  struct hecmwST_local_mesh** mesh;

  if (get_dist_fname(name_ID, fheader, &fg_single, refine, nrank, 0)) return NULL;

  if (fg_single) {
    fstr_out_log("mesh file type is NOT HECMW_DIST.\n");
    area_n = 1;
    fstr_out_log("area number is %d\n", area_n);
    mesh    = HECMW_malloc(area_n * sizeof(struct hecmwST_local_mesh*));
    mesh[0] = HECMW_get_mesh(name_ID);
    if (!mesh[0]) return NULL;
  } else {
    fstr_out_log("mesh file type is HECMW_DIST.\n");
    if (nrank == 0) {
      area_n = get_area_n(fheader);
    } else {
      area_n = nrank;
    }
    if (area_n == 0) return NULL;
    mesh = HECMW_malloc(area_n * sizeof(struct hecmwST_local_mesh*));
    for (i = 0; i < area_n; i++) {
      if (nrank == 0) {
        sprintf(fname, "%s.%d", fheader, i);
      } else {
        get_dist_fname(name_ID, fheader, &fg_single, refine, nrank, i);
        sprintf(fname, "%s.%d", fheader, i);
      }
      fstr_out_log("loading dist mesh from %s\n", fname);
      mesh[i] = HECMW_get_dist_mesh(fname);
      if (!mesh[i]) return NULL;
    }
  }
  *area_number = area_n;

  return mesh;
}

/**
 * @brief Delete mesh
 */

void fstr_free_mesh(struct hecmwST_local_mesh** mesh, int area_n) {
  int i;

  if (!mesh) return;

  for (i = 0; i < area_n; i++) {
    HECMW_dist_free(mesh[i]);
  }
  HECMW_free(mesh);
  return;
}

/**
 * @brief Check the number of steps (check for the existence of files)
 */

int fstr_get_step_n(char* name_ID, int nrank) {
  FILE* fp;
  int step, fg_text;
  char* fheader;
  char fname[HECMW_FILENAME_LEN + 1];

  if (nrank == 0) {
    if ((fheader = HECMW_ctrl_get_result_fileheader(name_ID, 1, &fg_text)) ==
        NULL)
      return 0;
  } else {
    if ((fheader = HECMW_ctrl_get_result_fileheader_sub(name_ID, 1, nrank, 0,
                                                        &fg_text)) == NULL)
      return 0;
  }

  step = 1;
  while (1) {
    sprintf(fname, "%s.0.%d", fheader, step);
    fstr_out_log("try open : %s  ... ", fname);
    fp = fopen(fname, "r");
    if (!fp) {
      fstr_out_log("fail\n");
      fstr_out_log("step number is %d\n", step - 1);
      return step - 1;
    } else {
      fstr_out_log("success\n");
      fclose(fp);
    }
    step++;
  }
}

/**
 * @brief Read all area data of step
 */

fstr_res_info** fstr_get_all_result(char* name_ID, int step, int area_n,
                                    int refine, int nrank) {
  char* fheader;
  char fname[HECMW_FILENAME_LEN + 1];
  fstr_res_info** res;
  struct hecmwST_result_data* data;
  int i, j, k, count, num, nnode, nelem, flag, fg_text;
  int *node_gid, *elem_gid;
  int refine_nnode = 0;

  if (nrank == 0) {
    if ((fheader = HECMW_ctrl_get_result_fileheader(name_ID, step,
                                                    &fg_text)) == NULL)
      return 0;
  }

  res = HECMW_malloc(area_n * sizeof(fstr_res_info*));
  if (!res) return NULL;

  for (i = 0; i < area_n; i++) {
    res[i] = HECMW_malloc(sizeof(fstr_res_info));
    if (!res[i]) return NULL;

    if (nrank != 0) {
      if ((fheader = HECMW_ctrl_get_result_fileheader_sub(
               name_ID, step, nrank, i, &fg_text)) == NULL)
        return 0;
    }
    sprintf(fname, "%s.%d.%d", fheader, i, step);
    data = HECMW_result_read_by_fname(fname);
    if (!data) return NULL;
    nnode    = HECMW_result_get_nnode();
    node_gid = HECMW_malloc(nnode * sizeof(int));
    if (!node_gid) return NULL;
    HECMW_result_get_nodeID(node_gid);
    nelem    = HECMW_result_get_nelem();
    elem_gid = HECMW_malloc(nelem * sizeof(int));
    if (!elem_gid) return NULL;
    if (data->ne_component) HECMW_result_get_elemID(elem_gid);

    if (refine) {
      res[i]->result = HECMW_malloc(sizeof(struct hecmwST_result_data));
      if (!res[i]->result) return NULL;

      res[i]->result->ng_component = data->ng_component;
      res[i]->result->ng_dof = HECMW_malloc(data->ng_component * sizeof(int));
      if (!res[i]->result->ng_dof) return NULL;
      res[i]->result->global_label =
          HECMW_malloc(data->ng_component * sizeof(char*));
      if (!res[i]->result->global_label) return NULL;
      for (j = 0; j < data->ng_component; j++) {
        res[i]->result->ng_dof[j]     = data->ng_dof[j];
        res[i]->result->global_label[j] = HECMW_strdup(data->global_label[j]);
      }

      res[i]->result->nn_component = data->nn_component;
      res[i]->result->nn_dof = HECMW_malloc(data->nn_component * sizeof(int));
      if (!res[i]->result->nn_dof) return NULL;
      res[i]->result->node_label =
          HECMW_malloc(data->nn_component * sizeof(char*));
      if (!res[i]->result->node_label) return NULL;
      num = 0;
      for (j = 0; j < data->nn_component; j++) {
        res[i]->result->nn_dof[j]     = data->nn_dof[j];
        res[i]->result->node_label[j] = HECMW_strdup(data->node_label[j]);
        num += data->nn_dof[j];
      }

      count = 1;
      flag  = 0;
      for (j = 1; j < nnode; j++) {
        if (flag == refine) break;
        count++;
        if (node_gid[j] > 0 && node_gid[j - 1] < 0) flag++;
        if (node_gid[j] < 0 && node_gid[j] > node_gid[j - 1]) flag++;
      }
      count--;
      fstr_out_log("\narea:%d -- refined_nn_internal:%d", i, count);
      refine_nnode += count;

      count = 0;
      for (j = 0; j < nnode; j++)
        if (node_gid[j] > 0) count++;
      res[i]->nnode_gid = count;
      res[i]->result->node_val_item =
          HECMW_malloc(num * count * sizeof(double));
      if (!res[i]->result->node_val_item) return NULL;
      count = 0;
      for (j = 0; j < nnode; j++) {
        if (node_gid[j] > 0) {
          for (k = 0; k < num; k++) {
            res[i]->result->node_val_item[count++] =
                data->node_val_item[num * j + k];
          }
        }
      }

      res[i]->result->ne_component = data->ne_component;
      res[i]->result->ne_dof = HECMW_malloc(data->ne_component * sizeof(int));
      if (!res[i]->result->ne_dof) return NULL;
      res[i]->result->elem_label =
          HECMW_malloc(data->ne_component * sizeof(char*));
      if (!res[i]->result->elem_label) return NULL;
      num = 0;
      for (j = 0; j < data->ne_component; j++) {
        res[i]->result->ne_dof[j]     = data->ne_dof[j];
        res[i]->result->elem_label[j] = HECMW_strdup(data->elem_label[j]);
        num += data->ne_dof[j];
      }

      count = 0;
      for (j = 0; j < nelem; j++)
        if (elem_gid[j] > 0) count++;
      fstr_out_log("\narea:%d -- ne_original from result:%d", i, count);
      res[i]->nelem_gid = count;
      res[i]->result->elem_val_item =
          HECMW_malloc(num * count * sizeof(double));
      if (!res[i]->result->elem_val_item) return NULL;
      count = 0;
      for (j = 0; j < nelem; j++) {
        if (elem_gid[j] > 0) {
          for (k = 0; k < num; k++) {
            res[i]->result->elem_val_item[count++] =
                data->elem_val_item[num * j + k];
          }
        }
      }

      HECMW_result_free(data);
      HECMW_result_free_nodeID();
      HECMW_result_free_elemID();

      res[i]->node_gid = HECMW_malloc(res[i]->nnode_gid * sizeof(int));
      if (!res[i]->node_gid) return NULL;
      count = 0;
      for (j = 0; j < nnode; j++) {
        if (node_gid[j] > 0) res[i]->node_gid[count++] = node_gid[j];
      }
      free(node_gid);

      res[i]->elem_gid = HECMW_malloc(res[i]->nelem_gid * sizeof(int));
      if (!res[i]->elem_gid) return NULL;
      count = 0;
      for (j = 0; j < nelem; j++) {
        if (elem_gid[j] > 0) res[i]->elem_gid[count++] = elem_gid[j];
      }
      free(elem_gid);
    } else {
      res[i]->result    = data;
      res[i]->nnode_gid = nnode;
      res[i]->node_gid  = node_gid;
      res[i]->nelem_gid = nelem;
      res[i]->elem_gid  = elem_gid;
    }
  }

  if (refine) fstr_out_log("\ntotal refined_nn_internal:%d\n", refine_nnode);

  return res;
}

/**
 * @brief Combine data in all areas of the step
 */

struct hecmwST_result_data* fstr_all_result(fstr_glt* glt, fstr_res_info** res,
                                            int refine) {
  int i, j, k, l_id, area;
  int gitem, nitem, eitem, count, irec;
  struct hecmwST_result_data* data;

  gitem = 0;
  for (i = 0; i < res[0]->result->ng_component; i++)
    gitem += res[0]->result->ng_dof[i];
  nitem = 0;
  for (i = 0; i < res[0]->result->nn_component; i++)
    nitem += res[0]->result->nn_dof[i];
  eitem = 0;
  for (i = 0; i < res[0]->result->ne_component; i++)
    eitem += res[0]->result->ne_dof[i];

  data             = HECMW_malloc(sizeof(struct hecmwST_result_data));
  data->ng_dof     = HECMW_malloc(res[0]->result->ng_component * sizeof(int));
  data->global_label = HECMW_malloc(res[0]->result->ng_component * sizeof(char*));
  data->global_val_item = HECMW_malloc(gitem * sizeof(double));
  data->nn_dof     = HECMW_malloc(res[0]->result->nn_component * sizeof(int));
  data->node_label = HECMW_malloc(res[0]->result->nn_component * sizeof(char*));
  data->node_val_item = HECMW_malloc(nitem * glt->node_n * sizeof(double));
  data->ne_dof     = HECMW_malloc(res[0]->result->ne_component * sizeof(int));
  data->elem_label = HECMW_malloc(res[0]->result->ne_component * sizeof(char*));
  data->elem_val_item = HECMW_malloc(eitem * glt->elem_n * sizeof(double));

  data->ng_component = res[0]->result->ng_component;
  for (i = 0; i < res[0]->result->ng_component; i++) {
    data->ng_dof[i]     = res[0]->result->ng_dof[i];
    data->global_label[i] = HECMW_strdup(res[0]->result->global_label[i]);
  }
  for (i = 0; i < gitem; i++) {
    data->global_val_item[i] = res[0]->result->global_val_item[i];
  }
  count = 0;
  for (i = 0; i < glt->node_n; i++) {
    l_id = glt->nrec[i].local;
    area = glt->nrec[i].area;
    irec = nitem * l_id;
    for (j = 0; j < res[0]->result->nn_component; j++) {
      for (k = 0; k < res[0]->result->nn_dof[j]; k++) {
        data->node_val_item[count++] = res[area]->result->node_val_item[irec++];
      }
    }
  }
  data->nn_component = res[0]->result->nn_component;
  for (i = 0; i < res[0]->result->nn_component; i++) {
    data->nn_dof[i]     = res[0]->result->nn_dof[i];
    data->node_label[i] = HECMW_strdup(res[0]->result->node_label[i]);
  }

  count = 0;
  for (i = 0; i < glt->elem_n; i++) {
    l_id = glt->erec[i].local;
    area = glt->erec[i].area;
    irec = eitem * l_id;
    for (j = 0; j < res[0]->result->ne_component; j++) {
      for (k = 0; k < res[0]->result->ne_dof[j]; k++) {
        if (refine) {
          data->elem_val_item[count++] = 0.0;
        } else {
          data->elem_val_item[count++] =
              res[area]->result->elem_val_item[irec++];
        }
      }
    }
  }
  data->ne_component = res[0]->result->ne_component;
  for (i = 0; i < res[0]->result->ne_component; i++) {
    data->ne_dof[i]     = res[0]->result->ne_dof[i];
    data->elem_label[i] = HECMW_strdup(res[0]->result->elem_label[i]);
  }

  return data;
}

/**
 * @biref Delete fstr_res_info
 */

void fstr_free_result(fstr_res_info** res, int area_n) {
  int i;

  if (!res) return;

  for (i = 0; i < area_n; i++) {
    HECMW_result_free(res[i]->result);
    HECMW_free(res[i]->node_gid);
    HECMW_free(res[i]->elem_gid);
    HECMW_free(res[i]);
  }
  HECMW_free(res);
  return;
}

/**
 * @brief Create table for global ID, local ID and belonging area records fstr_glt
 */

static int cmp_global_glt(const fstr_gl_rec* g1, const fstr_gl_rec* g2) {
  return (g1->global - g2->global);
}

typedef int (*cmp_func)(const void*, const void*);

fstr_glt* fstr_create_glt(struct hecmwST_local_mesh** mesh, int area_n) {
  int i, j, k, eid, count;
  int all_n, all_e;
  int area_e;
  fstr_gl_rec* nrec;
  fstr_gl_rec* erec;
  fstr_glt* glt;
  int* area_etype_list;

  all_n = 0;
  for (i = 0; i < area_n; i++) {
    fstr_out_log("area:%d -- nn_internal:%d\n", i, mesh[i]->nn_internal);
    all_n += mesh[i]->nn_internal;
  }
  fstr_out_log("total nn_internal:%d\n", all_n);

  nrec = HECMW_malloc(sizeof(fstr_gl_rec) * all_n);
  if (!nrec) return NULL;

  count = 0;
  for (i = 0; i < area_n; i++) {
    for (j = 0; j < mesh[i]->nn_internal; j++) {
      nrec[count].global = mesh[i]->global_node_ID[j];
      nrec[count].local  = j;
      nrec[count].area   = i;
      count++;
    }
  }

  qsort(nrec, all_n, sizeof(fstr_gl_rec), (cmp_func)cmp_global_glt);

  all_e = 0;
  for (i = 0; i < area_n; i++) {
    area_e = 0;
    area_etype_list = HECMW_malloc(sizeof(int) * mesh[i]->n_elem);
    for (j = 0; j < mesh[i]->n_elem_type; j++) {
      for (k = mesh[i]->elem_type_index[j]; k < mesh[i]->elem_type_index[j+1]; k++) {
         area_etype_list[k] = mesh[i]->elem_type_item[j];
      }
    }

    for (j = 0; j < mesh[i]->ne_internal; j++) {
      eid = mesh[i]->elem_internal_list[j] - 1;
      if ( HECMW_is_etype_patch(area_etype_list[eid]) ) continue;
      if ( HECMW_is_etype_link(area_etype_list[eid]) ) continue;
      area_e++;
    }

    HECMW_free(area_etype_list);
    all_e += area_e;
    fstr_out_log("area:%d -- ne_internal:%d\n", i, area_e);
  }
  fstr_out_log("total ne_internal:%d\n", all_e);

  erec = HECMW_malloc(sizeof(fstr_gl_rec) * all_e);
  if (!erec) return NULL;

  count = 0;
  for (i = 0; i < area_n; i++) {

    area_etype_list = HECMW_malloc(sizeof(int) * mesh[i]->n_elem);
    for (j = 0; j < mesh[i]->n_elem_type; j++) {
      for (k = mesh[i]->elem_type_index[j]; k < mesh[i]->elem_type_index[j+1]; k++) {
         area_etype_list[k] = mesh[i]->elem_type_item[j];
      }
    }

    for (j = 0; j < mesh[i]->ne_internal; j++) {
      eid = mesh[i]->elem_internal_list[j] - 1;
      if ( HECMW_is_etype_patch(area_etype_list[eid]) ) continue;
      if ( HECMW_is_etype_link(area_etype_list[eid]) ) continue;

      erec[count].global = mesh[i]->global_elem_ID[eid];
      erec[count].local = eid;
      erec[count].area  = i;
      count++;
    }

    HECMW_free(area_etype_list);
  }

  qsort(erec, all_e, sizeof(fstr_gl_rec), (cmp_func)cmp_global_glt);

  glt = HECMW_malloc(sizeof(fstr_glt));
  if (!glt) return NULL;
  glt->nrec   = nrec;
  glt->erec   = erec;
  glt->node_n = all_n;
  glt->elem_n = all_e;

  return glt;
}

/**
 * @brief Delete fstr_glt
 */

void fstr_free_glt(fstr_glt* glt) {
  if (!glt) return;

  HECMW_free(glt->nrec);
  HECMW_free(glt->erec);
  HECMW_free(glt);
  return;
}

/**
 * @brief Create a single region mesh
 */

struct hecmwST_local_mesh* fstr_create_glmesh(fstr_glt* glt) {
  struct hecmwST_local_mesh* mesh;
  int i;

  mesh                 = HECMW_calloc(1, sizeof(struct hecmwST_local_mesh));
  mesh->global_node_ID = HECMW_malloc(glt->node_n * sizeof(int));
  mesh->global_elem_ID = HECMW_malloc(glt->elem_n * sizeof(int));

  for (i = 0; i < glt->node_n; i++) {
    mesh->global_node_ID[i] = glt->nrec[i].global;
  }
  mesh->n_node = glt->node_n;

  for (i = 0; i < glt->elem_n; i++) {
    mesh->global_elem_ID[i] = glt->erec[i].global;
  }
  mesh->n_elem = glt->elem_n;

  return mesh;
}

/**
 * @brief Delete a single region mesh
 */

void fstr_free_glmesh(struct hecmwST_local_mesh* mesh) {
  if (!mesh) return;

  HECMW_free(mesh->global_node_ID);
  HECMW_free(mesh->global_elem_ID);
  HECMW_free(mesh);
  return;
}
