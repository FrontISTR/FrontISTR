/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "hecmw_struct.h"
#include "hecmw_util.h"
#include "hecmw_result.h"

static struct hecmwST_result_data *Result;
static int NNode, NElem;

/*-----------------------------------------------------------------------------
 * SetFunc
 */

static int set_ng_component(void *dst) {
  *((int *)dst) = Result->ng_component;
  return 0;
}

static int set_ng_dof(void *dst) {
  void *src;
  int size;

  if (Result->ng_component <= 0) return 0;

  src  = Result->ng_dof;
  size = sizeof(*Result->ng_dof) * Result->ng_component;
  memcpy(dst, src, size);

  return 0;
}

static int set_global_label(void *dst) {
  int i;

  if (Result->ng_component <= 0) return 0;

  for (i = 0; i < Result->ng_component; i++) {
    char *dst_point = (char *)dst + HECMW_NAME_LEN * i;
    char *src       = Result->global_label[i];
    HECMW_strcpy_c2f(src, dst_point, HECMW_NAME_LEN);
  }

  return 0;
}

static int set_global_val_item(void *dst) {
  void *src;
  int i, n, size;

  if (Result->ng_component <= 0) return 0;

  n = 0;
  for (i = 0; i < Result->ng_component; i++) {
    n += Result->ng_dof[i];
  }
  src  = Result->global_val_item;
  size = sizeof(*Result->global_val_item) * n;
  memcpy(dst, src, size);

  return 0;
}

static int set_nn_component(void *dst) {
  *((int *)dst) = Result->nn_component;
  return 0;
}

static int set_nn_dof(void *dst) {
  void *src;
  int size;

  if (Result->nn_component <= 0) return 0;

  src  = Result->nn_dof;
  size = sizeof(*Result->nn_dof) * Result->nn_component;
  memcpy(dst, src, size);

  return 0;
}

static int set_node_label(void *dst) {
  int i;

  if (Result->nn_component <= 0) return 0;

  for (i = 0; i < Result->nn_component; i++) {
    char *dst_point = (char *)dst + HECMW_NAME_LEN * i;
    char *src       = Result->node_label[i];
    HECMW_strcpy_c2f(src, dst_point, HECMW_NAME_LEN);
  }

  return 0;
}

static int set_node_val_item(void *dst) {
  void *src;
  int i, n, size;

  if (Result->nn_component <= 0) return 0;

  n = 0;
  for (i = 0; i < Result->nn_component; i++) {
    n += Result->nn_dof[i];
  }
  src  = Result->node_val_item;
  size = sizeof(*Result->node_val_item) * n * NNode;
  memcpy(dst, src, size);

  return 0;
}

static int set_ne_component(void *dst) {
  *((int *)dst) = Result->ne_component;
  return 0;
}

static int set_ne_dof(void *dst) {
  void *src;
  int size;

  if (Result->ne_component <= 0) return 0;

  src  = Result->ne_dof;
  size = sizeof(*Result->ne_dof) * Result->ne_component;
  memcpy(dst, src, size);

  return 0;
}

static int set_elem_label(void *dst) {
  int i;

  if (Result->ne_component <= 0) return 0;

  for (i = 0; i < Result->ne_component; i++) {
    char *dst_point = (char *)dst + HECMW_NAME_LEN * i;
    char *src       = Result->elem_label[i];
    HECMW_strcpy_c2f(src, dst_point, HECMW_NAME_LEN);
  }

  return 0;
}

static int set_elem_val_item(void *dst) {
  void *src;
  int i, n, size;

  if (Result->ne_component <= 0) return 0;

  n = 0;
  for (i = 0; i < Result->ne_component; i++) {
    n += Result->ne_dof[i];
  }
  src  = Result->elem_val_item;
  size = sizeof(*Result->elem_val_item) * n * NElem;
  memcpy(dst, src, size);

  return 0;
}

/*-----------------------------------------------------------------------------
 * SetFunc table
 */

typedef int (*SetFunc)(void *);

static struct func_table {
  char *struct_name;
  char *var_name;
  SetFunc set_func;
} functions[] =
    {
        /*  { Struct name, Variable name, memcpy function, check allocation
           function }*/
        {"hecmwST_result_data", "ng_component", set_ng_component},
        {"hecmwST_result_data", "ng_dof", set_ng_dof},
        {"hecmwST_result_data", "global_label", set_global_label},
        {"hecmwST_result_data", "global_val_item", set_global_val_item},

        {"hecmwST_result_data", "nn_component", set_nn_component},
        {"hecmwST_result_data", "nn_dof", set_nn_dof},
        {"hecmwST_result_data", "node_label", set_node_label},
        {"hecmwST_result_data", "node_val_item", set_node_val_item},

        {"hecmwST_result_data", "ne_component", set_ne_component},
        {"hecmwST_result_data", "ne_dof", set_ne_dof},
        {"hecmwST_result_data", "elem_label", set_elem_label},
        {"hecmwST_result_data", "elem_val_item", set_elem_val_item},
};

static const int NFUNC = sizeof(functions) / sizeof(functions[0]);

static SetFunc get_set_func(char *struct_name, char *var_name) {
  int i;

  for (i = 0; i < NFUNC; i++) {
    if (strcmp(functions[i].struct_name, struct_name) == 0 &&
        strcmp(functions[i].var_name, var_name) == 0) {
      return functions[i].set_func;
    }
  }
  return NULL;
}

/*----------------------------------------------------------------------------*/

int HECMW_result_copy_c2f_init(struct hecmwST_result_data *result_data,
                               int n_node, int n_elem) {
  Result = result_data;
  NNode  = n_node;
  NElem  = n_elem;
  return 0;
}

int HECMW_result_copy_c2f_finalize(void) {
  Result = NULL;
  return 0;
}

/*----------------------------------------------------------------------------*/

void hecmw_result_copy_c2f_set_if(char *struct_name, char *var_name, void *dst,
                                  int *err, int len_struct, int len_var) {
  SetFunc func;
  char sname[HECMW_NAME_LEN + 1];
  char vname[HECMW_NAME_LEN + 1];

  *err = 1;

  if (Result == NULL) {
    HECMW_set_error(
        HECMW_ALL_E0102,
        "hecmw_result_copy_c2f_set_if(): 'result' has not initialized yet");
    return;
  }
  if (struct_name == NULL) {
    HECMW_set_error(HECMW_ALL_E0101,
                    "hecmw_result_copy_c2f_set_if(): 'sname' is NULL");
    return;
  }
  if (var_name == NULL) {
    HECMW_set_error(HECMW_ALL_E0101,
                    "hecmw_result_copy_c2f_set_if(): 'vname' is NULL");
    return;
  }
  if (dst == NULL) {
    HECMW_set_error(HECMW_ALL_E0101,
                    "hecmw_result_copy_c2f_set_if(): 'dst' is NULL");
    return;
  }

  if (HECMW_strcpy_f2c_r(struct_name, len_struct, sname, sizeof(sname)) ==
      NULL) {
    return;
  }
  if (HECMW_strcpy_f2c_r(var_name, len_var, vname, sizeof(vname)) == NULL) {
    return;
  }

  func = get_set_func(sname, vname);
  if (func == NULL) {
    HECMW_set_error(HECMW_ALL_E0102,
                    "hecmw_result_copy_c2f_set_if(): SetFunc not found");
    return;
  }

  if ((*func)(dst)) {
    return;
  }

  *err = 0;
}

void hecmw_result_copy_c2f_set_if_(char *struct_name, char *var_name, void *dst,
                                   int *err, int len_struct, int len_var) {
  hecmw_result_copy_c2f_set_if(struct_name, var_name, dst, err, len_struct,
                               len_var);
}

void hecmw_result_copy_c2f_set_if__(char *struct_name, char *var_name,
                                    void *dst, int *err, int len_struct,
                                    int len_var) {
  hecmw_result_copy_c2f_set_if(struct_name, var_name, dst, err, len_struct,
                               len_var);
}

void HECMW_RESULT_COPY_C2F_SET_IF(char *struct_name, char *var_name, void *dst,
                                  int *err, int len_struct, int len_var) {
  hecmw_result_copy_c2f_set_if(struct_name, var_name, dst, err, len_struct,
                               len_var);
}

/*----------------------------------------------------------------------------*/

void hecmw_result_read_by_name_if(char *name_ID, int *i_step,
                                  int *n_node, int *n_elem, int *err, int len) {
  char name_ID_str[HECMW_NAME_LEN + 1];

  *err = 1;

  if (HECMW_strcpy_f2c_r(name_ID, len, name_ID_str, sizeof(name_ID_str)) ==
      NULL)
    return;

  Result = HECMW_result_read_by_name(name_ID_str, *i_step);
  if (Result == NULL) return;

  NNode   = HECMW_result_get_nnode();
  NElem   = HECMW_result_get_nelem();
  *n_node = NNode;
  *n_elem = NElem;

  *err = 0;
}

void hecmw_result_read_by_name_if_(char *name_ID, int *i_step,
                                   int *n_node, int *n_elem, int *err,
                                   int len) {
  hecmw_result_read_by_name_if(name_ID, i_step, n_node, n_elem, err,
                               len);
}

void hecmw_result_read_by_name_if__(char *name_ID, int *i_step,
                                    int *n_node, int *n_elem, int *err,
                                    int len) {
  hecmw_result_read_by_name_if(name_ID, i_step, n_node, n_elem, err,
                               len);
}

void HECMW_RESULT_READ_BY_NAME_IF(char *name_ID, int *i_step,
                                  int *n_node, int *n_elem, int *err, int len) {
  hecmw_result_read_by_name_if(name_ID, i_step, n_node, n_elem, err,
                               len);
}

/*----------------------------------------------------------------------------*/

void hecmw_result_read_finalize_if(int *err) {
  *err = 1;
  HECMW_result_free(Result);
  HECMW_result_free_nodeID();
  HECMW_result_free_elemID();
  *err = 0;
}

void hecmw_result_read_finalize_if_(int *err) {
  hecmw_result_read_finalize_if(err);
}

void hecmw_result_read_finalize_if__(int *err) {
  hecmw_result_read_finalize_if(err);
}

void HECMW_RESULT_READ_FINALIZE_IF(int *err) {
  hecmw_result_read_finalize_if(err);
}
