/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include "hecmw_util.h"
#include "hecmw_config.h"
#include "hecmw_bin_io.h"
#include "hecmw_result.h"
#include "hecmw_result_io.h"
#include "hecmw_result_io_bin.h"
#include "hecmw_result_io_txt.h"

struct fortran_remainder {
  double *ptr;
  struct fortran_remainder *next;
};

static struct fortran_remainder *remainder; /* for Fortran */

void HECMW_result_free(struct hecmwST_result_data *result) {
  int i;

  if (result == NULL) return;

  if (result->ng_component > 0) {
    HECMW_free(result->ng_dof);
    HECMW_free(result->global_val_item);
    for (i = 0; i < result->ng_component; i++) {
      HECMW_free(result->global_label[i]);
    }
    HECMW_free(result->global_label);
  }

  if (result->nn_component > 0) {
    HECMW_free(result->nn_dof);
    HECMW_free(result->node_val_item);
    for (i = 0; i < result->nn_component; i++) {
      HECMW_free(result->node_label[i]);
    }
    HECMW_free(result->node_label);
  }

  if (result->ne_component > 0) {
    HECMW_free(result->ne_dof);
    HECMW_free(result->elem_val_item);
    for (i = 0; i < result->ne_component; i++) {
      HECMW_free(result->elem_label[i]);
    }
    HECMW_free(result->elem_label);
  }

  HECMW_free(result);
}

int HECMW_result_init(struct hecmwST_local_mesh *hecMESH,
                      int i_step, char *header, char *comment) {
  return HECMW_result_io_init(
      hecMESH->n_node, hecMESH->n_elem, hecMESH->global_node_ID,
      hecMESH->global_elem_ID, i_step, header, comment);
}

int HECMW_result_finalize(void) {
  HECMW_result_io_finalize();
  return 0;
}


/*---------------------------------------------------------------------------*/
/* UNIVERSAL I/O                                                             */
/*---------------------------------------------------------------------------*/

int HECMW_result_write_by_name(char *name_ID) {
  char *basename, filename[HECMW_FILENAME_LEN + 1];
  int fg_text, ret;

  if ((basename =
           HECMW_ctrl_get_result_file(name_ID, istep, &fg_text)) == NULL)
    return -1;

  ret = snprintf(filename, HECMW_FILENAME_LEN + 1, "%s.%d", basename, istep);
  HECMW_free(basename);
  if (ret > HECMW_FILENAME_LEN) return -1;

  if (fg_text) {
    if (HECMW_result_io_txt_write_by_fname(filename)) return -1;
  } else {
    if (HECMW_result_io_bin_write_by_fname(filename)) return -1;
  }

  return 0;
}

int HECMW_result_write_ST_by_name(char *name_ID,
                                  struct hecmwST_result_data *result,
                                  int n_node, int n_elem, char *header, char *comment) {
  char *basename, filename[HECMW_FILENAME_LEN + 1];
  int fg_text, ret;

  if ((basename =
           HECMW_ctrl_get_result_file(name_ID, istep, &fg_text)) == NULL)
    return -1;

  ret = snprintf(filename, HECMW_FILENAME_LEN + 1, "%s.%d", basename, istep);
  HECMW_free(basename);
  if (ret > HECMW_FILENAME_LEN) return -1;

  if (fg_text) {
    if (HECMW_result_io_txt_write_ST_by_fname(filename, result, n_node, n_elem,
                                           header, comment))
      return -1;
  } else {
    if (HECMW_result_io_bin_write_ST_by_fname(filename, result, n_node, n_elem,
                                           header, comment))
      return -1;
  }

  return 0;
}

int HECMW_result_write_by_addfname(char *name_ID, char *addfname) {
  char *basename, filename[HECMW_FILENAME_LEN + 1];
  int fg_text, myrank, ret;

  if ((basename = HECMW_ctrl_get_result_fileheader(name_ID, istep,
                                                   &fg_text)) == NULL)
    return -1;

  myrank = HECMW_comm_get_rank();
  ret    = snprintf(filename, HECMW_FILENAME_LEN + 1, "%s%s.%d.%d", basename,
                 addfname, myrank, istep);
  HECMW_free(basename);
  if (ret > HECMW_FILENAME_LEN) return -1;

  if (fg_text) {
    if (HECMW_result_io_txt_write_by_fname(filename)) return -1;
  } else {
    if (HECMW_result_io_bin_write_by_fname(filename)) return -1;
  }

  return 0;
}

int HECMW_result_checkfile_by_name(char *name_ID, int i_step) {
  char *basename, filename[HECMW_FILENAME_LEN + 1];
  int fg_text, ret;
  FILE *fp;

  if ((basename = HECMW_ctrl_get_result_file(name_ID, i_step,
                                             &fg_text)) == NULL)
    return -1;

  ret = snprintf(filename, HECMW_FILENAME_LEN + 1, "%s.%d", basename, i_step);
  HECMW_free(basename);
  if (ret > HECMW_FILENAME_LEN) return -1;

  fp = fopen(filename, "r");
  if (fp == NULL) return -1;
  fclose(fp);

  return 0;
}

struct hecmwST_result_data *HECMW_result_read_by_fname(char *filename) {
  struct hecmwST_result_data *result;

  if (HECMW_result_io_bin_judge_file(filename)) {
    result = HECMW_result_io_bin_read_by_fname(filename);
  } else {
    result = HECMW_result_io_txt_read_by_fname(filename);
  }

  return result;
}

struct hecmwST_result_data *HECMW_result_read_by_name(char *name_ID,
                                                      int i_step) {
  char *basename, filename[HECMW_FILENAME_LEN + 1];
  struct hecmwST_result_data *result;
  int fg_text, ret;

  if ((basename = HECMW_ctrl_get_result_file(name_ID, i_step,
                                             &fg_text)) == NULL)
    return NULL;

  ret = snprintf(filename, HECMW_FILENAME_LEN + 1, "%s.%d", basename, i_step);
  HECMW_free(basename);
  if (ret > HECMW_FILENAME_LEN) return NULL;

  if ((result = HECMW_result_read_by_fname(filename)) == NULL) return NULL;

  return result;
}

/*---------------------------------------------------------------------------*/
/* etc.                                                                      */
/*---------------------------------------------------------------------------*/

int HECMW_result_get_nnode(void) { return nnode; }

int HECMW_result_get_nelem(void) { return nelem; }

char *HECMW_result_get_header(char *buff) {
  strcpy(buff, head);
  return buff;
}

char *HECMW_result_get_comment(char *buff) {
  strcpy(buff, comment_line);
  return buff;
}

int *HECMW_result_get_nodeID(int *buff) {
  int i;
  for (i = 0; i < nnode; i++) {
    buff[i] = node_global_ID[i];
  }
  return buff;
}

int *HECMW_result_get_elemID(int *buff) {
  int i;
  for (i = 0; i < nelem; i++) {
    buff[i] = elem_global_ID[i];
  }
  return buff;
}

void HECMW_result_free_nodeID(void) {
  HECMW_free(node_global_ID);
  node_global_ID = NULL;
}

void HECMW_result_free_elemID(void) {
  HECMW_free(elem_global_ID);
  elem_global_ID = NULL;
}

/*---------------------------------------------------------------------------*/
/* FORTRAN INTERFACE                                                         */
/*---------------------------------------------------------------------------*/

void hecmw_result_init_if(int *n_node, int *n_elem, int *nodeID, int *elemID,
                          int *i_step, char *header, char *comment, int *err,
                          int len) {
  char header_str[HECMW_HEADER_LEN + 1];
  char comment_str[HECMW_MSG_LEN + 1];

  *err = 1;
  if (HECMW_strcpy_f2c_r(header, len, header_str, sizeof(header_str)) == NULL)
    return;
  if (HECMW_strcpy_f2c_r(comment, len, comment_str, sizeof(comment_str)) == NULL)
    return;
  if (HECMW_result_io_init(*n_node, *n_elem, nodeID, elemID, *i_step,
                             header_str, comment_str))
    return;
  *err = 0;
}

void hecmw_result_init_if_(int *n_node, int *n_elem, int *nodeID, int *elemID,
                           int *i_step, char *header, char *comment, int *err,
                           int len) {
  hecmw_result_init_if(n_node, n_elem, nodeID, elemID, i_step, header, comment,
                       err, len);
}

void hecmw_result_init_if__(int *n_node, int *n_elem, int *nodeID, int *elemID,
                            int *i_step, char *header, char *comment, int *err,
                            int len) {
  hecmw_result_init_if(n_node, n_elem, nodeID, elemID, i_step, header, comment,
                       err, len);
}

void HECMW_RESULT_INIT_IF(int *n_node, int *n_elem, int *nodeID, int *elemID,
                          int *i_step, char *header, char *comment, int *err,
                          int len) {
  hecmw_result_init_if(n_node, n_elem, nodeID, elemID, i_step, header, comment,
                       err, len);
}

/*---------------------------------------------------------------------------*/

void hecmw_result_finalize_if(int *err) {
  *err = 1;
  if (HECMW_result_finalize()) return;
  node_global_ID = NULL;
  elem_global_ID = NULL;
  *err           = 0;
}

void hecmw_result_finalize_if_(int *err) { hecmw_result_finalize_if(err); }

void hecmw_result_finalize_if__(int *err) { hecmw_result_finalize_if(err); }

void HECMW_RESULT_FINALIZE_IF(int *err) { hecmw_result_finalize_if(err); }

/*---------------------------------------------------------------------------*/

void hecmw_result_add_if(int *dtype, int *n_dof, char *label,
                         double *ptr, int *err, int len) {
  char label_str[HECMW_NAME_LEN + 1];
  int n, size;
  double *data;
  struct fortran_remainder *remain;

  *err = 1;

  if (HECMW_strcpy_f2c_r(label, len, label_str, sizeof(label_str)) == NULL)
    return;

  if (*dtype == 1) { //node
    n = nnode;
  } else if (*dtype == 2) { //element
    n = nelem;
  } else { // global
    n = 1;
  }
  size = sizeof(double) * n * (*n_dof);
  data = HECMW_malloc(size);
  if (data == NULL) {
    HECMW_set_error(errno, "");
    return;
  }
  memcpy(data, ptr, size);

  remain = HECMW_malloc(sizeof(*remain));
  if (remain == NULL) {
    HECMW_set_error(errno, "");
    return;
  }
  remain->ptr  = data;
  remain->next = remainder;
  remainder    = remain;

  if (HECMW_result_io_add(*dtype, *n_dof, label_str, data)) return;

  *err = 0;
}

void hecmw_result_add_if_(int *dtype, int *n_dof, char *label,
                          double *ptr, int *err, int len) {
  hecmw_result_add_if(dtype, n_dof, label, ptr, err, len);
}

void hecmw_result_add_if__(int *dtype, int *n_dof, char *label,
                           double *ptr, int *err, int len) {
  hecmw_result_add_if(dtype, n_dof, label, ptr, err, len);
}

void HECMW_RESULT_ADD_IF(int *dtype, int *n_dof, char *label,
                         double *ptr, int *err, int len) {
  hecmw_result_add_if(dtype, n_dof, label, ptr, err, len);
}

/*----------------------------------------------------------------------------*/

void hecmw_result_write_by_name_if(char *name_ID, int *err, int len) {
  char name_ID_str[HECMW_NAME_LEN + 1];
  struct fortran_remainder *p, *q;

  *err = 1;

  if (HECMW_strcpy_f2c_r(name_ID, len, name_ID_str, sizeof(name_ID_str)) ==
      NULL)
    return;

  if (HECMW_result_write_by_name(name_ID_str)) return;

  for (p = remainder; p; p = q) {
    q = p->next;
    HECMW_free(p->ptr);
    HECMW_free(p);
  }
  remainder = NULL;

  *err = 0;
}

void hecmw_result_write_by_name_if_(char *name_ID, int *err, int len) {
  hecmw_result_write_by_name_if(name_ID, err, len);
}

void hecmw_result_write_by_name_if__(char *name_ID, int *err, int len) {
  hecmw_result_write_by_name_if(name_ID, err, len);
}

void HECMW_RESULT_WRITE_BY_NAME_IF(char *name_ID, int *err, int len) {
  hecmw_result_write_by_name_if(name_ID, err, len);
}

/*---------------------------------------------------------------------------*/

void hecmw_result_write_by_addfname_if(char *name_ID, char *addfname, int *err,
                                       int len1, int len2) {
  char name_ID_str[HECMW_NAME_LEN + 1];
  char addfname_str[HECMW_NAME_LEN + 1];
  struct fortran_remainder *p, *q;

  *err = 1;

  if (HECMW_strcpy_f2c_r(name_ID, len1, name_ID_str, sizeof(name_ID_str)) ==
      NULL)
    return;
  if (HECMW_strcpy_f2c_r(addfname, len2, addfname_str, sizeof(name_ID_str)) ==
      NULL)
    return;

  if (HECMW_result_write_by_addfname(name_ID_str, addfname_str)) return;

  for (p = remainder; p; p = q) {
    q = p->next;
    HECMW_free(p->ptr);
    HECMW_free(p);
  }
  remainder = NULL;

  *err = 0;
}

void hecmw_result_write_by_addfname_if_(char *name_ID, char *addfname, int *err,
                                        int len1, int len2) {
  hecmw_result_write_by_addfname_if(name_ID, addfname, err, len1, len2);
}

void hecmw_result_write_by_addfname_if__(char *name_ID, char *addfname,
                                         int *err, int len1, int len2) {
  hecmw_result_write_by_addfname_if(name_ID, addfname, err, len1, len2);
}

void HECMW_RESULT_WRITE_BY_ADDFNAME_IF(char *name_ID, char *addfname, int *err,
                                       int len1, int len2) {
  hecmw_result_write_by_addfname_if(name_ID, addfname, err, len1, len2);
}

/*----------------------------------------------------------------------------*/

void hecmw_result_checkfile_by_name_if(char *name_ID, int *i_step, int *err, int len) {
  char name_ID_str[HECMW_NAME_LEN + 1];

  *err = 1;

  if (HECMW_strcpy_f2c_r(name_ID, len, name_ID_str, sizeof(name_ID_str)) ==
      NULL)
    return;

  if (HECMW_result_checkfile_by_name(name_ID_str, *i_step)) return;

  *err = 0;
}

void hecmw_result_checkfile_by_name_if_(char *name_ID, int *i_step, int *err, int len) {
  hecmw_result_checkfile_by_name_if(name_ID, i_step, err, len);
}

void hecmw_result_checkfile_by_name_if__(char *name_ID, int *i_step, int *err, int len) {
  hecmw_result_checkfile_by_name_if(name_ID, i_step, err, len);
}

void HECMW_RESULT_CHECKFILE_BY_NAME_IF(char *name_ID, int *i_step, int *err, int len) {
  hecmw_result_checkfile_by_name_if(name_ID, i_step, err, len);
}
