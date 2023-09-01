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
#include "hecmw_log.h"
#include "hecmw_malloc.h"
#include "hecmw_util.h"
#include "hecmw_config.h"
#include "hecmw_etype.h"
#include "hecmw_result_io.h"

struct hecmwST_result_io_data ResIO;

static int is_valid_label(char *label) {
#define ALLOW_CHAR_FIRST "_" /* and alphabet */
#define ALLOW_CHAR "_-+"     /* and alphabet, digit */
  int c;
  char *p, *q;

  if (label == NULL) return 0;

  c = label[0];

  /* check first character */
  if (!isalpha(c)) {
    q = ALLOW_CHAR_FIRST;
    while (*q) {
      if (*q == c) break;
      q++;
    }
    if (!*q) return 0;
  }

  /* check 2nd character or later */
  p = &label[1];
  while (*p) {
    if (!isalnum(*p)) {
      q = ALLOW_CHAR;
      while (*q) {
        if (*q == *p) break;
        q++;
      }
      if (!*q) return 0;
    }
    p++;
  }
  return 1;
}

void HECMW_result_io_finalize() {
  struct result_list *p, *q;

  for (p = ResIO.global_list; p; p = q) {
    q = p->next;
    HECMW_free(p->label);
    HECMW_free(p->ptr);
    HECMW_free(p);
  }
  ResIO.global_list = NULL;

  for (p = ResIO.node_list; p; p = q) {
    q = p->next;
    HECMW_free(p->label);
    HECMW_free(p->ptr);
    HECMW_free(p);
  }
  ResIO.node_list = NULL;

  for (p = ResIO.elem_list; p; p = q) {
    q = p->next;
    HECMW_free(p->label);
    HECMW_free(p->ptr);
    HECMW_free(p);
  }
  ResIO.elem_list = NULL;

  ResIO.nnode = ResIO.nelem = 0;
  strcpy(ResIO.head, "");

  if (ResIO.MPC_exist) {
    ResIO.MPC_exist = 0;
    HECMW_free(ResIO.eid_wo_MPC);
    HECMW_free(ResIO.elem_global_ID);
  }

  ResIO.node_global_ID = NULL;
  ResIO.elem_global_ID = NULL;
}

static int setup_MPC(int n_elem_type, int *elem_type_index, int *elem_type_item,
                     int *elemID) {
  int itype, ic_type, is, ie, icel;
  int nelem_wo_MPC;
  int *elemID_wo_MPC;

  ResIO.MPC_exist = 0;
  ResIO.eid_wo_MPC = NULL;

  for (itype = 0; itype < n_elem_type; itype++) {
    ic_type = elem_type_item[itype];
    if (HECMW_is_etype_link(ic_type) ||
        HECMW_is_etype_patch(ic_type) ||
        HECMW_is_etype_smoothing(ic_type)) {
      ResIO.MPC_exist = 1;
      break;
    }
  }

  if (ResIO.MPC_exist) {
    ResIO.eid_wo_MPC = (int *) HECMW_calloc(ResIO.nelem, sizeof(int));
    if (ResIO.eid_wo_MPC == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    }
    elemID_wo_MPC = (int *) HECMW_calloc(ResIO.nelem, sizeof(int));
    if (elemID_wo_MPC == NULL) {
      HECMW_set_error(errno, "");
      goto error;
    }

    nelem_wo_MPC = 0;
    for (itype = 0; itype < n_elem_type; itype++) {
      ic_type = elem_type_item[itype];
      if (HECMW_is_etype_link(ic_type) ||
          HECMW_is_etype_patch(ic_type) ||
          HECMW_is_etype_smoothing(ic_type)) {
        continue;
      }
      is = elem_type_index[itype];
      ie = elem_type_index[itype + 1];
      for (icel = is; icel < ie; icel++) {
        if (icel >= ResIO.nelem) {
          /*
           * Safeguard for inconsistent n_elem and elem_type_index
           */
          HECMW_log(HECMW_LOG_WARN, "result output: ignoring elements type=%d, %d..%d (n_elem=%d)\n",
                  ic_type, icel+1, ie, ResIO.nelem);
          break;
        }
        ResIO.eid_wo_MPC[nelem_wo_MPC] = icel;
        elemID_wo_MPC[nelem_wo_MPC] = elemID[icel];
        nelem_wo_MPC++;
      }
    }
    ResIO.nelem = nelem_wo_MPC;
    ResIO.elem_global_ID = elemID_wo_MPC;
  }

  return 0;
error:
  return -1;
}

int HECMW_result_io_init(int n_node, int n_elem, int *nodeID, int *elemID,
                         int n_elem_type, int *elem_type_index, int *elem_type_item,
                         int i_step, char *header, char *comment) {
  size_t len;
  int rtc;
  char *p, *q;

  ResIO.nnode = n_node;
  ResIO.nelem = n_elem;
  ResIO.istep = i_step;

  ResIO.node_global_ID = nodeID;
  ResIO.elem_global_ID = elemID;

  if (header == NULL) {
    ResIO.head[0] = '\0';
    return 0;
  }

  len = 0;
  p   = header;
  q   = ResIO.head;
  while (len < sizeof(ResIO.head) - 1 && *p && *p != '\n') {
    *q++ = *p++;
    len++;
  }
  *q++ = '\0';

  if (comment == NULL) {
    ResIO.comment_line[0] = '\0';
    return 0;
  }

  len = 0;
  p   = comment;
  q   = ResIO.comment_line;
  while (len < sizeof(ResIO.comment_line) - 1 && *p && *p != '\n') {
    *q++ = *p++;
    len++;
  }
  *q++ = '\0';

  rtc = setup_MPC(n_elem_type, elem_type_index, elem_type_item, elemID);
  if (rtc != 0) {
    goto error;
  }

  return 0;
error:
  return -1;
}

static int add_to_global_list(struct result_list *result) {
  struct result_list *p, *q;

  q = NULL;
  for (p = ResIO.global_list; p; p = (q = p)->next)
    ;

  if (q == NULL) {
    ResIO.global_list = result;
  } else {
    q->next = result;
  }
  return 0;
}

static int add_to_node_list(struct result_list *result) {
  struct result_list *p, *q;

  q = NULL;
  for (p = ResIO.node_list; p; p = (q = p)->next)
    ;

  if (q == NULL) {
    ResIO.node_list = result;
  } else {
    q->next = result;
  }
  return 0;
}

static int add_to_elem_list(struct result_list *result) {
  struct result_list *p, *q;

  q = NULL;
  for (p = ResIO.elem_list; p; p = (q = p)->next)
    ;

  if (q == NULL) {
    ResIO.elem_list = result;
  } else {
    q->next = result;
  }
  return 0;
}

static struct result_list *make_result_list(int n_dof, char *label,
                                            double *ptr) {
  struct result_list *result;
  char *new_label;

  result = HECMW_malloc(sizeof(*result));
  if (result == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  new_label = HECMW_strdup(label);
  if (new_label == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  result->label = new_label;
  result->ptr   = ptr;
  result->n_dof = n_dof;
  result->next  = NULL;

  return result;
error:
  HECMW_free(result);
  HECMW_free(new_label);
  return NULL;
}

int HECMW_result_io_add(int dtype, int n_dof, char *label, double *ptr) {
  struct result_list *result;
  size_t n, size;
  double *data;
  int i, icel, idof;

  if (dtype < HECMW_RESULT_DTYPE_MIN && dtype > HECMW_RESULT_DTYPE_MAX) {
    HECMW_set_error(HECMW_UTIL_E0206, "");
    goto error;
  }

  if (!is_valid_label(label)) {
    HECMW_set_error(HECMW_UTIL_E0207, "");
    goto error;
  }

  if (dtype == HECMW_RESULT_DTYPE_NODE) {
    n = ResIO.nnode;
  } else if (dtype == HECMW_RESULT_DTYPE_ELEM) {
    n = ResIO.nelem;
  } else { // dtype == HECMW_RESULT_DTYPE_GLOBAL
    n = 1;
  }
  size = sizeof(double) * n * n_dof;
  data = (double *) HECMW_calloc(n * n_dof, sizeof(double));
  if (data == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  if (dtype == HECMW_RESULT_DTYPE_ELEM && ResIO.MPC_exist) {
    for (i = 0; i < ResIO.nelem; i++) {
      icel = ResIO.eid_wo_MPC[i];
      for (idof = 0; idof < n_dof; idof++) {
        data[n_dof * i + idof] = ptr[n_dof * icel + idof];
      }
    }
  } else {
    memcpy(data, ptr, size);
  }

  result = make_result_list(n_dof, label, data);
  if (result == NULL) {
    goto error;
  }

  if (dtype == HECMW_RESULT_DTYPE_NODE) {
    if (add_to_node_list(result)) goto error;
  } else if (dtype == HECMW_RESULT_DTYPE_ELEM)  {
    if (add_to_elem_list(result)) goto error;
  } else { /* dtype == HECMW_RESULT_DTYPE_GLOBAL */
    if (add_to_global_list(result)) goto error;
  }

  return 0;
error:
  return -1;
}

int HECMW_result_io_count_ng_comp(void) {
  int ng_comp;
  struct result_list *p;

  ng_comp = 0;
  for (p = ResIO.global_list; p; p = p->next) {
    ng_comp++;
  }
  return ng_comp;
}

int HECMW_result_io_count_nn_comp(void) {
  int nn_comp;
  struct result_list *p;

  nn_comp = 0;
  for (p = ResIO.node_list; p; p = p->next) {
    nn_comp++;
  }
  return nn_comp;
}

int HECMW_result_io_count_ne_comp(void) {
  int ne_comp;
  struct result_list *p;

  ne_comp = 0;
  for (p = ResIO.elem_list; p; p = p->next) {
    ne_comp++;
  }
  return ne_comp;
}
