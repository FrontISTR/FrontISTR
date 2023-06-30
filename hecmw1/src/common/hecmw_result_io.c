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
#include "hecmw_result_io.h"

int istep;
int nnode;
int nelem;
#ifdef OLD_RES_FORMAT
int filever_major=1;
int filever_minor=0;
#else
int filever_major=2;
int filever_minor=0;
#endif // OLD_RES_FORMAT
char head[HECMW_HEADER_LEN + 1];
char comment_line[HECMW_MSG_LEN + 1];
char line_buf[LINEBUF_SIZE + 1];

struct result_list *global_list;
struct result_list *node_list;
struct result_list *elem_list;

int *node_global_ID = NULL;
int *elem_global_ID = NULL;

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

  for (p = global_list; p; p = q) {
    q = p->next;
    HECMW_free(p->label);
    HECMW_free(p);
  }
  global_list = NULL;

  for (p = node_list; p; p = q) {
    q = p->next;
    HECMW_free(p->label);
    HECMW_free(p);
  }
  node_list = NULL;

  for (p = elem_list; p; p = q) {
    q = p->next;
    HECMW_free(p->label);
    HECMW_free(p);
  }
  elem_list = NULL;

  nnode = nelem = 0;
  strcpy(head, "");

  node_global_ID = NULL;
  elem_global_ID = NULL;
}

int HECMW_result_io_init(int n_node, int n_elem, int *nodeID, int *elemID,
                           int i_step, char *header, char *comment) {
  int len;
  char *p, *q;

  nnode = n_node;
  nelem = n_elem;
  istep = i_step;

  node_global_ID = nodeID;
  elem_global_ID = elemID;

  if (header == NULL) {
    head[0] = '\0';
    return 0;
  }

  len = 0;
  p   = header;
  q   = head;
  while (len < sizeof(head) - 1 && *p && *p != '\n') {
    *q++ = *p++;
    len++;
  }
  *q++ = '\0';

  if (comment == NULL) {
    comment_line[0] = '\0';
    return 0;
  }

  len = 0;
  p   = comment;
  q   = comment_line;
  while (len < sizeof(comment_line) - 1 && *p && *p != '\n') {
    *q++ = *p++;
    len++;
  }
  *q++ = '\0';

  return 0;
}

static int add_to_global_list(struct result_list *result) {
  struct result_list *p, *q;

  q = NULL;
  for (p = global_list; p; p = (q = p)->next)
    ;

  if (q == NULL) {
    global_list = result;
  } else {
    q->next = result;
  }
  return 0;
}

static int add_to_node_list(struct result_list *result) {
  struct result_list *p, *q;

  q = NULL;
  for (p = node_list; p; p = (q = p)->next)
    ;

  if (q == NULL) {
    node_list = result;
  } else {
    q->next = result;
  }
  return 0;
}

static int add_to_elem_list(struct result_list *result) {
  struct result_list *p, *q;

  q = NULL;
  for (p = elem_list; p; p = (q = p)->next)
    ;

  if (q == NULL) {
    elem_list = result;
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

  if (dtype < 1 && dtype > 3) {
    HECMW_set_error(HECMW_UTIL_E0206, "");
    goto error;
  }

  if (!is_valid_label(label)) {
    HECMW_set_error(HECMW_UTIL_E0207, "");
    goto error;
  }

  result = make_result_list(n_dof, label, ptr);
  if (result == NULL) {
    goto error;
  }

  if (dtype == 1) {
    /* node */
    if (add_to_node_list(result)) goto error;
  } else if (dtype == 2)  {
    /* elem */
    if (add_to_elem_list(result)) goto error;
  } else {
    /* global */
    if (add_to_global_list(result)) goto error;
  }

  return 0;
error:
  HECMW_free(result);
  return -1;
}

int HECMW_result_io_count_ng_comp(void) {
  int ng_comp;
  struct result_list *p;

  ng_comp = 0;
  for (p = global_list; p; p = p->next) {
    ng_comp++;
  }
  return ng_comp;
}

int HECMW_result_io_count_nn_comp(void) {
  int nn_comp;
  struct result_list *p;

  nn_comp = 0;
  for (p = node_list; p; p = p->next) {
    nn_comp++;
  }
  return nn_comp;
}

int HECMW_result_io_count_ne_comp(void) {
  int ne_comp;
  struct result_list *p;

  ne_comp = 0;
  for (p = elem_list; p; p = p->next) {
    ne_comp++;
  }
  return ne_comp;
}
