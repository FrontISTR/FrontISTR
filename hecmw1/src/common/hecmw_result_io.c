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

int IStep;
int NNode;
int NElem;
#ifdef OLD_RES_FORMAT
int filever_major=1;
int filever_minor=0;
#else
int FileVer_Major=2;
int FileVer_Minor=0;
#endif // OLD_RES_FORMAT
char Head[HECMW_HEADER_LEN + 1];
char Comment_Line[HECMW_MSG_LEN + 1];
char Line_Buf[LINEBUF_SIZE + 1];

struct result_list *Global_List;
struct result_list *Node_List;
struct result_list *Elem_List;

int *Node_Global_ID = NULL;
int *Elem_Global_ID = NULL;

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

  for (p = Global_List; p; p = q) {
    q = p->next;
    HECMW_free(p->label);
    HECMW_free(p);
  }
  Global_List = NULL;

  for (p = Node_List; p; p = q) {
    q = p->next;
    HECMW_free(p->label);
    HECMW_free(p);
  }
  Node_List = NULL;

  for (p = Elem_List; p; p = q) {
    q = p->next;
    HECMW_free(p->label);
    HECMW_free(p);
  }
  Elem_List = NULL;

  NNode = NElem = 0;
  strcpy(Head, "");

  Node_Global_ID = NULL;
  Elem_Global_ID = NULL;
}

int HECMW_result_io_init(int n_node, int n_elem, int *nodeID, int *elemID,
                           int i_step, char *header, char *comment) {
  int len;
  char *p, *q;

  NNode = n_node;
  NElem = n_elem;
  IStep = i_step;

  Node_Global_ID = nodeID;
  Elem_Global_ID = elemID;

  if (header == NULL) {
    Head[0] = '\0';
    return 0;
  }

  len = 0;
  p   = header;
  q   = Head;
  while (len < sizeof(Head) - 1 && *p && *p != '\n') {
    *q++ = *p++;
    len++;
  }
  *q++ = '\0';

  if (comment == NULL) {
    Comment_Line[0] = '\0';
    return 0;
  }

  len = 0;
  p   = comment;
  q   = Comment_Line;
  while (len < sizeof(Comment_Line) - 1 && *p && *p != '\n') {
    *q++ = *p++;
    len++;
  }
  *q++ = '\0';

  return 0;
}

static int add_to_global_list(struct result_list *result) {
  struct result_list *p, *q;

  q = NULL;
  for (p = Global_List; p; p = (q = p)->next)
    ;

  if (q == NULL) {
    Global_List = result;
  } else {
    q->next = result;
  }
  return 0;
}

static int add_to_node_list(struct result_list *result) {
  struct result_list *p, *q;

  q = NULL;
  for (p = Node_List; p; p = (q = p)->next)
    ;

  if (q == NULL) {
    Node_List = result;
  } else {
    q->next = result;
  }
  return 0;
}

static int add_to_elem_list(struct result_list *result) {
  struct result_list *p, *q;

  q = NULL;
  for (p = Elem_List; p; p = (q = p)->next)
    ;

  if (q == NULL) {
    Elem_List = result;
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

  if (dtype < HECMW_RESULT_DTYPE_MIN && dtype > HECMW_RESULT_DTYPE_MAX) {
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

  if (dtype == HECMW_RESULT_DTYPE_NODE) {
    if (add_to_node_list(result)) goto error;
  } else if (dtype == HECMW_RESULT_DTYPE_ELEM)  {
    if (add_to_elem_list(result)) goto error;
  } else { /* dtype == HECMW_RESULT_DTYPE_GLOBAL */
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
  for (p = Global_List; p; p = p->next) {
    ng_comp++;
  }
  return ng_comp;
}

int HECMW_result_io_count_nn_comp(void) {
  int nn_comp;
  struct result_list *p;

  nn_comp = 0;
  for (p = Node_List; p; p = p->next) {
    nn_comp++;
  }
  return nn_comp;
}

int HECMW_result_io_count_ne_comp(void) {
  int ne_comp;
  struct result_list *p;

  ne_comp = 0;
  for (p = Elem_List; p; p = p->next) {
    ne_comp++;
  }
  return ne_comp;
}
