/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "hecmw_util.h"
#include "hecmw_result.h"
#include "hecmw_result_io.h"

#define COL_INT 10
#define COL_DOUBLE 5

#define LINEBUF_SIZE 1023
static char Line_Buf[LINEBUF_SIZE + 1];

/*---------------------------------------------------------------------------*/
/* TEXT MODE I/O ---- output_result                                          */
/*---------------------------------------------------------------------------*/


static int output_result_header(FILE *fp) {
  int rc;

  /* header */
  if( HECMW_RESULT_FILEVER_MAJOR > 1 ){
    sprintf(ResIO.head,"%s %d.%d",ResIO.head,HECMW_RESULT_FILEVER_MAJOR,HECMW_RESULT_FILEVER_MINOR);
  }
  rc = fprintf(fp, "%s\n", ResIO.head);
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "head");
    return -1;
  }

  return 0;
}


static int output_result_global(FILE *fp) {
  int i,j,k,n,rc,ng_comp;
  struct result_list *p,**data;

  /* comment */
  rc = fprintf(fp, "*comment\n");
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "*comment");
    return -1;
  }
  rc = fprintf(fp, "%s\n", ResIO.comment_line);
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "comment");
    return -1;
  }

  /* global header */
  rc = fprintf(fp, "*global\n");
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "*global");
    return -1;
  }

  /* ng_component */
  rc = fprintf(fp, "%d\n", HECMW_result_io_count_ng_comp());
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "ng_comp");
    return -1;
  }

  /* ng_dof */
  n = 0;
  for(p=ResIO.global_list; p; p=p->next) {
    rc = fprintf(fp, "%d%c", p->n_dof, (n+1)%COL_INT ? ' ' : '\n');
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "ng_dof");
      return -1;
    }
    n++;
  }
  if(n % COL_INT) {
    rc = fprintf(fp, "\n");
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "");
      return -1;
    }
  }

  /* global_label */
  for(p=ResIO.global_list; p; p=p->next) {
    rc = fprintf(fp, "%s\n", p->label);
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "global_label");
      return -1;
    }
  }

  /* global_val_item */
  ng_comp = HECMW_result_io_count_ng_comp();
  if(ng_comp == 0) return 0;
  data = HECMW_malloc(sizeof(*data) * ng_comp);
  if(data == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }
  i = 0;
  for(p=ResIO.global_list; p; p=p->next) {
    data[i++] = p;
  }
  n = 0;
  for(j=0; j < ng_comp; j++) {
    p = data[j];
    for(k=0; k < p->n_dof; k++) {
      rc = fprintf(fp, "%.16E%c", p->ptr[k], (n+1)%COL_DOUBLE ? ' ' : '\n');
      if(rc < 0) {
        HECMW_set_error(HECMW_UTIL_E0205, "global_val_item");
        return -1;
      }
      n++;
    }
  }
  if(n % COL_DOUBLE) {
    rc = fprintf(fp, "\n");
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "");
      return -1;
    }
  }
  HECMW_free(data);

  return 0;
}


static int output_result_dataheader(FILE *fp) {
  int rc;

  /* n_node, n_elem */
  rc = fprintf(fp, "%d %d\n", ResIO.nnode, ResIO.nelem);
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "nnode,nelem");
    return -1;
  }

  /* nn_component, ne_component */
  rc = fprintf(fp, "%d %d\n", HECMW_result_io_count_nn_comp(), HECMW_result_io_count_ne_comp());
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "nn_comp,ne_comp");
    return -1;
  }

  return 0;
}


static int output_result_node(FILE *fp) {
  int i,j,k,n,rc,nn_comp;
  struct result_list *p,**data;

  /* nn_dof */
  n = 0;
  for(p=ResIO.node_list; p; p=p->next) {
    rc = fprintf(fp, "%d%c", p->n_dof, (n+1)%COL_INT ? ' ' : '\n');
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "nn_dof");
      return -1;
    }
    n++;
  }
  if(n % COL_INT) {
    rc = fprintf(fp, "\n");
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "");
      return -1;
    }
  }

  /* node_label */
  for(p=ResIO.node_list; p; p=p->next) {
    rc = fprintf(fp, "%s\n", p->label);
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "node_label");
      return -1;
    }
  }

  /* node_val_item */
  nn_comp = HECMW_result_io_count_nn_comp();
  if(nn_comp == 0) return 0;
  data = HECMW_malloc(sizeof(*data) * nn_comp);
  if(data == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }
  i = 0;
  for(p=ResIO.node_list; p; p=p->next) {
    data[i++] = p;
  }
  for(i=0; i < ResIO.nnode; i++) {
    rc = fprintf(fp, "%d \n", ResIO.node_global_ID[i]);
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "node_global_ID");
      return -1;
    }
    n = 0;
    for(j=0; j < nn_comp; j++) {
      p = data[j];
      for(k=0; k < p->n_dof; k++) {
        rc = fprintf(fp, "%.16E%c", p->ptr[i*p->n_dof+k], (n+1)%COL_DOUBLE ? ' ' : '\n');
        if(rc < 0) {
          HECMW_set_error(HECMW_UTIL_E0205, "node_val_item");
          return -1;
        }
        n++;
      }
    }
    if(n % COL_DOUBLE) {
      rc = fprintf(fp, "\n");
      if(rc < 0) {
        HECMW_set_error(HECMW_UTIL_E0205, "");
        return -1;
      }
    }
  }
  HECMW_free(data);

  return 0;
}


static int output_result_elem(FILE *fp) {
  int i,j,k,n,rc,ne_comp;
  struct result_list *p,**data;

  /* ne_dof */
  n = 0;
  for(p=ResIO.elem_list; p; p=p->next) {
    rc = fprintf(fp, "%d%c", p->n_dof, (n+1)%COL_INT ? ' ' : '\n');
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "ne_dof");
      return -1;
    }
    n++;
  }
  if(n % COL_INT) {
    rc = fprintf(fp, "\n");
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "");
      return -1;
    }
  }

  /* elem_label */
  for(p=ResIO.elem_list; p; p=p->next) {
    rc = fprintf(fp, "%s\n", p->label);
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "elem_label");
      return -1;
    }
  }

  /* elem_val_item */
  ne_comp = HECMW_result_io_count_ne_comp();
  if(ne_comp == 0) return 0;
  data = HECMW_malloc(sizeof(*data) * ne_comp);
  if(data == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }
  i = 0;
  for(p=ResIO.elem_list; p; p=p->next) {
    data[i++] = p;
  }
  for(i=0; i < ResIO.nelem; i++) {
    rc = fprintf(fp, "%d\n", ResIO.elem_global_ID[i]);
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "elem_global_ID");
      return -1;
    }
    n = 0;
    for(j=0; j < ne_comp; j++) {
      p = data[j];
      for(k=0; k < p->n_dof; k++) {
        rc = fprintf(fp, "%.16E%c", p->ptr[i*p->n_dof+k], (n+1)%COL_DOUBLE ? ' ' : '\n');
        if(rc < 0) {
          HECMW_set_error(HECMW_UTIL_E0205, "elem_val_item");
          return -1;
        }
        n++;
      }
    }
    if(n % COL_DOUBLE) {
      rc = fprintf(fp, "\n");
      if(rc < 0) {
        HECMW_set_error(HECMW_UTIL_E0205, "");
        return -1;
      }
    }
  }
  HECMW_free(data);

  return 0;
}


static int output_result_data(FILE *fp) {
  int rc;
  HECMW_assert(fp);

  if(output_result_header(fp)) {
    return -1;
  }
  if( HECMW_RESULT_FILEVER_MAJOR > 1 ){
    if(output_result_global(fp)) {
      return -1;
    }
    /* data header */
    rc = fprintf(fp, "*data\n");
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "*data");
      return -1;
    }
  }
  if(output_result_dataheader(fp)) {
    return -1;
  }
  if(output_result_node(fp)) {
    return -1;
  }
  if(output_result_elem(fp)) {
    return -1;
  }

  return 0;
}

/*---------------------------------------------------------------------------*/

int HECMW_result_io_txt_write_by_fname(char *filename) {
  FILE *fp = NULL;

  if (HECMW_ctrl_is_subdir()) {
    if (HECMW_ctrl_make_subdir(filename)) {
      HECMW_set_error(HECMW_UTIL_E0201, "File: %s, %s", filename,
                      HECMW_strmsg(errno));
      goto error;
    }
  }

  if ((fp = fopen(filename, "w")) == NULL) {
    HECMW_set_error(HECMW_UTIL_E0201, "File: %s, %s", filename,
                    HECMW_strmsg(errno));
    goto error;
  }

  if (output_result_data(fp)) {
    goto error;
  }

  if (fclose(fp)) {
    HECMW_set_error(HECMW_UTIL_E0202, HECMW_strmsg(errno));
    goto error;
  }
  fp = NULL;

  return 0;
error:
  if (fp) fclose(fp);
  return -1;
}


/*---------------------------------------------------------------------------*/
/* TEXT MODE I/O ---- output_result_ST                                       */
/*---------------------------------------------------------------------------*/


static int output_result_header_ST(struct hecmwST_result_data *result, char *header, FILE *fp) {
  size_t len;
  int rc;
  char *p,*q;
  char head[HECMW_HEADER_LEN+1];

  if(header == NULL) {
    head[0] = '\0';
  } else {
    len = 0;
    p = header;
    q = head;
    while(len < sizeof(head)-1 && *p && *p != '\n') {
      *q++ = *p++;
      len++;
    }
    *q++ = '\0';
  }

  /* header */
  if( HECMW_RESULT_FILEVER_MAJOR > 1 ){
    sprintf(head,"%s %d.%d",head,HECMW_RESULT_FILEVER_MAJOR,HECMW_RESULT_FILEVER_MINOR);
  }
  rc = fprintf(fp, "%s\n", head);
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "header");
    return -1;
  }

  return 0;
}


static int output_result_global_ST(struct hecmwST_result_data *result, char *comment, FILE *fp) {
  size_t len;
  int i,j,k,n,rc;
  char *p,*q;
  char comment_line[HECMW_MSG_LEN+1];

  if(comment == NULL) {
    comment_line[0] = '\0';
  } else {
    len = 0;
    p = comment;
    q = comment_line;
    while(len < sizeof(comment_line)-1 && *p && *p != '\n') {
      *q++ = *p++;
      len++;
    }
    *q++ = '\0';
  }

  /* comment */
  rc = fprintf(fp, "*comment\n");
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "*comment");
    return -1;
  }
  rc = fprintf(fp, "%s\n", comment);
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "comment");
    return -1;
  }


  /* global header */
  rc = fprintf(fp, "*global\n");
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "*global");
    return -1;
  }

  /* ng_component */
  rc = fprintf(fp, "%d\n", result->ng_component);
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "ng_comp");
    return -1;
  }

  /* ng_dof */
  n = 0;
  for(i=0; i < result->ng_component; i++) {
    rc = fprintf(fp, "%d%c", result->ng_dof[i], (n+1)%COL_INT ? ' ' : '\n');
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "ng_dof");
      return -1;
    }
    n++;
  }
  if(n % COL_INT) {
    rc = fprintf(fp, "\n");
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "global_label");
      return -1;
    }
  }

  /* global_label */
  for(i=0; i < result->ng_component; i++) {
    rc = fprintf(fp, "%s\n", result->global_label[i]);
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "");
      return -1;
    }
  }

  /* global_val_item */
  if(result->ng_component == 0) return 0;
  n = 0;
  for(j=0; j < result->ng_component; j++) {
    for(k=0; k < result->ng_dof[j]; k++) {
      rc = fprintf(fp, "%.16E%c", result->global_val_item[n], (n+1)%COL_DOUBLE ? ' ' : '\n');
      if(rc < 0) {
        HECMW_set_error(HECMW_UTIL_E0205, "global_val_item");
        return -1;
      }
      n++;
    }
  }
  if(n % COL_DOUBLE) {
    rc = fprintf(fp, "\n");
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "");
      return -1;
    }
  }

  /* dataheader */
  rc = fprintf(fp, "*data\n");
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "*data");
    return -1;
  }

  return 0;
}


static int output_result_dataheader_ST(struct hecmwST_result_data *result,
                                       int n_node, int n_elem, FILE *fp) {
  int rc;

  /* n_node, n_elem */
  rc = fprintf(fp, "%d %d\n", n_node, n_elem);
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "n_node,n_elem");
    return -1;
  }

  /* nn_component, ne_component */
  rc = fprintf(fp, "%d %d\n", result->nn_component, result->ne_component);
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "nn_comp,ne_comp");
    return -1;
  }

  return 0;
}


static int output_result_node_ST(struct hecmwST_result_data *result, int n_node, FILE *fp) {
  int i,j,k,n,m,rc;

  /* nn_dof */
  n = 0;
  for(i=0; i < result->nn_component; i++) {
    rc = fprintf(fp, "%d%c", result->nn_dof[i], (n+1)%COL_INT ? ' ' : '\n');
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "nn_dof");
      return -1;
    }
    n++;
  }
  if(n % COL_INT) {
    rc = fprintf(fp, "\n");
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "node_label");
      return -1;
    }
  }

  /* node_label */
  for(i=0; i < result->nn_component; i++) {
    rc = fprintf(fp, "%s\n", result->node_label[i]);
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "");
      return -1;
    }
  }

  /* node_val_item */
  if(result->nn_component == 0) return 0;
  m = 0;
  for(i=0; i < n_node; i++) {
    rc = fprintf(fp, "%d \n", ResIO.node_global_ID[i]);
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "node_global_ID");
      return -1;
    }
    n = 0;
    for(j=0; j < result->nn_component; j++) {
      for(k=0; k < result->nn_dof[j]; k++) {
        rc = fprintf(fp, "%.16E%c", result->node_val_item[m], (n+1)%COL_DOUBLE ? ' ' : '\n');
        if(rc < 0) {
          HECMW_set_error(HECMW_UTIL_E0205, "node_val_item");
          return -1;
        }
        n++;
        m++;
      }
    }
    if(n % COL_DOUBLE) {
      rc = fprintf(fp, "\n");
      if(rc < 0) {
        HECMW_set_error(HECMW_UTIL_E0205, "");
        return -1;
      }
    }
  }

  return 0;
}


static int output_result_elem_ST(struct hecmwST_result_data *result, int n_elem, FILE *fp) {
  int i,j,k,n,m,rc;

  /* ne_dof */
  n = 0;
  for(i=0; i < result->ne_component; i++) {
    rc = fprintf(fp, "%d%c", result->ne_dof[i], (n+1)%COL_INT ? ' ' : '\n');
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "ne_dof");
      return -1;
    }
    n++;
  }
  if(n % COL_INT) {
    rc = fprintf(fp, "\n");
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "");
      return -1;
    }
  }

  /* elem_label */
  for(i=0; i < result->ne_component; i++) {
    rc = fprintf(fp, "%s\n", result->elem_label[i]);
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "elem_label");
      return -1;
    }
  }

  /* elem_val_item */
  if(result->ne_component == 0) return 0;
  m = 0;
  for(i=0; i < n_elem; i++) {
    rc = fprintf(fp, "%d\n", ResIO.elem_global_ID[i]);
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "elem_global_ID");
      return -1;
    }
    n = 0;
    for(j=0; j < result->ne_component; j++) {
      for(k=0; k < result->ne_dof[j]; k++) {
        rc = fprintf(fp, "%.16E%c", result->elem_val_item[m], (n+1)%COL_DOUBLE ? ' ' : '\n');
        if(rc < 0) {
          HECMW_set_error(HECMW_UTIL_E0205, "elem_val_item");
          return -1;
        }
        n++;
        m++;
      }
    }
    if(n % COL_DOUBLE) {
      rc = fprintf(fp, "\n");
      if(rc < 0) {
        HECMW_set_error(HECMW_UTIL_E0205, "");
        return -1;
      }
    }
  }

  return 0;
}


static int output_result_data_ST(struct hecmwST_result_data *result, int n_node, int n_elem,
                                 char *header, char *comment, FILE *fp) {
  HECMW_assert(fp);

  if(output_result_header_ST(result, header, fp)) {
    return -1;
  }
  if( HECMW_RESULT_FILEVER_MAJOR > 1 ){
    if(output_result_global_ST(result, comment, fp)) {
      return -1;
    }
  }
  if(output_result_dataheader_ST(result, n_node, n_elem, fp)) {
    return -1;
  }
  if(output_result_node_ST(result, n_node, fp)) {
    return -1;
  }
  if(output_result_elem_ST(result, n_elem, fp)) {
    return -1;
  }

  return 0;
}

/*---------------------------------------------------------------------------*/

int HECMW_result_io_txt_write_ST_by_fname(char *filename,
                                       struct hecmwST_result_data *result,
                                       int n_node, int n_elem, char *header, char *comment) {
  FILE *fp = NULL;

  if (HECMW_ctrl_is_subdir()) {
    if (HECMW_ctrl_make_subdir(filename)) {
      HECMW_set_error(HECMW_UTIL_E0201, "File: %s, %s", filename,
                      HECMW_strmsg(errno));
      goto error;
    }
  }

  if ((fp = fopen(filename, "w")) == NULL) {
    HECMW_set_error(HECMW_UTIL_E0201, "File: %s, %s", filename,
                    HECMW_strmsg(errno));
    goto error;
  }

  if (output_result_data_ST(result, n_node, n_elem, header, comment, fp)) {
    goto error;
  }

  if (fclose(fp)) {
    HECMW_set_error(HECMW_UTIL_E0202, HECMW_strmsg(errno));
    goto error;
  }
  fp = NULL;

  return 0;
error:
  if (fp) fclose(fp);
  return -1;
}


/*---------------------------------------------------------------------------*/
/* TEXT MODE I/O ---- input_result                                           */
/*---------------------------------------------------------------------------*/


static int get_line(char *buf, int bufsize, FILE *fp) {
  if(fgets(buf, bufsize, fp) == NULL) {
    HECMW_set_error(HECMW_UTIL_E0205, "get_line");
    return -1;
  }
  return strlen(Line_Buf);
}


static int input_result_header(struct hecmwST_result_data *result, FILE *fp) {
  char *ptr;

  /* header */
  if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
    return -1;
  }
  Line_Buf[ strlen(Line_Buf)-1 ] = 0;/* remove CR/LF*/
  if( HECMW_RESULT_FILEVER_MAJOR > 1 ){
    ptr = strtok(Line_Buf, " ");
    sprintf(Line_Buf, "%s", ptr);
  }
  strcpy( ResIO.head, Line_Buf );

  return 0;
}


static int input_result_global(struct hecmwST_result_data *result, FILE *fp) {
#define DELIM " \n"
  int i,rc,n;
  char *p,*buf;

  /* comment */
  if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) { //skip comment header
    return -1;
  }
  if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
    return -1;
  }
  Line_Buf[ strlen(Line_Buf)-1 ] = 0;/* remove CR/LF*/
  strcpy( ResIO.comment_line, Line_Buf );


  /* skip global header */
  if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
    return -1;
  }

  /* ng_component */
  if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
    return -1;
  }
  if(sscanf(Line_Buf, "%d", &result->ng_component) != 1) {
    HECMW_set_error(HECMW_UTIL_E0205, "ng_comp");
    return -1;
  }

  if(result->ng_component <= 0) {
    return 0;
  }

  /* ng_dof */
  result->ng_dof = HECMW_malloc(sizeof(*result->ng_dof)*result->ng_component);
  if(result->ng_dof == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
    return -1;
  }
  i = n = 0;
  buf = Line_Buf;
  while(i < result->ng_component) {
    p = strtok(buf, DELIM);
    if(p == NULL) {
      if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
        return -1;
      }
      buf = Line_Buf;
      continue;
    }
    buf = NULL;
    rc = sscanf(p, "%d", &result->ng_dof[i]);
    if(rc == EOF) {
      HECMW_set_error(HECMW_UTIL_E0204, "");
      return -1;
    }
    if(rc != 1) {
      HECMW_set_error(HECMW_UTIL_E0205, "ng_dof");
      return -1;
    }
    n += result->ng_dof[i];
    i++;
  }

  /* global_label */
  result->global_label = HECMW_malloc(sizeof(*result->global_label)*result->ng_component);
  if(result->global_label == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  for(i=0; i < result->ng_component; i++) {
    char label[HECMW_NAME_LEN+1];
    if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
      return -1;
    }
    rc = sscanf(Line_Buf, "%s", label);
    if(rc == EOF) {
      HECMW_set_error(HECMW_UTIL_E0204, "");
      return -1;
    }
    if(rc != 1) {
      HECMW_set_error(HECMW_UTIL_E0205, "global_label");
      return -1;
    }
    result->global_label[i] = HECMW_strdup(label);
    if(result->global_label[i] == NULL) {
      HECMW_set_error(errno, "");
      return -1;
    }
  }
  /* global_val_item */
  result->global_val_item = HECMW_malloc(sizeof(*result->global_val_item)*n);
  if(result->global_val_item == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
    return -1;
  }
  i = 0;
  buf = Line_Buf;
  while(i < n) {
    p = strtok(buf, DELIM);
    if(p == NULL) {
      if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
        return -1;
      }
      buf = Line_Buf;
      continue;
    }
    buf = NULL;
    rc = sscanf(p, "%lf", &result->global_val_item[i]);
    if(rc == EOF) {
      HECMW_set_error(HECMW_UTIL_E0204, "");
      return -1;
    }
    if(rc != 1) {
      HECMW_set_error(HECMW_UTIL_E0205, "global_val_item");
      return -1;
    }
    i++;
  }

  /* skip data header */
  if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
    return -1;
  }

  return 0;
}


static int input_result_dataheader(struct hecmwST_result_data *result,
                                   int *n_node, int *n_elem, FILE *fp) {

  /* n_node, n_elem */
  if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
    return -1;
  }
  if(sscanf(Line_Buf, "%d%d", n_node, n_elem) != 2) {
    HECMW_set_error(HECMW_UTIL_E0205, "n_node,n_elem");
    return -1;
  }

  /* nn_component, ne_component */
  if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
    return -1;
  }
  if(sscanf(Line_Buf, "%d%d", &result->nn_component, &result->ne_component) != 2) {
    HECMW_set_error(HECMW_UTIL_E0205, "nn_comp,ne_comp");
    return -1;
  }

  return 0;
}


static int input_result_node(struct hecmwST_result_data *result, int n_node, FILE *fp) {
#define DELIM " \n"
  int i,rc,n;
  int label_counter;
  char *p,*buf;

  if(result->nn_component <= 0) {
    return 0;
  }

  /* nn_dof */
  result->nn_dof = HECMW_malloc(sizeof(*result->nn_dof)*result->nn_component);
  if(result->nn_dof == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
    return -1;
  }
  i = n = 0;
  buf = Line_Buf;
  while(i < result->nn_component) {
    p = strtok(buf, DELIM);
    if(p == NULL) {
      if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
        return -1;
      }
      buf = Line_Buf;
      continue;
    }
    buf = NULL;
    rc = sscanf(p, "%d", &result->nn_dof[i]);
    if(rc == EOF) {
      HECMW_set_error(HECMW_UTIL_E0204, "");
      return -1;
    }
    if(rc != 1) {
      HECMW_set_error(HECMW_UTIL_E0205, "nn_dof");
      return -1;
    }
    n += result->nn_dof[i];
    i++;
  }

  /* node_label */
  result->node_label = HECMW_malloc(sizeof(*result->node_label)*result->nn_component);
  if(result->node_label == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  for(i=0; i < result->nn_component; i++) {
    char label[HECMW_NAME_LEN+1];
    if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
      return -1;
    }
    rc = sscanf(Line_Buf, "%s", label);
    if(rc == EOF) {
      HECMW_set_error(HECMW_UTIL_E0204, "");
      return -1;
    }
    if(rc != 1) {
      HECMW_set_error(HECMW_UTIL_E0205, "node_label");
      return -1;
    }
    result->node_label[i] = HECMW_strdup(label);
    if(result->node_label[i] == NULL) {
      HECMW_set_error(errno, "");
      return -1;
    }
  }
  /* node_val_item */
  ResIO.node_global_ID = HECMW_malloc(sizeof(*ResIO.node_global_ID)*n_node);
  if(ResIO.node_global_ID == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }
  result->node_val_item = HECMW_malloc(sizeof(*result->node_val_item)*n*n_node);
  if(result->node_val_item == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
    return -1;
  }
  i = 0;
  label_counter = 0;
  n++;         /**** For global node ID ****/
  buf = Line_Buf;
  while(i < n*n_node) {
    p = strtok(buf, DELIM);
    if(p == NULL) {
      if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
        return -1;
      }
      buf = Line_Buf;
      continue;
    }
    buf = NULL;
    if ( (i%n) == 0 ) {
      rc = sscanf(p, "%d", &ResIO.node_global_ID[label_counter]);
      label_counter++;
    } else {
      rc = sscanf(p, "%lf", &result->node_val_item[i-label_counter]);
    }
    if(rc == EOF) {
      HECMW_set_error(HECMW_UTIL_E0204, "");
      return -1;
    }
    if(rc != 1) {
      HECMW_set_error(HECMW_UTIL_E0205, "node_val_item");
      return -1;
    }
    i++;
  }

  return 0;
}


static int input_result_elem(struct hecmwST_result_data *result, int n_elem, FILE *fp) {
#define DELIM " \n"
  int i,rc,n;
  int label_counter;
  char *p,*buf;

  if(result->ne_component <= 0) {
    return 0;
  }

  /* ne_dof */
  result->ne_dof = HECMW_malloc(sizeof(*result->ne_dof)*result->ne_component);
  if(result->ne_dof == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
    return -1;
  }
  i = n = 0;
  buf = Line_Buf;
  while(i < result->ne_component) {

    p = strtok(buf, DELIM);
    if(p == NULL) {
      if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
        return -1;
      }
      buf = Line_Buf;
      continue;
    }
    buf = NULL;
    rc = sscanf(p, "%d", &result->ne_dof[i]);
    if(rc == EOF) {
      HECMW_set_error(HECMW_UTIL_E0204, "");
      return -1;
    }
    if(rc != 1) {
      HECMW_set_error(HECMW_UTIL_E0205, "ne_dof");
      return -1;
    }
    n += result->ne_dof[i];
    i++;
  }

  /* elem_label */
  result->elem_label = HECMW_malloc(sizeof(*result->elem_label)*result->ne_component);
  if(result->elem_label == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  for(i=0; i < result->ne_component; i++) {
    char label[HECMW_NAME_LEN+1];
    if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
      return -1;
    }
    rc = sscanf(Line_Buf, "%s", label);
    if(rc == EOF) {
      HECMW_set_error(HECMW_UTIL_E0204, "");
      return -1;
    }
    if(rc != 1) {
      HECMW_set_error(HECMW_UTIL_E0205, "elem_label");
      return -1;
    }
    result->elem_label[i] = HECMW_strdup(label);
    if(result->elem_label[i] == NULL) {
      HECMW_set_error(errno, "");
      return -1;
    }
  }

  /* elem_val_item */
  ResIO.elem_global_ID = HECMW_malloc(sizeof(*ResIO.elem_global_ID)*n_elem);
  if(ResIO.elem_global_ID == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }
  result->elem_val_item = HECMW_malloc(sizeof(*result->elem_val_item)*n*n_elem);
  if(result->elem_val_item == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
    return -1;
  }
  i = 0;
  label_counter = 0;
  n++;          /**** For global element ID ****/
  buf = Line_Buf;
  while(i < n*n_elem) {
    p = strtok(buf, DELIM);
    if(p == NULL) {
      if(get_line(Line_Buf, sizeof(Line_Buf), fp) < 0) {
        return -1;
      }
      buf = Line_Buf;
      continue;
    }
    buf = NULL;
    if ( (i%n) == 0 ) {
      rc = sscanf(p, "%d", &ResIO.elem_global_ID[label_counter]);
      label_counter++;
    } else {
      rc = sscanf(p, "%lf", &result->elem_val_item[i-label_counter]);
    }
    if(rc == EOF) {
      HECMW_set_error(HECMW_UTIL_E0204, "");
      return -1;
    }
    if(rc != 1) {
      HECMW_set_error(HECMW_UTIL_E0205, "elem_val_item");
      return -1;
    }
    i++;
  }

  return 0;
}


static struct hecmwST_result_data *input_result_data(FILE *fp) {
  int n_node, n_elem;
  struct hecmwST_result_data *result;

  HECMW_assert(fp);

  result = HECMW_calloc(1, sizeof(*result));
  if(result == NULL) {
    HECMW_set_error(errno, "");
    return NULL;
  }
  if(input_result_header(result, fp)) {
    return NULL;
  }
  if( HECMW_RESULT_FILEVER_MAJOR > 1 ){
    if(input_result_global(result, fp)) {
      return NULL;
    }
  }
  if(input_result_dataheader(result, &n_node, &n_elem, fp)) {
    return NULL;
  }
  ResIO.nnode = n_node;
  ResIO.nelem = n_elem;
  if(input_result_node(result, n_node, fp)) {
    return NULL;
  }
  if(input_result_elem(result, n_elem, fp)) {
    return NULL;
  }

  return result;
}

/*---------------------------------------------------------------------------*/

struct hecmwST_result_data *HECMW_result_io_txt_read_by_fname(char *filename) {
  FILE *fp;
  struct hecmwST_result_data *result;

  if((fp = fopen(filename, "r")) == NULL) {
    HECMW_set_error(HECMW_UTIL_E0201, "File: %s, %s", filename, HECMW_strmsg(errno));
    return NULL;
  }

  result = input_result_data(fp);
  if(result == NULL) {
    return NULL;
  }

  if(fclose(fp)) {
    HECMW_set_error(HECMW_UTIL_E0202, "");
    return NULL;
  }

  return result;
}
