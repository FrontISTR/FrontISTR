/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "hecmw_util.h"
#include "hecmw_bin_io.h"
#include "hecmw_result.h"
#include "hecmw_result_io.h"

#define RES_BIN_HEADER "HECMW_BINARY_RESULT"

#define LINEBUF_SIZE 1023
static char Line_Buf[LINEBUF_SIZE + 1];


/*---------------------------------------------------------------------------*/
/* BINARY MODE I/O --- bin_header                                            */
/*---------------------------------------------------------------------------*/


static int write_bin_header(FILE* fp) {
  char* s = (char*)RES_BIN_HEADER;
  size_t n;
  char nbyte[3];

  n = strlen(s);
  if( fwrite( s, sizeof(char), n, fp) != n ) return -1;
  n = sizeof(long);
  sprintf( nbyte, "%2zd", n );
  if( fwrite( nbyte, sizeof(char), 2, fp) != 2 ) return -1;
  return 0;
}


static int check_bin_header(FILE* fp) {
  char* s = (char*)RES_BIN_HEADER;
  size_t n = strlen(s);
  char buff[256], nbyte[3];

  if( fread( buff, sizeof(char), n, fp) != n ) return 0;
  if( fread( nbyte, sizeof(char), 2, fp) != 2 ) return 0;

  buff[n] = 0;
  return ( strcmp( buff, s ) == 0 );
}


int HECMW_result_io_bin_judge_file(char *filename) {
  int rcode;
  FILE* fp;

  if((fp = fopen(filename, "rb")) == NULL) {
    HECMW_set_error(HECMW_UTIL_E0201, "File: %s, %s", filename, HECMW_strmsg(errno));
    return 0;
  }

  hecmw_set_endian_info();
  rcode = check_bin_header(fp);
  fclose(fp);

  return rcode;
}


/*---------------------------------------------------------------------------*/
/* BINARY MODE I/O --- bin_output_result                                     */
/*---------------------------------------------------------------------------*/


static int bin_output_result_header(FILE *fp) {
  int rc;

  /* header */
  if( HECMW_RESULT_FILEVER_MAJOR > 1 ){
    sprintf(ResIO.head,"%s %d.%d",ResIO.head,HECMW_RESULT_FILEVER_MAJOR,HECMW_RESULT_FILEVER_MINOR);
  }
  rc = hecmw_write_bin(fp,"S", ResIO.head);
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "head");
    return -1;
  }

  return 0;
}


static int bin_output_result_global(FILE *fp) {
  int i,j,k,n,rc,ng_comp;
  struct result_list *p,**data;

  /* comment */
  rc = hecmw_write_bin(fp,"S","*comment");
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "*comment");
    return -1;
  }
  rc = hecmw_write_bin(fp,"S", ResIO.comment_line);
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "head");
    return -1;
  }

  /* global header */
  rc = hecmw_write_bin(fp,"S","*global");
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "*global");
    return -1;
  }

  /* ng_component */
  rc = hecmw_write_bin(fp, "II", HECMW_result_io_count_ng_comp());
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "ng_comp");
    return -1;
  }

  /* ng_dof */
  n = 0;
  for(p=ResIO.global_list; p; p=p->next) {
    rc = hecmw_write_bin(fp, "I", p->n_dof );
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "ng_dof");
      return -1;
    }
    n++;
  }

  /* global_label */
  for(p=ResIO.global_list; p; p=p->next) {
    rc = hecmw_write_bin(fp, "S", p->label);
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
  for(j=0; j < ng_comp; j++) {
    p = data[j];
    for(k=0; k < p->n_dof; k++) {
      rc = hecmw_write_bin(fp, "F", p->ptr[k] );
      if(rc < 0) {
        HECMW_set_error(HECMW_UTIL_E0205, "global_val_item");
        return -1;
      }
    }
  }
  HECMW_free(data);

  return 0;
}


static int bin_output_result_dataheader(FILE *fp) {
  int rc;

  /* n_node, n_elem */
  rc = hecmw_write_bin(fp, "II", ResIO.nnode, ResIO.nelem);
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "nnode,nelem");
    return -1;
  }

  /* nn_component, ne_component */
  rc = hecmw_write_bin(fp, "II", HECMW_result_io_count_nn_comp(), HECMW_result_io_count_ne_comp());
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "nn_comp,ne_comp");
    return -1;
  }

  return 0;
}


static int bin_output_result_node(FILE *fp) {
  int i,j,k,n,rc,nn_comp;
  struct result_list *p,**data;

  /* nn_dof */
  n = 0;
  for(p=ResIO.node_list; p; p=p->next) {
    rc = hecmw_write_bin(fp, "I", p->n_dof );
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "nn_dof");
      return -1;
    }
    n++;
  }

  /* node_label */
  for(p=ResIO.node_list; p; p=p->next) {
    rc = hecmw_write_bin(fp, "S", p->label);
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
    rc = hecmw_write_bin(fp, "I", ResIO.node_global_ID[i] );
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "node_global_ID");
      return -1;
    }
    for(j=0; j < nn_comp; j++) {
      p = data[j];
      for(k=0; k < p->n_dof; k++) {
        rc = hecmw_write_bin(fp, "F", p->ptr[i*p->n_dof+k] );
        if(rc < 0) {
          HECMW_set_error(HECMW_UTIL_E0205, "node_val_item");
          return -1;
        }
      }
    }
  }
  HECMW_free(data);

  return 0;
}


static int bin_output_result_elem(FILE *fp) {
  int i,j,k,n,rc,ne_comp;
  struct result_list *p,**data;

  /* ne_dof */
  n = 0;
  for(p=ResIO.elem_list; p; p=p->next) {
    rc = hecmw_write_bin(fp, "I", p->n_dof );
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "ne_dof");
      return -1;
    }
    n++;
  }

  /* elem_label */
  for(p=ResIO.elem_list; p; p=p->next) {
    rc = hecmw_write_bin(fp, "S", p->label);
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
    rc = hecmw_write_bin(fp, "I", ResIO.elem_global_ID[i] );
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "elem_global_ID");
      return -1;
    }
    for(j=0; j < ne_comp; j++) {
      p = data[j];
      for(k=0; k < p->n_dof; k++) {
        rc = hecmw_write_bin(fp, "F", p->ptr[i*p->n_dof+k] );
        if(rc < 0) {
          HECMW_set_error(HECMW_UTIL_E0205, "elem_val_item");
          return -1;
        }
      }
    }
  }
  HECMW_free(data);

  return 0;
}


static int bin_output_result_data(FILE *fp) {
  int rc;
  HECMW_assert(fp);

  if(bin_output_result_header(fp)) {
    return -1;
  }
  if( HECMW_RESULT_FILEVER_MAJOR > 1 ){
    if(bin_output_result_global(fp)) {
      return -1;
    }
    /* data header */
    rc = hecmw_write_bin(fp,"S","*data");
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "head");
      return -1;
    }
  }
  if(bin_output_result_dataheader(fp)) {
    return -1;
  }
  if(bin_output_result_node(fp)) {
    return -1;
  }
  if(bin_output_result_elem(fp)) {
    return -1;
  }

  return 0;
}

/*---------------------------------------------------------------------------*/

int HECMW_result_io_bin_write_by_fname(char *filename) {
  FILE *fp = NULL;

  if (HECMW_ctrl_is_subdir()) {
    if (HECMW_ctrl_make_subdir(filename)) {
      HECMW_set_error(HECMW_UTIL_E0201, "File: %s, %s", filename,
                      HECMW_strmsg(errno));
      goto error;
    }
  }

  if ((fp = fopen(filename, "wb")) == NULL) {
    HECMW_set_error(HECMW_UTIL_E0201, "File: %s, %s", filename,
                    HECMW_strmsg(errno));
    goto error;
  }

  hecmw_set_endian_info();
  if (write_bin_header(fp)) goto error;
  if (bin_output_result_data(fp)) goto error;

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
/* BINARY MODE I/O --- bin_output_result_ST                                  */
/*---------------------------------------------------------------------------*/


static int bin_output_result_header_ST(struct hecmwST_result_data *result,
                                       char *header, FILE *fp) {
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
  rc = hecmw_write_bin(fp, "S", header);
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "header");
    return -1;
  }

  return 0;
}


static int bin_output_result_global_ST(struct hecmwST_result_data *result,
                                       char *comment, FILE *fp) {
  size_t len;
  int i,j,k,n,m,rc;
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
  rc = hecmw_write_bin(fp,"S","*comment");
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "*comment");
    return -1;
  }
  rc = hecmw_write_bin(fp, "S", comment);
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "comment");
    return -1;
  }

  /* global header */
  rc = hecmw_write_bin(fp,"S","*global");
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "*global");
    return -1;
  }

  /* ng_component */
  rc = hecmw_write_bin(fp, "II", result->ng_component);
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "ng_comp");
    return -1;
  }

  /* ng_dof */
  n = 0;
  for(i=0; i < result->ng_component; i++) {
    rc = hecmw_write_bin(fp, "I", result->ng_dof[i] );
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "ng_dof");
      return -1;
    }
    n++;
  }

  /* global_label */
  for(i=0; i < result->ng_component; i++) {
    rc = hecmw_write_bin(fp, "S", result->global_label[i]);
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "global_label");
      return -1;
    }
  }

  /* global_val_item */
  if(result->ng_component == 0) return 0;
  m = 0;
  for(j=0; j < result->ng_component; j++) {
    for(k=0; k < result->ng_dof[j]; k++) {
      rc = hecmw_write_bin(fp, "F", result->global_val_item[m] );
      if(rc < 0) {
        HECMW_set_error(HECMW_UTIL_E0205, "global_val_item");
        return -1;
      }
      m++;
    }
  }

  /* data header */
  rc = hecmw_write_bin(fp,"S","*data");
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "*data");
    return -1;
  }

  return 0;
}


static int bin_output_result_dataheader_ST(struct hecmwST_result_data *result,
                                           int n_node, int n_elem, FILE *fp) {
  int rc;

  /* n_node, n_elem */
  rc = hecmw_write_bin(fp, "II", n_node, n_elem);
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "n_node,n_elem");
    return -1;
  }

  /* nn_component, ne_component */
  rc = hecmw_write_bin(fp, "II", result->nn_component, result->ne_component);
  if(rc < 0) {
    HECMW_set_error(HECMW_UTIL_E0205, "nn_comp,ne_comp");
    return -1;
  }

  return 0;
}


static int bin_output_result_node_ST(struct hecmwST_result_data *result,
                                     int n_node, FILE *fp) {
  int i,j,k,n,m,rc;

  /* nn_dof */
  n = 0;
  for(i=0; i < result->nn_component; i++) {
    rc = hecmw_write_bin(fp, "I", result->nn_dof[i] );
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "nn_dof");
      return -1;
    }
    n++;
  }

  /* node_label */
  for(i=0; i < result->nn_component; i++) {
    rc = hecmw_write_bin(fp, "S", result->node_label[i]);
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "node_label");
      return -1;
    }
  }

  /* node_val_item */
  if(result->nn_component == 0) return 0;
  m = 0;
  for(i=0; i < n_node; i++) {
    rc = hecmw_write_bin(fp, "I", ResIO.node_global_ID[i] );
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "node_global_ID");
      return -1;
    }
    for(j=0; j < result->nn_component; j++) {
      for(k=0; k < result->nn_dof[j]; k++) {
        rc = hecmw_write_bin(fp, "F", result->node_val_item[m] );
        if(rc < 0) {
          HECMW_set_error(HECMW_UTIL_E0205, "node_val_item");
          return -1;
        }
        m++;
      }
    }
  }

  return 0;
}


static int bin_output_result_elem_ST(struct hecmwST_result_data *result,
                                     int n_elem, FILE *fp) {
  int i,j,k,n,m,rc;

  /* ne_dof */
  n = 0;
  for(i=0; i < result->ne_component; i++) {
    rc = hecmw_write_bin(fp, "I", result->ne_dof[i] );
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "ne_dof");
      return -1;
    }
    n++;
  }

  /* elem_label */
  for(i=0; i < result->ne_component; i++) {
    rc = hecmw_write_bin(fp, "S", result->elem_label[i]);
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "elem_label");
      return -1;
    }
  }

  /* elem_val_item */
  if(result->ne_component == 0) return 0;
  m = 0;
  for(i=0; i < n_elem; i++) {
    rc = hecmw_write_bin(fp, "I", ResIO.elem_global_ID[i] );
    if(rc < 0) {
      HECMW_set_error(HECMW_UTIL_E0205, "elem_global_ID");
      return -1;
    }
    for(j=0; j < result->ne_component; j++) {
      for(k=0; k < result->ne_dof[j]; k++) {
        rc = hecmw_write_bin(fp, "F", result->elem_val_item[m]);
        if(rc < 0) {
          HECMW_set_error(HECMW_UTIL_E0205, "elem_val_item");
          return -1;
        }
        m++;
      }
    }
  }

  return 0;
}


static int bin_output_result_data_ST(struct hecmwST_result_data *result,
                                     int n_node, int n_elem, char *header,
                                     char *comment, FILE *fp) {
  HECMW_assert(fp);

  if(bin_output_result_header_ST(result, header, fp)) {
    return -1;
  }
  if( HECMW_RESULT_FILEVER_MAJOR > 1 ){
    if(bin_output_result_global_ST(result, comment, fp)) {
      return -1;
    }
  }
  if(bin_output_result_dataheader_ST(result, n_node, n_elem, fp)) {
    return -1;
  }
  if(bin_output_result_node_ST(result, n_node, fp)) {
    return -1;
  }
  if(bin_output_result_elem_ST(result, n_elem, fp)) {
    return -1;
  }

  return 0;
}

/*---------------------------------------------------------------------------*/

int HECMW_result_io_bin_write_ST_by_fname(char *filename,
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

  if ((fp = fopen(filename, "wb")) == NULL) {
    HECMW_set_error(HECMW_UTIL_E0201, "File: %s, %s", filename,
                    HECMW_strmsg(errno));
    goto error;
  }

  hecmw_set_endian_info();
  if (write_bin_header(fp)) goto error;
  if (bin_output_result_data_ST(result, n_node, n_elem, header, comment, fp)) goto error;

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
/* BINARY MODE I/O --- bin_input_result                                      */
/*---------------------------------------------------------------------------*/


static int bin_input_result_header(struct hecmwST_result_data *result, FILE *fp) {
  char *ptr;

  /* header */
  if(hecmw_read_bin(fp, "S", Line_Buf)) {
    HECMW_set_error(HECMW_UTIL_E0205, "header");
    return -1;
  }
  if( HECMW_RESULT_FILEVER_MAJOR > 1 ){
    ptr = strtok(Line_Buf, " ");
    sprintf(Line_Buf, "%s", ptr);
  }
  strcpy( ResIO.head, Line_Buf );

  return 0;
}


static int bin_input_result_global(struct hecmwST_result_data *result, FILE *fp) {
  int i,j,k,n,m;
  /* comment */
  if(hecmw_read_bin(fp, "S", Line_Buf)) { //skip comment header
    HECMW_set_error(HECMW_UTIL_E0205, "comment");
    return -1;
  }
  if(hecmw_read_bin(fp, "S", Line_Buf)) {
    HECMW_set_error(HECMW_UTIL_E0205, "comment");
    return -1;
  }
  strcpy( ResIO.comment_line, Line_Buf );

  /* skip global header */
  if(hecmw_read_bin(fp, "S", Line_Buf)) {
    HECMW_set_error(HECMW_UTIL_E0205, "*global");
    return -1;
  }

  /* ng_component */
  if(hecmw_read_bin(fp, "II", &result->ng_component)) {
    HECMW_set_error(HECMW_UTIL_E0205, "ng_component");
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

  n = 0;
  for(i=0; i<result->ng_component; i++){
    if(hecmw_read_bin( fp, "I", &result->ng_dof[i] )) {
      HECMW_set_error(HECMW_UTIL_E0205, "ng_dof");
      return -1;
    }
    n += result->ng_dof[i];
  }

  /* global_label */
  result->global_label = HECMW_malloc(sizeof(*result->global_label)*result->ng_component);
  if(result->global_label == NULL) {
    HECMW_set_error(errno, "(global_label)");
    return -1;
  }

  for(i=0; i < result->ng_component; i++) {
    char label[HECMW_NAME_LEN+1];
    if(hecmw_read_bin( fp, "S", label )) {
      HECMW_set_error(HECMW_UTIL_E0205, "global_label");
      return -1;
    }
    result->global_label[i] = HECMW_strdup(label);
    if(result->global_label[i] == NULL) {
      HECMW_set_error(errno, "(label)");
      return -1;
    }
  }

  /* global_val_item */
  result->global_val_item = HECMW_malloc(sizeof(*result->global_val_item)*n);
  if(result->global_val_item == NULL) {
    HECMW_set_error(errno, "(global_val_item)");
    return -1;
  }

  m = 0;
  for(j=0; j < result->ng_component; j++) {
    for(k=0; k < result->ng_dof[j]; k++) {
      if(hecmw_read_bin( fp, "F", &result->global_val_item[m] )) {
        HECMW_set_error(HECMW_UTIL_E0205, "global_val_item");
        return -1;
      }
      m++;
    }
  }

  /* skip data header */
  if(hecmw_read_bin(fp, "S", Line_Buf)) {
    HECMW_set_error(HECMW_UTIL_E0205, "*data");
    return -1;
  }

  return 0;
}


static int bin_input_result_dataheader(struct hecmwST_result_data *result,
                                       int *n_node, int *n_elem, FILE *fp) {
  int nn, ne;

  /* n_node, n_elem */
  if(hecmw_read_bin(fp, "II", &nn, &ne)) {
    HECMW_set_error(HECMW_UTIL_E0205, "n_node,n_elem");
    return -1;
  }
  *n_node = nn;
  *n_elem = ne;

  /* nn_component, ne_component */
  if(hecmw_read_bin(fp, "II", &result->nn_component, &result->ne_component)) {
    HECMW_set_error(HECMW_UTIL_E0205, "nn_comp,ne_comp");
    return -1;
  }

  return 0;
}


static int bin_input_result_node(struct hecmwST_result_data *result, int n_node, FILE *fp) {
  int i,j,k,n,m;

  if(result->nn_component <= 0) {
    return 0;
  }

  /* nn_dof */
  result->nn_dof = HECMW_malloc(sizeof(*result->nn_dof)*result->nn_component);
  if(result->nn_dof == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  n = 0;
  for(i=0; i<result->nn_component; i++){
    if(hecmw_read_bin( fp, "I", &result->nn_dof[i] )) {
      HECMW_set_error(HECMW_UTIL_E0205, "nn_dof");
      return -1;
    }
    n += result->nn_dof[i];
  }

  /* node_label */
  result->node_label = HECMW_malloc(sizeof(*result->node_label)*result->nn_component);
  if(result->node_label == NULL) {
    HECMW_set_error(errno, "(node_label)");
    return -1;
  }

  for(i=0; i < result->nn_component; i++) {
    char label[HECMW_NAME_LEN+1];
    if(hecmw_read_bin( fp, "S", label )) {
      HECMW_set_error(HECMW_UTIL_E0205, "node_label");
      return -1;
    }
    result->node_label[i] = HECMW_strdup(label);
    if(result->node_label[i] == NULL) {
      HECMW_set_error(errno, "(label)");
      return -1;
    }
  }

  /* node_val_item */
  ResIO.node_global_ID = HECMW_malloc(sizeof(*ResIO.node_global_ID)*n_node);
  if(ResIO.node_global_ID == NULL) {
    HECMW_set_error(errno, "(node_global_ID)");
    return -1;
  }
  result->node_val_item = HECMW_malloc(sizeof(*result->node_val_item)*n*n_node);
  if(result->node_val_item == NULL) {
    HECMW_set_error(errno, "(node_val_item)");
    return -1;
  }

  m = 0;
  for(i=0; i < n_node; i++) {
    if(hecmw_read_bin( fp, "I", &ResIO.node_global_ID[i] )) {
      HECMW_set_error(HECMW_UTIL_E0205, "node_global_ID");
      return -1;
    }
    for(j=0; j < result->nn_component; j++) {
      for(k=0; k < result->nn_dof[j]; k++) {
        if(hecmw_read_bin( fp, "F", &result->node_val_item[m] )) {
          HECMW_set_error(HECMW_UTIL_E0205, "node_val_item");
          return -1;
        }
        m++;
      }
    }

  }

  return 0;
}


static int bin_input_result_elem(struct hecmwST_result_data *result, int n_elem, FILE *fp) {
  int i,j,k,n,m;

  if(result->ne_component <= 0) {
    return 0;
  }

  /* ne_dof */
  result->ne_dof = HECMW_malloc(sizeof(*result->ne_dof)*result->ne_component);
  if(result->ne_dof == NULL) {
    HECMW_set_error(errno, "(ne_dof)");
    return -1;
  }

  n = 0;
  for(i=0; i<result->ne_component;i++ ){
    if(hecmw_read_bin( fp, "I", &result->ne_dof[i] )) {
      HECMW_set_error(HECMW_UTIL_E0205, "ne_dof");
      return -1;
    }
    n += result->ne_dof[i];
  }

  /* elem_label */
  result->elem_label = HECMW_malloc(sizeof(*result->elem_label)*result->ne_component);
  if(result->elem_label == NULL) {
    HECMW_set_error(errno, "(elem_label)");
    return -1;
  }

  for(i=0; i < result->ne_component; i++) {
    char label[HECMW_NAME_LEN+1];
    if(hecmw_read_bin( fp, "S", label )) {
      HECMW_set_error(HECMW_UTIL_E0205, "elem_label");
      return -1;
    }
    result->elem_label[i] = HECMW_strdup(label);
    if(result->elem_label[i] == NULL) {
      HECMW_set_error(errno, "(label)");
      return -1;
    }
  }

  /* elem_val_item */
  ResIO.elem_global_ID = HECMW_malloc(sizeof(*ResIO.elem_global_ID)*n_elem);
  if(ResIO.elem_global_ID == NULL) {
    HECMW_set_error(errno, "(elem_global_ID)");
    return -1;
  }
  result->elem_val_item = HECMW_malloc(sizeof(*result->elem_val_item)*n*n_elem);
  if(result->elem_val_item == NULL) {
    HECMW_set_error(errno, "(elem_val_item)");
    return -1;
  }

  m = 0;
  for(i=0; i < n_elem; i++) {
    if(hecmw_read_bin( fp, "I", &ResIO.elem_global_ID[i] )) {
      HECMW_set_error(HECMW_UTIL_E0205, "elem_global_ID");
      return -1;
    }
    for(j=0; j < result->ne_component; j++) {
      for(k=0; k < result->ne_dof[j]; k++) {
        if(hecmw_read_bin( fp, "F", &result->elem_val_item[m] )) {
          HECMW_set_error(HECMW_UTIL_E0205, "elem_val_item");
          return -1;
        }
        m++;
      }
    }

  }

  return 0;
}


static struct hecmwST_result_data *bin_input_result_data(FILE *fp) {
  int n_node, n_elem;
  struct hecmwST_result_data *result;

  HECMW_assert(fp);

  result = HECMW_calloc(1, sizeof(*result));
  if(result == NULL) {
    HECMW_set_error(errno, "");
    return NULL;
  }

  if(bin_input_result_header(result, fp)) {
    return NULL;
  }

  if( HECMW_RESULT_FILEVER_MAJOR > 1 ){
    if(bin_input_result_global(result, fp)) {
      return NULL;
    }
  }
  if(bin_input_result_dataheader(result, &n_node, &n_elem, fp)) {
    return NULL;
  }
  ResIO.nnode = n_node;
  ResIO.nelem = n_elem;

  if(bin_input_result_node(result, n_node, fp)) {
    return NULL;
  }

  if(bin_input_result_elem(result, n_elem, fp)) {
    return NULL;
  }

  return result;
}


/*---------------------------------------------------------------------------*/


struct hecmwST_result_data *HECMW_result_io_bin_read_by_fname(char *filename) {
  FILE *fp;
  struct hecmwST_result_data *result;

  if((fp = fopen(filename, "rb")) == NULL) {
    HECMW_set_error(HECMW_UTIL_E0201, "File: %s, %s", filename, HECMW_strmsg(errno));
    return NULL;
  }

  hecmw_set_endian_info();

  if(!check_bin_header(fp)) {
    fclose(fp);
    HECMW_set_error(HECMW_UTIL_E0202, "%s is not binary result file", filename);
    return NULL;
  }
  result = bin_input_result_data(fp);
  if(result == NULL) {
    return NULL;
  }
  if(fclose(fp)) {
    HECMW_set_error(HECMW_UTIL_E0202, "");
    return NULL;
  }

  return result;
}
