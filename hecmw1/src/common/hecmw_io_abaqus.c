/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "hecmw_util.h"
#include "hecmw_ablex.h"
#include "hecmw_io_abaqus.h"
#include "hecmw_io_mesh.h"
#include "hecmw_io_struct.h"
#include "hecmw_struct.h"
#include "hecmw_config.h"
#include "hecmw_system.h"
#include "hecmw_dist.h"
#include "hecmw_dist_print.h"
#include "hecmw_common.h"
#include "hecmw_path.h"
#include "hecmw_conn_conv.h"
#include "hecmw_map_int.h"

/*----------------------------------------------------------------------------*/

static char grid_filename[HECMW_FILENAME_LEN + 1]    = "Unknown";
static char include_filename[HECMW_FILENAME_LEN + 1] = "Unknown";

/*----------------------------------------------------------------------------*/

static void do_logging(int loglv, int msgno, int add_location, const char *fmt,
                       va_list ap) {
  char line[100] = "";
  char msg[HECMW_MSG_LEN + 1];
  char *p;

  HECMW_vsnprintf(msg, sizeof(msg), fmt, ap);
  if (add_location) {
    char *s                = "";
    if (strlen(msg) > 0) s = ": ";
    p = HECMW_ablex_is_including() ? include_filename : grid_filename;
    HECMW_snprintf(line, sizeof(line), "%s:%d%s", p, HECMW_ablex_get_lineno(),
                   s);
  }
  if (loglv == HECMW_LOG_ERROR) {
    HECMW_set_error(msgno, "%s%s", line, msg);
  } else {
    HECMW_print_msg(loglv, msgno, "%s%s", line, msg);
  }
}

static void set_err_noloc(int msgno, const char *fmt, ...) {
  va_list ap;

  va_start(ap, fmt);
  do_logging(HECMW_LOG_ERROR, msgno, 0, fmt, ap);
  va_end(ap);
}

static void set_err(int msgno, const char *fmt, ...) {
  va_list ap;

  va_start(ap, fmt);
  do_logging(HECMW_LOG_ERROR, msgno, 1, fmt, ap);
  va_end(ap);
}

static void set_err_token(int token, int msgno, const char *fmt, ...) {
  int msg_no;
  va_list ap;

  if (!token) {
    msg_no = HECMW_IO_ABAQUS_E0003;
  } else {
    msg_no = msgno;
  }
  va_start(ap, fmt);
  do_logging(HECMW_LOG_ERROR, msg_no, 1, fmt, ap);
  va_end(ap);
}

static void log_warn(int msgno, const char *fmt, ...) {
  va_list ap;

  va_start(ap, fmt);
  do_logging(HECMW_LOG_WARN, msgno, 1, fmt, ap);
  va_end(ap);
}

/*----------------------------------------------------------------------------*/

typedef struct { int i; } Integer;

static struct hecmw_map_int *elem_secopt;

static void free_Integer(void *v) {
  Integer *i;
  i = (Integer *)v;
  /* do nothing on members */
  HECMW_free(i);
}

/*----------------------------------------------------------------------------*/

static struct material_keywords {
  int keyword;
  char *string;
} material_keywords[] = {
    {HECMW_ABLEX_H_CONDUCTIVITY, "*CONDUCTIVITY"},
    {HECMW_ABLEX_H_DENSITY, "*DENSITY"},
    {HECMW_ABLEX_H_ELASTIC, "*ELASTIC"},
    {HECMW_ABLEX_H_SPECIFIC_HEAT, "*SPECIFIC HEAT"},
};

#define N_MAT_KEYS (sizeof(material_keywords) / sizeof(material_keywords[0]))

static int is_material_keyword(int keyword) {
  int i;

  for (i = 0; i < N_MAT_KEYS; i++) {
    if (keyword == material_keywords[i].keyword) return 1;
  }
  return 0;
}

static char *get_material_string(int keyword) {
  int i;

  for (i = 0; i < N_MAT_KEYS; i++) {
    if (keyword == material_keywords[i].keyword) {
      return material_keywords[i].string;
    }
  }
  return NULL;
}

/*----------------------------------------------------------------------------*/

static int flag_material_zone;
/* true: in *MATERIAL to next *MATERIAL or not material keyword */
static char matname[HECMW_NAME_LEN + 1] = "";
static struct material_data {
  int keyword;
  struct hecmw_io_matitem *matitem;
  struct material_data *next;
} * matdata;

static int is_material_zone(void) { return flag_material_zone; }

static void set_material_zone(int flag) { flag_material_zone = flag ? 1 : 0; }

static int add_mat_data(int keyword, struct hecmw_io_matitem *matitem) {
  int i;
  struct material_data *p, *q, *mdata;

  for (p = matdata; p; p = p->next) {
    if (p->keyword == keyword) break;
  }
  if (p) {
    struct hecmw_io_matsubitem *sip, *siq;
    struct hecmw_io_matitem *item = p->matitem;
    p->matitem                    = matitem; /* update */
    matitem->item                 = item->item;
    for (sip = item->subitem; sip; sip = siq) {
      siq = sip->next;
      HECMW_free(sip->val);
      HECMW_free(sip);
    }
    HECMW_free(item);
    log_warn(HECMW_IO_ABAQUS_W0095, "%s updated for *MATERIAL %s",
             get_material_string(keyword), matname);
    return 0;
  }

  mdata = HECMW_malloc(sizeof(*mdata));
  if (mdata == NULL) {
    set_err(errno, "");
    return -1;
  }

  mdata->keyword = keyword;
  mdata->matitem = matitem;
  mdata->next    = NULL;

  q = NULL;
  for (i = 0, p = matdata; p; p = (q = p)->next, i++)
    ;
  matitem->item = i + 1;
  if (q == NULL) {
    matdata = mdata;
  } else {
    q->next = mdata;
  }

  return 0;
}

static int count_mat_item(void) {
  int n;
  struct material_data *p, *q;

  for (n = 0, p = matdata; p; p = (q = p)->next, n++)
    ;
  return n;
}

static int regist_material(void) {
  int i, n;
  struct hecmw_io_material *mat = NULL;
  struct material_data *p, *q;

  n = count_mat_item();
  if (n == 0) return 0;

  mat = HECMW_calloc(1, sizeof(*mat) * n);
  if (mat == NULL) {
    set_err(errno, "");
    goto error;
  }

  HECMW_assert(strlen(matname) > 0);

  mat->nitem = n;
  strcpy(mat->name, matname);

  mat->item = HECMW_malloc(sizeof(*mat->item) * n);
  if (mat->item == NULL) {
    set_err(errno, "");
    goto error;
  }

  for (i = 0, p = matdata; p; p = q, i++) {
    q = p->next;
    HECMW_assert(p->matitem->item == i + 1);
    mat->item[i] = *p->matitem;
    HECMW_free(p->matitem);
    HECMW_free(p);
  }

  if (HECMW_io_add_mat(matname, mat) == NULL) goto error;

  strcpy(matname, "");
  matdata = NULL;

  return 0;
error:
  if (mat) {
    HECMW_free(mat->item);
    HECMW_free(mat);
  }
  return -1;
}

static int read_mat_data_line_common(int nval, int nval_line, double *val,
                                     int errmsgno) {
  int i, token;

  HECMW_assert(nval > 0);
  HECMW_assert(nval_line > 0);
  HECMW_assert(val);

  for (i = 0; i < nval; i++) {
    val[i] = 0.0;
  }
  i = 0;
  while (1) {
    token = HECMW_ablex_next_token();
    if (token != HECMW_ABLEX_DOUBLE && token != HECMW_ABLEX_INT &&
        token != ',' && token != HECMW_ABLEX_NL) {
      set_err_token(token, errmsgno, "VAL or ',' or NL reuqired");
      return -1;
    }
    if (token == HECMW_ABLEX_NL) {
      HECMW_ablex_unput_token();
    } else if (token == ',') {
      i++;
      HECMW_ablex_unput_token();
    } else {
      val[i++] = HECMW_ablex_get_number();
    }

    token = HECMW_ablex_next_token();
    if (token != ',' && token != HECMW_ABLEX_NL) {
      set_err_token(token, errmsgno, "',' or NL required after VAL");
      return -1;
    }

    if (token == HECMW_ABLEX_NL) {
      if (i % nval_line) {
        i += nval_line - (i % nval_line);
      }
    }
    if (i >= nval) {
      if (token == ',') {
        token = HECMW_ablex_next_token();
        if (token != HECMW_ABLEX_NL) {
          set_err_token(token, errmsgno, "NL required");
          return -1;
        }
      }
      break;
    }
  }
  return 0;
}

static int read_mat_data_common(int nval, int nval_line, int temp_col,
                                struct hecmw_io_matitem **matitem, int msgno) {
  int token;
  struct hecmw_io_matitem *item       = NULL;
  struct hecmw_io_matsubitem *subitem = NULL;
  struct hecmw_io_matsubitem *q, *p = NULL;

  item = HECMW_malloc(sizeof(*item));
  if (item == NULL) {
    set_err(errno, "");
    goto error;
  }

  item->item    = -1; /* I don't know */
  item->nval    = nval;
  item->subitem = NULL;

  while (1) {
    subitem = HECMW_calloc(1, sizeof(*subitem));
    if (subitem == NULL) {
      set_err(errno, "");
      goto error;
    }

    subitem->val = HECMW_malloc(sizeof(*subitem->val) * nval);
    if (subitem->val == NULL) {
      set_err(errno, "");
      goto error;
    }

    if (read_mat_data_line_common(nval, nval_line, subitem->val, msgno))
      goto error;

    subitem->temp = subitem->val[temp_col - 1];
    subitem->next = NULL;

    if (p == NULL) {
      item->subitem = subitem;
    } else {
      p->next = subitem;
    }
    p = subitem;

    token = HECMW_ablex_next_token();
    HECMW_ablex_unput_token();
    if (token != HECMW_ABLEX_DOUBLE && token != HECMW_ABLEX_INT &&
        token != ',') {
      break;
    }
  }
  *matitem = item;

  return 0;
error:
  if (item) {
    for (p = item->subitem; p; p = q) {
      q = p->next;
      HECMW_free(p->val);
      HECMW_free(p);
    }
    HECMW_free(item);
  }
  return -1;
}

/*----------------------------------------------------------------------------*/

static int read_param_dependencies(int *dependencies, int msgno) {
  int token;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, msgno, "'=' required after DEPENDENCIES");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_INT) {
    set_err_token(token, msgno, "Invalid DEPENDENCIES");
    return -1;
  }
  *dependencies = HECMW_ablex_get_number();
  if (*dependencies <= 0) {
    set_err_token(token, msgno, "DEPENDENCIES must be positive integer");
    return -1;
  }
  return 0;
}

/*----------------------------------------------------------------------------*/

static struct etype_conv {
  int abaqus_etype;
  int hecmw_etype;
  int secopt;
} etype_conv[] = {
    {HECMW_ABLEX_E_B31, 611, 0},     {HECMW_ABLEX_E_B32, 612, 0},
    {HECMW_ABLEX_E_C3D4, 341, 0},    {HECMW_ABLEX_E_C3D6, 351, 0},
    {HECMW_ABLEX_E_C3D8, 361, 0},    {HECMW_ABLEX_E_C3D8I, 361, 0},
    {HECMW_ABLEX_E_C3D10, 342, 0},   {HECMW_ABLEX_E_C3D15, 352, 0},
    {HECMW_ABLEX_E_C3D20, 362, 0},   {HECMW_ABLEX_E_CAX3, 231, 2},
    {HECMW_ABLEX_E_CAX4, 241, 2},    {HECMW_ABLEX_E_CAX4I, 241, 2},
    {HECMW_ABLEX_E_CAX4R, 241, 12},  {HECMW_ABLEX_E_CAX6, 232, 2},
    {HECMW_ABLEX_E_CAX8, 242, 2},    {HECMW_ABLEX_E_CAX8R, 242, 12},
    {HECMW_ABLEX_E_CPE3, 231, 1},    {HECMW_ABLEX_E_CPE4, 241, 1},
    {HECMW_ABLEX_E_CPE4I, 241, 1},   {HECMW_ABLEX_E_CPE4R, 241, 11},
    {HECMW_ABLEX_E_CPE6, 232, 1},    {HECMW_ABLEX_E_CPE8, 242, 1},
    {HECMW_ABLEX_E_CPE8R, 242, 11},  {HECMW_ABLEX_E_CPS3, 231, 0},
    {HECMW_ABLEX_E_CPS4, 241, 0},    {HECMW_ABLEX_E_CPS4I, 241, 0},
    {HECMW_ABLEX_E_CPS4R, 241, 10},  {HECMW_ABLEX_E_CPS6, 232, 0},
    {HECMW_ABLEX_E_CPS8, 242, 0},    {HECMW_ABLEX_E_CPS8R, 242, 10},
    {HECMW_ABLEX_E_DC1D2, 111, 0},   {HECMW_ABLEX_E_DC1D3, 112, 0},
    {HECMW_ABLEX_E_DC2D3, 231, 0},   {HECMW_ABLEX_E_DC2D4, 241, 0},
    {HECMW_ABLEX_E_DC2D6, 232, 0},   {HECMW_ABLEX_E_DC2D8, 242, 0},
    {HECMW_ABLEX_E_DC3D4, 341, 0},   {HECMW_ABLEX_E_DC3D6, 351, 0},
    {HECMW_ABLEX_E_DC3D8, 361, 0},   {HECMW_ABLEX_E_DC3D10, 342, 0},
    {HECMW_ABLEX_E_DC3D15, 352, 0},  {HECMW_ABLEX_E_DC3D20, 362, 0},
    {HECMW_ABLEX_E_DCAX3, 231, 2},   {HECMW_ABLEX_E_DCAX4, 241, 2},
    {HECMW_ABLEX_E_DCAX6, 232, 0},   {HECMW_ABLEX_E_DCAX8, 242, 0},
    {HECMW_ABLEX_E_DINTER4, 541, 0}, {HECMW_ABLEX_E_DINTER8, 542, 0},
    {HECMW_ABLEX_E_DS4, 741, 0},     {HECMW_ABLEX_E_DS8, 742, 0},
    {HECMW_ABLEX_E_INTER4, 541, 0},  {HECMW_ABLEX_E_INTER8, 542, 0},
    {HECMW_ABLEX_E_S3R, 731, 0},     {HECMW_ABLEX_E_S4R, 741, 0},
    {HECMW_ABLEX_E_S8R, 742, 0},     {HECMW_ABLEX_E_T3D2, 111, 0},
    {HECMW_ABLEX_E_T3D3, 112, 0},
};

static int get_HECMW_etype(int abaqus_etype) {
  int i;

  for (i = 0; i < sizeof(etype_conv) / sizeof(etype_conv[0]); i++) {
    if (etype_conv[i].abaqus_etype == abaqus_etype) {
      return etype_conv[i].hecmw_etype;
    }
  }
  return -1;
}

static int get_secopt_abaqus(int abaqus_etype) {
  int i;

  for (i = 0; i < sizeof(etype_conv) / sizeof(etype_conv[0]); i++) {
    if (etype_conv[i].abaqus_etype == abaqus_etype) {
      return etype_conv[i].secopt;
    }
  }
  return -1;
}

static int get_secopt(char *elset, int msgno) {
  int i, secopt_prev;
  struct hecmw_io_id_array *elem = NULL;

  elem = HECMW_io_get_elem_in_egrp(elset);
  if (elem == NULL) return -1;

  secopt_prev = -1;
  for (i = 0; i < elem->n; i++) {
    Integer *secopt = HECMW_map_int_get(elem_secopt, elem->id[i]);
    if (secopt == NULL) {
      set_err_noloc(msgno, "");
      goto error;
    }
    if (i == 0) {
      secopt_prev = secopt->i;
    } else {
      if (secopt_prev != secopt->i) goto error;
    }
  }
  HECMW_free(elem->id);
  HECMW_free(elem);
  return secopt_prev;
error:
  if (elem) {
    HECMW_free(elem->id);
    HECMW_free(elem);
  }
  return -1;
}

/*------------------------------------------------------------------------------
  ReadFunc
*/

static int read_input(int msgno_invalid_token) {
  int token;
  char *p;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, msgno_invalid_token, "'=' required after INPUT");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_FILENAME && token != HECMW_ABLEX_NAME) {
    set_err_token(token, msgno_invalid_token, "Invalid filename for INPUT");
    return -1;
  }
  p = HECMW_ablex_get_text();
  if (strlen(p) > HECMW_FILENAME_LEN) {
    set_err(HECMW_IO_E0002, "");
    return -1;
  }
  if (HECMW_is_absolute_path(p)) {
    strcpy(include_filename, p);
  } else {
    char separator[10];
    char *dname = HECMW_dirname(grid_filename);
    sprintf(separator, "%c", HECMW_get_path_separator());
    if (strlen(dname) + strlen(separator) + strlen(p) > HECMW_FILENAME_LEN) {
      set_err(HECMW_IO_E0002, "");
      return -1;
    }
    sprintf(include_filename, "%s%s%s", dname, separator, p);
  }
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_amplitude_keyword(void) {
  int token;

  /* *AMPLITUDE */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_H_AMPLITUDE) {
    set_err_token(token, HECMW_IO_ABAQUS_E0100, "*AMPLITUDE required");
    return -1;
  }

  token = HECMW_ablex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_ABAQUS_E0100,
                  "',' required after *AMPLITUDE");
    return -1;
  }
  return 0;
}

static int read_amplitude_param_name(char *name) {
  int token;
  char *p;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E0100, "'=' required after NAME");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_NAME) {
    set_err_token(token, HECMW_IO_ABAQUS_E0100,
                  "NAME must begin with a letter");
    return -1;
  }
  p = HECMW_ablex_get_text();
  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }
  strcpy(name, p);
  HECMW_toupper(name);
  if (HECMW_io_is_reserved_name(name)) {
    set_err(HECMW_IO_E0003, "");
    return -1;
  }
  return 0;
}

static int read_amplitude_param_definition(int *definition) {
  int token;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E0100,
                  "'=' required after DEFINITION");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_K_TABULAR) {
    set_err_token(token, HECMW_IO_ABAQUS_E0100, "Invalid DEFINITION");
    return -1;
  }
  *definition = HECMW_AMP_TYPEDEF_TABULAR;
  return 0;
}

static int read_amplitude_param_time(int *time) {
  int token;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E0100, "'=' after TIME required");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_K_STEP_TIME) {
    set_err_token(token, HECMW_IO_ABAQUS_E0100, "Invalid TIME");
    return -1;
  }
  *time = HECMW_AMP_TYPETIME_STEP;
  return 0;
}

static int read_amplitude_param_value(int *value) {
  int token;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E0100, "'=' required after VALUE");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token == HECMW_ABLEX_K_RELATIVE) {
    *value = HECMW_AMP_TYPEVAL_RELATIVE;
  } else if (token == HECMW_ABLEX_K_ABSOLUTE) {
    *value = HECMW_AMP_TYPEVAL_ABSOLUTE;
  } else {
    set_err_token(token, HECMW_IO_ABAQUS_E0100, "Invalid VALUE");
    return -1;
  }
  return 0;
}

static int read_amplitude_data(char *name, int definition, int time,
                               int value) {
  int i, token;
  enum { NITEM = 4 };

  i = 0;
  while (1) {
    double val, t;

    token = HECMW_ablex_next_token();
    if (i != 0 && token == HECMW_ABLEX_NL) break;
    /* T */
    if (token != HECMW_ABLEX_DOUBLE && token != HECMW_ABLEX_INT &&
        token != ',') {
      set_err_token(token, HECMW_IO_ABAQUS_E0100, "T required");
      return -1;
    }
    if (token == ',') {
      t = 0.0;
      HECMW_ablex_unput_token();
    } else {
      t = HECMW_ablex_get_number();
    }

    /* ',' */
    token = HECMW_ablex_next_token();
    if (token != ',') {
      set_err_token(token, HECMW_IO_ABAQUS_E0100, "',' required after T");
      return -1;
    }

    /* VAL */
    token = HECMW_ablex_next_token();
    if (token != HECMW_ABLEX_DOUBLE && token != HECMW_ABLEX_INT &&
        token != ',' && token != HECMW_ABLEX_NL) {
      set_err_token(token, HECMW_IO_ABAQUS_E0100, "VAL required");
      return -1;
    }
    if (token == ',' || token == HECMW_ABLEX_NL) {
      val = 0.0;
      HECMW_ablex_unput_token();
    } else {
      val = HECMW_ablex_get_number();
    }

    /* add */
    if (HECMW_io_add_amp(name, definition, time, value, val, t) == NULL)
      return -1;

    i++;

    /* ',' or NL */
    token = HECMW_ablex_next_token();
    if (token != ',' && token != HECMW_ABLEX_NL) {
      set_err_token(token, HECMW_IO_ABAQUS_E0100, "',' or NL required");
      return -1;
    }
    if (token == ',' && i == NITEM) {
      token = HECMW_ablex_next_token();
      if (token != HECMW_ABLEX_NL) {
        set_err_token(token, HECMW_IO_ABAQUS_E0100,
                      "Only %d items allow per line", NITEM);
        return -1;
      }
      break;
    }
    if (token == HECMW_ABLEX_NL) break;
  }

  return 0;
}

static int read_amplitude(void) {
  int token, state;
  int definition                = HECMW_AMP_TYPEDEF_TABULAR;
  int time                      = HECMW_AMP_TYPETIME_STEP;
  int value                     = HECMW_AMP_TYPEVAL_RELATIVE;
  int flag_name                 = 0; /* flag for NAME */
  int flag_definition           = 0; /* flag for DEFINITION */
  int flag_time                 = 0; /* flag for TIME */
  int flag_value                = 0; /* flag for VALUE */
  int flag_input                = 0; /* flag for INPUT */
  char name[HECMW_NAME_LEN + 1] = "";
  enum {
    ST_FINISHED,
    ST_KEYWORD_LINE,
    ST_KEYWORD_LINE_PARAM,
    ST_DATA_INCLUDE,
    ST_DATA_LINE
  };

  state = ST_KEYWORD_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_KEYWORD_LINE) {
      if (read_amplitude_keyword()) return -1;
      state = ST_KEYWORD_LINE_PARAM;
    } else if (state == ST_KEYWORD_LINE_PARAM) {
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_K_NAME) {
        /* must */
        if (read_amplitude_param_name(name)) return -1;
        flag_name = 1;
      } else if (token == HECMW_ABLEX_K_DEFINITION) {
        /* optional */
        if (read_amplitude_param_definition(&definition)) return -1;
        flag_definition = 1;
      } else if (token == HECMW_ABLEX_K_TIME) {
        /* optional */
        if (read_amplitude_param_time(&time)) return -1;
        flag_time = 1;
      } else if (token == HECMW_ABLEX_K_VALUE) {
        /* optional */
        if (read_amplitude_param_value(&value)) return -1;
        flag_value = 1;
      } else if (token == HECMW_ABLEX_K_INPUT) {
        /* optional */
        if (read_input(HECMW_IO_ABAQUS_E0100)) return -1;
        flag_input = 1;
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E0100, "Unknown parameter");
        return -1;
      }

      /* check next state */
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_NL) {
        /* check NAME */
        if (!flag_name) {
          set_err(HECMW_IO_ABAQUS_E0101, "");
          return -1;
        }
        state = flag_input ? ST_DATA_INCLUDE : ST_DATA_LINE;
      } else if (token == ',') {
        ; /* continue this state */
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E0100, "Unknown parameter");
        return -1;
      }
    } else if (state == ST_DATA_INCLUDE) {
      HECMW_assert(flag_input);
      HECMW_assert(flag_name);
      if (HECMW_ablex_switch_to_include(include_filename)) return -1;
      state = ST_DATA_LINE;
    } else if (state == ST_DATA_LINE) {
      HECMW_assert(flag_name);

      if (read_amplitude_data(name, definition, time, value)) return -1;

      token = HECMW_ablex_next_token();
      if (token != HECMW_ABLEX_DOUBLE && token != HECMW_ABLEX_INT &&
          token != ',') {
        state = ST_FINISHED;
      }
      HECMW_ablex_unput_token();
    } else {
      HECMW_assert(0);
    }
  }
  return 0;
}

/*----------------------------------------------------------------------------*/
#if 0
static int
read_ecopy(void)
{
	fprintf(stderr, "*ECOPY has not implemented yet\n");
	HECMW_abort(HECMW_comm_get_comm());
	return 0;
}


/*----------------------------------------------------------------------------*/

static int
read_egen(void)
{
	fprintf(stderr, "*EGEN has not implemented yet\n");
	HECMW_abort(HECMW_comm_get_comm());
	return 0;
}
#endif

/*----------------------------------------------------------------------------*/

static int read_elset_keyword(void) {
  int token;

  /* *ELSET */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_H_ELSET) {
    set_err_token(token, HECMW_IO_ABAQUS_E0500, "*ELSET required");
    return -1;
  }

  token = HECMW_ablex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_ABAQUS_E0500, "',' required after *ELSET");
    return -1;
  }
  return 0;
}

static int read_elset_param_elset(char *elset, int *isAll) {
  int token;
  char *p;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E0500, "'=' required after ELSET");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_NAME) {
    set_err_token(token, HECMW_IO_ABAQUS_E0500,
                  "ELSET must begin with a letter");
    return -1;
  }
  p = HECMW_ablex_get_text();
  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }
  strcpy(elset, p);
  HECMW_toupper(elset);
  if (HECMW_io_is_reserved_name(elset)) {
    set_err(HECMW_IO_E0003, "");
    return -1;
  }
  if (strcmp(elset, "ALL") == 0) {
    *isAll = 1;
    HECMW_print_msg(HECMW_LOG_WARN, HECMW_IO_W1030, "");
    // return -1;
  }
  return 0;
}

static int read_elset_data(int *nelem, int **elem_array) {
  int token, i, n, *elem;
  struct hecmw_io_id *head, *prev, *p, *q;

  n    = 0;
  prev = NULL;
  head = NULL;
  elem = NULL;
  while (1) {
    struct hecmw_io_id *id;

    token = HECMW_ablex_next_token();
    if (n != 0 && token == HECMW_ABLEX_NL) break;

    id = HECMW_malloc(sizeof(*id));
    if (id == NULL) {
      set_err(errno, "");
      goto error;
    }

    /* elemX */
    if (token != HECMW_ABLEX_INT) {
      set_err_token(token, HECMW_IO_ABAQUS_E0500, "Element ID required");
      goto error;
    }
    id->id   = HECMW_ablex_get_number();
    id->next = NULL;
    if (head == NULL) {
      head = id;
    } else {
      prev->next = id;
    }
    prev = id;
    n++;

    /* ',' or NL */
    token = HECMW_ablex_next_token();
    if (token != ',' && token != HECMW_ABLEX_NL) {
      set_err_token(token, HECMW_IO_ABAQUS_E0500,
                    "',' or NL required after element ID");
      goto error;
    }
    if (token == HECMW_ABLEX_NL) break;
  }
  HECMW_assert(head);
  HECMW_assert(n > 0);

  /* add elem to group */
  elem = HECMW_malloc(sizeof(*elem) * n);
  if (elem == NULL) {
    set_err(errno, "");
    goto error;
  }
  i = 0;
  for (p = head; p; p = q) {
    q         = p->next;
    elem[i++] = p->id;
    HECMW_free(p);
  }
  head = NULL;

  *nelem      = n;
  *elem_array = elem;
  return 0;
error:
  for (p = head; p; p = q) {
    q = p->next;
    HECMW_free(p);
  }
  HECMW_free(elem);
  return -1;
}

static int read_elset_data_generate(int *nelem, int **elem_array) {
  int token, i, n, id, *elem;
  int elem1, elem2, elem3;

  /* elem1 */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_INT) {
    set_err_token(token, HECMW_IO_ABAQUS_E0500, "elem1 required");
    return -1;
  }
  elem1 = HECMW_ablex_get_number();
  if (elem1 <= 0) {
    set_err(HECMW_IO_ABAQUS_E0502, "");
    return -1;
  }

  /* ',' */
  token = HECMW_ablex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_ABAQUS_E0500, "',' required after elem1");
    return -1;
  }

  /* elem2 */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_INT) {
    set_err_token(token, HECMW_IO_ABAQUS_E0500, "elem2 required");
    return -1;
  }
  elem2 = HECMW_ablex_get_number();
  if (elem2 <= 0) {
    set_err(HECMW_IO_ABAQUS_E0502, "");
    return -1;
  }

  /* ',' or NL */
  token = HECMW_ablex_next_token();
  if (token == ',') {
    /* elem3 */
    token = HECMW_ablex_next_token();
    if (token != HECMW_ABLEX_INT) {
      set_err_token(token, HECMW_IO_ABAQUS_E0500, "Increment required");
      return -1;
    }
    elem3 = HECMW_ablex_get_number();
    if (elem3 <= 0) {
      set_err(HECMW_IO_ABAQUS_E0502, "");
      return -1;
    }

    /* NL */
    token = HECMW_ablex_next_token();
    if (token != HECMW_ABLEX_NL) {
      set_err_token(token, HECMW_IO_ABAQUS_E0500,
                    "NL required after increment");
      return -1;
    }
  } else if (token == HECMW_ABLEX_NL) {
    elem3 = 1;
  } else {
    set_err_token(token, HECMW_IO_ABAQUS_E0500,
                  "',' or NL required after elem2");
    return -1;
  }
  HECMW_assert(token == HECMW_ABLEX_NL);

  /* make element */
  if (elem1 > elem2) {
    set_err(HECMW_IO_ABAQUS_E0503,
            "Cannot generate between %d and %d with an increment of %d", elem1,
            elem2, elem3);
    return -1;
  }
  if ((elem2 - elem1) % elem3) {
    set_err(HECMW_IO_ABAQUS_E0503,
            "Cannot generate between %d and %d with an increment of %d", elem1,
            elem2, elem3);
    return -1;
  }

  n    = (elem2 - elem1) / elem3 + 1;
  elem = HECMW_malloc(sizeof(*elem) * n);
  if (elem == NULL) {
    set_err(errno, "");
    return -1;
  }

  i = 0;
  for (id = elem1; id <= elem2; id += elem3) {
    elem[i++] = id;
  }
  HECMW_assert(i == n);

  *nelem      = n;
  *elem_array = elem;
  return 0;
}

static int read_elset(void) {
  int token, state;
  int flag_elset                 = 0; /* flag for ELSET */
  int flag_generate              = 0; /* flag for GENERATE */
  int isAll                      = 0;
  char elset[HECMW_NAME_LEN + 1] = "";
  enum {
    ST_FINISHED,
    ST_KEYWORD_LINE,
    ST_KEYWORD_LINE_PARAM,
    ST_DATA_LINE,
    ST_DATA_LINE_GENERATE
  };

  state = ST_KEYWORD_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_KEYWORD_LINE) {
      if (read_elset_keyword()) return -1;
      state = ST_KEYWORD_LINE_PARAM;
    } else if (state == ST_KEYWORD_LINE_PARAM) {
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_K_ELSET) {
        /* must */
        if (read_elset_param_elset(elset, &isAll)) return -1;
        flag_elset = 1;
      } else if (token == HECMW_ABLEX_K_GENERATE) {
        /* oprtional */
        flag_generate = 1;
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E0500, "Unknown parameter");
        return -1;
      }

      /* check next parameter */
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_NL) {
        /* check */
        if (!flag_elset) {
          set_err(HECMW_IO_ABAQUS_E0501, "");
          return -1;
        }
        state = flag_generate ? ST_DATA_LINE_GENERATE : ST_DATA_LINE;
      } else if (token == ',') {
        ; /* continue this state */
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E0500, "Unknown parameter");
        return -1;
      }
    } else if (state == ST_DATA_LINE) {
      int n, *elem;

      HECMW_assert(!flag_generate);
      HECMW_assert(flag_elset);

      if (read_elset_data(&n, &elem)) return -1;

      if (HECMW_io_add_egrp(elset, n, elem) < 0) return -1;
      HECMW_free(elem);

      /* check next state */
      token = HECMW_ablex_next_token();
      if (token != HECMW_ABLEX_INT) {
        state = ST_FINISHED;
      } else {
        state = ST_DATA_LINE;
      }
      HECMW_ablex_unput_token();
    } else if (state == ST_DATA_LINE_GENERATE) {
      int n, *elem;

      HECMW_assert(flag_generate);
      HECMW_assert(flag_elset);

      if (read_elset_data_generate(&n, &elem)) return -1;

      if (HECMW_io_add_egrp(elset, n, elem) < 0) return -1;
      HECMW_free(elem);

      /* check next state */
      token = HECMW_ablex_next_token();
      if (token != HECMW_ABLEX_INT) {
        state = ST_FINISHED;
      }
      HECMW_ablex_unput_token();
    } else {
      HECMW_assert(0);
    }
  }
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_element_header(void) {
  int token;

  /* *ELEMENT */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_H_ELEMENT) {
    set_err_token(token, HECMW_IO_ABAQUS_E0600, "*ELEMENT required");
    return -1;
  }

  token = HECMW_ablex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_ABAQUS_E0600, "',' required after *ELEMENT");
    return -1;
  }
  return 0;
}

static int read_element_param_type(int *hecmw_etype, int *abaqus_etype) {
  int token;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E0600, "'=' required after TYPE");
    return -1;
  }
  token         = HECMW_ablex_next_token();
  *abaqus_etype = token;
  *hecmw_etype  = get_HECMW_etype(*abaqus_etype);
  if (*hecmw_etype == -1) {
    set_err(HECMW_IO_ABAQUS_E0601, "Invalid type: %s", HECMW_ablex_get_text());
    return -1;
  }
  if (HECMW_get_max_node(*hecmw_etype) == -1) {
    set_err(HECMW_IO_ABAQUS_E0601, "Invalid type: %s", HECMW_ablex_get_text());
    return -1;
  }
  return 0;
}

static int read_element_param_elset(char *elset) {
  char *p;
  int token;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E0600, "'=' required after ELSET");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_NAME) {
    set_err_token(token, HECMW_IO_ABAQUS_E0600,
                  "ELSET must begin with a letter");
    return -1;
  }
  p = HECMW_ablex_get_text();
  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }
  strcpy(elset, p);
  HECMW_toupper(elset);
  if (HECMW_io_is_reserved_name(elset)) {
    set_err(HECMW_IO_E0003, "");
    return -1;
  }
  if (strcmp(elset, "ALL") == 0) {
    HECMW_print_msg(HECMW_LOG_WARN, HECMW_IO_W1030, "");
    strcpy(elset, "ABAQUS_ESET_ALL");
    // return -1;
  }
  return 0;
}

static int read_element_data(int *id, int nnode, int *node) {
  int token, i;

  /* element ID */
  *id   = 0;
  token = HECMW_ablex_next_token();
  if (token == ',') {
    HECMW_ablex_unput_token();
  } else if (token == HECMW_ABLEX_INT) {
    *id = HECMW_ablex_get_number();
  } else {
    set_err_token(token, HECMW_IO_ABAQUS_E0600, "");
    return -1;
  }
  if (*id <= 0) {
    set_err_token(token, HECMW_IO_ABAQUS_E0603, "");
    return -1;
  }

  /* ',' */
  token = HECMW_ablex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_ABAQUS_E0600, "',' required after element ID");
    return -1;
  }

  /* connectivity */
  i = 0;
  while (1) {
    token = HECMW_ablex_next_token();
    if (i != 0 && token == HECMW_ABLEX_NL) continue;
    node[i] = 0;
    if (token == ',') {
      HECMW_ablex_unput_token();
    } else if (token == HECMW_ABLEX_INT) {
      node[i] = HECMW_ablex_get_number();
    } else {
      set_err(HECMW_IO_ABAQUS_E0600, "");
      return -1;
    }
    if (node[i] <= 0) {
      set_err(HECMW_IO_ABAQUS_E0604, "");
      return -1;
    }

    if (i == nnode - 1) break;

    /* ',' or NL */
    token = HECMW_ablex_next_token();
    if (token != ',' && token != HECMW_ABLEX_NL) {
      set_err_token(token, HECMW_IO_ABAQUS_E0600,
                    "',' or NL required after connectivity");
      return -1;
    }

    i++;
  }

  /* ',' */
  token = HECMW_ablex_next_token();
  if (token != ',') HECMW_ablex_unput_token();
  /* NL */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_NL) {
    set_err_token(token, HECMW_IO_ABAQUS_E0600, "NL required");
    return -1;
  }

  return 0;
}

static int read_element(void) {
  int token, state;
  int id;
  int nnode = 0;
  int node[HECMW_MAX_NODE_MAX];
  int type        = -1; /* ABAQUS element type by LEX parameter */
  int hecmw_etype = -1; /* HEC-MW element type */
  int flag_type   = 0;  /* flag for TYPE */
  int flag_elset  = 0;  /* flag for ELSET */
  int flag_input  = 0;  /* flag for INPUT */
  char elset[HECMW_NAME_LEN + 1] = "";
  enum {
    ST_FINISHED,
    ST_KEYWORD_LINE,
    ST_KEYWORD_LINE_PARAM,
    ST_DATA_INCLUDE,
    ST_DATA_LINE,
    ST_DATA_LINE_REGIST
  };

  state = ST_KEYWORD_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_KEYWORD_LINE) {
      if (read_element_header()) return -1;
      state = ST_KEYWORD_LINE_PARAM;
    } else if (state == ST_KEYWORD_LINE_PARAM) {
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_K_TYPE) {
        /* must */
        if (read_element_param_type(&hecmw_etype, &type)) return -1;
        flag_type = 1;

        /* get # of connectivity */
        nnode = HECMW_get_max_node(hecmw_etype);
        HECMW_assert(nnode > 0);
        HECMW_assert(nnode <= HECMW_MAX_NODE_MAX);
      } else if (token == HECMW_ABLEX_K_ELSET) {
        /* optional */
        if (read_element_param_elset(elset)) return -1;
        flag_elset = 1;
      } else if (token == HECMW_ABLEX_K_INPUT) {
        /* optional */
        if (read_input(HECMW_IO_ABAQUS_E0600)) return -1;
        flag_input = 1;
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E0600, "Unknown parameter");
        return -1;
      }

      /* check next state */
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_NL) {
        /* check TYPE */
        if (!flag_type) {
          set_err(HECMW_IO_ABAQUS_E0606, "");
          return -1;
        }
        state = flag_input ? ST_DATA_INCLUDE : ST_DATA_LINE;
      } else if (token == ',') {
        ; /* continue this state */
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E0600, "Unknown parameter");
        return -1;
      }
    } else if (state == ST_DATA_INCLUDE) {
      HECMW_assert(flag_input);
      HECMW_assert(flag_type);
      if (HECMW_ablex_switch_to_include(include_filename)) return -1;
      state = ST_DATA_LINE;
    } else if (state == ST_DATA_LINE) {
      HECMW_assert(flag_type);
      if (read_element_data(&id, nnode, node)) return -1;
      state = ST_DATA_LINE_REGIST;
    } else if (state == ST_DATA_LINE_REGIST) {
      Integer *secopt;

      HECMW_assert(flag_type);

      /* convert connectivity */
      if (HECMW_convert_connectivity(HECMW_CONNTYPE_ABAQUS, hecmw_etype, node))
        return -1;

      /* add element */
      if (HECMW_io_add_elem(id, hecmw_etype, node, 0, NULL) == NULL) return -1;

      /* save secopt */
      secopt = HECMW_malloc(sizeof(*secopt));
      if (secopt == NULL) {
        set_err(errno, "");
        return -1;
      }
      secopt->i = get_secopt_abaqus(type);
      HECMW_assert(secopt->i != -1);
      if (elem_secopt == NULL) {
        elem_secopt =
            (struct hecmw_map_int *)HECMW_malloc(sizeof(struct hecmw_map_int));
        if (elem_secopt == NULL) return -1;
        if (HECMW_map_int_init(elem_secopt, free_Integer)) return -1;
      }
      if (HECMW_map_int_add(elem_secopt, id, secopt) < 0) return -1;

      /* add element to egroup */
      if (HECMW_io_add_egrp("ALL", 1, &id) < 0) return -1;

      if (flag_elset) {
        if (HECMW_io_add_egrp(elset, 1, &id) < 0) return -1;
      }

      /* check next state */
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_INT) {
        state = ST_DATA_LINE;
      } else {
        state = ST_FINISHED;
      }
      HECMW_ablex_unput_token();
    } else {
      HECMW_assert(0);
    }
  }
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_equation_keyword(int *token) {
  /* *EQUATION */
  *token = HECMW_ablex_next_token();
  if (*token != HECMW_ABLEX_H_EQUATION) {
    set_err_token(*token, HECMW_IO_ABAQUS_E0700, "*EQUATION required");
    return -1;
  }
  *token = HECMW_ablex_next_token();
  if (*token != ',' && *token != HECMW_ABLEX_NL) {
    set_err_token(*token, HECMW_IO_ABAQUS_E0700,
                  "',' or NL required after *EQUATION");
    return -1;
  }
  return 0;
}

static int read_equation_data_line1(int *neq) {
  int token;

  /* NEQ */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_INT) {
    set_err_token(token, HECMW_IO_ABAQUS_E0700, "required NEQ");
    return -1;
  }
  *neq = HECMW_ablex_get_number();
  if (*neq < 2) {
    set_err(HECMW_IO_ABAQUS_E0701, "");
    return -1;
  }

  /* NL */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_NL) {
    set_err_token(token, HECMW_IO_ABAQUS_E0700, "NL required after NEQ");
    return -1;
  }
  return 0;
}

static int read_equation_data_line2(int neq) {
  int i, token;
  int is_node = 0;
  int is_ngrp = 0;
  enum { NITEM = 4 };
  struct hecmw_io_mpcitem *mpcitem = NULL;

  mpcitem = HECMW_malloc(sizeof(*mpcitem) * neq);
  if (mpcitem == NULL) {
    set_err(errno, "");
    goto error;
  }

  for (i = 0; i < neq; i++) {
    token = HECMW_ablex_next_token();
    if (i != 0 && token == HECMW_ABLEX_NL) break;

    /* nod */
    if (token == HECMW_ABLEX_INT) {
      if (is_ngrp) {
        set_err(HECMW_IO_ABAQUS_E0702, "");
        goto error;
      }
      mpcitem[i].node = HECMW_ablex_get_number();
      strcpy(mpcitem[i].ngrp, "");
      is_node = 1;
    } else if (token == HECMW_ABLEX_NAME) {
      char *p = HECMW_ablex_get_text();
      if (is_node) {
        set_err(HECMW_IO_ABAQUS_E0702, "");
        goto error;
      }
      if (strlen(p) > HECMW_NAME_LEN) {
        set_err(HECMW_IO_E0001, "");
        goto error;
      }
      strcpy(mpcitem[i].ngrp, p);
      HECMW_toupper(mpcitem[i].ngrp);
      if (HECMW_io_is_reserved_name(mpcitem[i].ngrp)) {
        set_err(HECMW_IO_E0003, "");
        goto error;
      }
      mpcitem[i].node = -1;
      is_ngrp         = 1;
    } else {
      set_err_token(token, HECMW_IO_ABAQUS_E0700, "Node ID or NGRP required");
      goto error;
    }

    /* ',' */
    token = HECMW_ablex_next_token();
    if (token != ',') {
      set_err_token(token, HECMW_IO_ABAQUS_E0700, "',' required after node");
      goto error;
    }

    /* DOF */
    token = HECMW_ablex_next_token();
    if (token != HECMW_ABLEX_INT) {
      set_err(HECMW_IO_ABAQUS_E0703, "");
      goto error;
    }
    mpcitem[i].dof = HECMW_ablex_get_number();
    if (HECMW_io_check_mpc_dof(mpcitem[i].dof)) {
      set_err(HECMW_IO_ABAQUS_E0703, "");
      goto error;
    }

    /* ',' */
    token = HECMW_ablex_next_token();
    if (token != ',') {
      set_err_token(token, HECMW_IO_ABAQUS_E0700, "',' required after DOF");
      goto error;
    }

    /* A */
    token = HECMW_ablex_next_token();
    if (token != HECMW_ABLEX_DOUBLE && token != HECMW_ABLEX_INT) {
      set_err_token(token, HECMW_IO_ABAQUS_E0700, "A(coefficient) required ");
      goto error;
    }
    mpcitem[i].a = HECMW_ablex_get_number();

    /* ',' or NL */
    token = HECMW_ablex_next_token();
    if (token != ',' && token != HECMW_ABLEX_NL) {
      set_err_token(token, HECMW_IO_ABAQUS_E0700,
                    "',' or NL required after coefficient");
      goto error;
    }
    if (token == ',' && i == NITEM - 1) {
      token = HECMW_ablex_next_token();
      if (token != HECMW_ABLEX_NL) {
        set_err_token(token, HECMW_IO_ABAQUS_E0700, "NL required");
        goto error;
      }
      continue;
    }
    if (token == HECMW_ABLEX_NL) continue;
  }

  /* add */
  if (HECMW_io_add_mpc(neq, mpcitem, 0.0) == NULL) goto error;
  HECMW_free(mpcitem);
  return 0;
error:
  HECMW_free(mpcitem);
  return -1;
}

static int read_equation(void) {
  int token, state;
  int neq        = -1;
  int flag_input = 0; /* flag for INPUT */
  enum {
    ST_FINISHED,
    ST_KEYWORD_LINE,
    ST_KEYWORD_LINE_PARAM,
    ST_DATA_INCLUDE,
    ST_DATA_LINE1,
    ST_DATA_LINE2
  };

  state = ST_KEYWORD_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_KEYWORD_LINE) {
      if (read_equation_keyword(&token)) return -1;
      if (token == ',') {
        state = ST_KEYWORD_LINE_PARAM;
      } else if (token == HECMW_ABLEX_NL) {
        state = ST_DATA_LINE1;
      } else {
        HECMW_assert(0);
      }
    } else if (state == ST_KEYWORD_LINE_PARAM) {
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_K_INPUT) {
        /* optional */
        if (read_input(HECMW_IO_ABAQUS_E0700)) return -1;
        flag_input = 1;
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E0700, "Unknown parameter");
        return -1;
      }

      /* check next state */
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_NL) {
        state = flag_input ? ST_DATA_INCLUDE : ST_DATA_LINE1;
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E0700, "NL required");
        return -1;
      }
    } else if (state == ST_DATA_INCLUDE) {
      HECMW_assert(flag_input);
      if (HECMW_ablex_switch_to_include(include_filename)) return -1;
      state = ST_DATA_LINE1;
    } else if (state == ST_DATA_LINE1) {
      if (read_equation_data_line1(&neq)) return -1;
      /* set next state */
      state = ST_DATA_LINE2;
    } else if (state == ST_DATA_LINE2) {
      HECMW_assert(neq != -1);
      if (read_equation_data_line2(neq)) return -1;
      /* check next state */
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_INT) {
        state = ST_DATA_LINE1;
      } else {
        state = ST_FINISHED;
      }
      HECMW_ablex_unput_token();
    } else {
      HECMW_assert(0);
    }
  }
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_heading(void) {
  int token, len;
  char *p;
  struct hecmw_io_header *header;

  header = HECMW_malloc(sizeof(struct hecmw_io_header));
  if (header == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* *HEADING */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_H_HEADING) {
    set_err_token(token, HECMW_IO_ABAQUS_E0800, "*HEADING required");
    return -1;
  }

  /* get header data */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_HEADER) {
    set_err_token(token, HECMW_IO_ABAQUS_E0800,
                  "TITLE required after *HEADING");
    return -1;
    /* set_err_token(token, HECMW_IO_ABAQUS_E0800, "TITLE ignored after
    *HEADING");
    header->header[0] = ' ';
    header->header[1] = '\0';*/
  } else {
    p = HECMW_ablex_get_text();
    while (*p && *p == ' ') p++;
    if (p == NULL) p                = "";
    len                             = strlen(p);
    if (len > HECMW_HEADER_LEN) len = HECMW_HEADER_LEN;
    strncpy(header->header, p, len);
    header->header[len] = '\0';
  }

  /* Note: * NL is ignored by LEX until the end of the header data. */

  /* Ignore the rest of the header data */
  while (HECMW_ablex_next_token() == HECMW_ABLEX_HEADER)
    ;
  HECMW_ablex_unput_token();

  /* set */
  HECMW_io_set_header(header);

  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_include(void) {
  int token;

  /* !INCLUDE */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_H_INCLUDE) {
    set_err_token(token, HECMW_IO_ABAQUS_E0900, "*INCLUDE required");
    return -1;
  }

  /* ',' */
  token = HECMW_ablex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_ABAQUS_E0900, "',' required after *INCLUDE");
    return -1;
  }

  /* INPUT */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_K_INPUT) {
    set_err_token(token, HECMW_IO_ABAQUS_E0901, "");
    return -1;
  }

  /* =filename */
  if (read_input(HECMW_IO_ABAQUS_E0900)) return -1;

  /* NL */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_NL) {
    set_err_token(token, HECMW_IO_ABAQUS_E0900,
                  "NL required after INPUT value");
    return -1;
  }

  /* include */
  if (HECMW_ablex_switch_to_include(include_filename)) return -1;

  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_initial_keyword(void) {
  int token;

  /* *INITIAL CONDITIONS */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_H_INITIAL) {
    set_err_token(token, HECMW_IO_ABAQUS_E1000, "*INITIAL CONDITIONS required");
    return -1;
  }

  token = HECMW_ablex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_ABAQUS_E1001, "");
    return -1;
  }

  return 0;
}

static int read_initial_param_type(int *type) {
  int token;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E1000, "'=' required after TYPE");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_K_TEMPERATURE) {
    set_err_token(token, HECMW_IO_ABAQUS_E1000, "TEMPERATURE required");
    return -1;
  }
  *type = HECMW_INITIAL_TYPE_TEMPERATURE;
  return 0;
}

static int read_initial_data(int type) {
  int node, token;
  char *ngrp;
  double val;

  /* node or ngrp */
  node  = -1;
  ngrp  = NULL;
  token = HECMW_ablex_next_token();
  if (token == ',') {
    set_err(HECMW_IO_ABAQUS_E1002, "");
    return -1;
  } else if (token == HECMW_ABLEX_INT) {
    node = HECMW_ablex_get_number();
    if (node <= 0) {
      set_err(HECMW_IO_ABAQUS_E1002, "");
      return -1;
    }
  } else if (token == HECMW_ABLEX_NAME) {
    ngrp = HECMW_ablex_get_text();
    if (strlen(ngrp) > HECMW_NAME_LEN) {
      set_err(HECMW_IO_E0001, "");
      return -1;
    }
    HECMW_toupper(ngrp);
    ngrp = HECMW_strdup(ngrp);
    if (ngrp == NULL) {
      set_err(errno, "");
      return -1;
    }
  } else {
    set_err_token(token, HECMW_IO_ABAQUS_E1000,
                  "Node ID or NGROUP name required");
    return -1;
  }

  /* ',' */
  token = HECMW_ablex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_ABAQUS_E1000, "',' required after node");
    return -1;
  }

  /* VAL */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_DOUBLE && token != HECMW_ABLEX_INT &&
      token != HECMW_ABLEX_NL) {
    set_err_token(token, HECMW_IO_ABAQUS_E1000, "VAL required");
    return -1;
  }
  if (token == HECMW_ABLEX_NL) {
    val = 0.0;
    HECMW_ablex_unput_token();
  } else {
    val = HECMW_ablex_get_number();
  }

  /* skip this line */
  while ((token = HECMW_ablex_next_token()) && token != HECMW_ABLEX_NL)
    ;
  if (token != HECMW_ABLEX_NL) {
    set_err_token(token, HECMW_IO_ABAQUS_E1000, "NL required");
    return -1;
  }

  /* add */
  HECMW_assert(type != -1);
  if (HECMW_io_add_initial(type, node, ngrp, val) == NULL) return -1;
  HECMW_free(ngrp);

  return 0;
}

static int read_initial(void) {
  int token, state;
  int type       = -1;
  int flag_type  = 0; /* flag for TYPE */
  int flag_input = 0; /* flag for INPUT */
  enum {
    ST_FINISHED,
    ST_KEYWORD_LINE,
    ST_KEYWORD_LINE_PARAM,
    ST_DATA_INCLUDE,
    ST_DATA_LINE
  };

  state = ST_KEYWORD_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_KEYWORD_LINE) {
      if (read_initial_keyword()) return -1;
      state = ST_KEYWORD_LINE_PARAM;
    } else if (state == ST_KEYWORD_LINE_PARAM) {
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_K_TYPE) {
        /* must */
        if (read_initial_param_type(&type)) return -1;
        flag_type = 1;
      } else if (token == HECMW_ABLEX_K_INPUT) {
        /* oprtional */
        if (read_input(HECMW_IO_ABAQUS_E1000)) return -1;
        flag_input = 1;
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E1000, "Unknown parameter");
        return -1;
      }

      /* check next parameter */
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_NL) {
        if (flag_input) {
          state = ST_DATA_INCLUDE;
        } else {
          state = ST_DATA_LINE;
        }
        /* check */
        if (!flag_type) {
          set_err(HECMW_IO_ABAQUS_E1001, "");
          return -1;
        }
      } else if (token == ',') {
        ; /* continue this state */
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E1000, "Unknown parameter");
        return -1;
      }
    } else if (state == ST_DATA_INCLUDE) {
      HECMW_assert(flag_input);
      HECMW_assert(flag_type);
      if (HECMW_ablex_switch_to_include(include_filename)) return -1;
      state = ST_DATA_LINE;
    } else if (state == ST_DATA_LINE) {
      HECMW_assert(flag_type);
      if (read_initial_data(type)) return -1;
      /* check next state */
      token = HECMW_ablex_next_token();
      if (token != HECMW_ABLEX_INT && token != HECMW_ABLEX_NAME) {
        state = ST_FINISHED;
      }
      HECMW_ablex_unput_token();
    } else {
      HECMW_assert(0);
    }
  }
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_conductivity_keyword(int *last_token) {
  int token;

  /* *CONDUCTIVITY */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_H_CONDUCTIVITY) {
    set_err_token(token, HECMW_IO_ABAQUS_E2500, "*CONDUCTIVITY required");
    return -1;
  }

  token = HECMW_ablex_next_token();
  if (token != ',' && token != HECMW_ABLEX_NL) {
    set_err_token(token, HECMW_IO_ABAQUS_E2500,
                  "',' or NL required after *CONDUCTIVITY");
    return -1;
  }
  *last_token = token;

  return 0;
}

static int read_conductivity_param_dependencies(int *dependencies) {
  return read_param_dependencies(dependencies, HECMW_IO_ABAQUS_E2500);
}

static int read_conductivity_param_type(int *type) {
  int token;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E2500, "'=' required after TYPE");
    return -1;
  }
  token = HECMW_ablex_next_token();
  switch (token) {
    case HECMW_ABLEX_K_ISOTROPIC:   /* fall through */
    case HECMW_ABLEX_K_ORTHOTROPIC: /* fall through */
    case HECMW_ABLEX_K_ANISOTROPIC:
      break;
    default:
      set_err_token(token, HECMW_IO_ABAQUS_E2500, "Invalid TYPE");
      return -1;
  }
  *type = token;

  return 0;
}

static int read_conductivity_data(int type, int dependencies,
                                  struct hecmw_io_matitem **matitem) {
  int n, temp_col;
  n        = 0;
  temp_col = 0;

  switch (type) {
    case HECMW_ABLEX_K_ISOTROPIC:
      n = temp_col = 2;
      break;
    case HECMW_ABLEX_K_ORTHOTROPIC:
      n = temp_col = 4;
      break;
    case HECMW_ABLEX_K_ANISOTROPIC:
      n = temp_col = 7;
      break;
    default:
      HECMW_assert(0);
  }
  return read_mat_data_common(n + dependencies, 8, temp_col, matitem,
                              HECMW_IO_ABAQUS_E2500);
}

static int read_conductivity(void) {
  int token, state;
  int flag_dependencies = 0; /* flag for DEPENDENCIES */
  int flag_type         = 0; /* flag for TYPE */
  int dependencies      = 0;
  int type              = HECMW_ABLEX_K_ISOTROPIC;
  enum { ST_FINISHED, ST_KEYWORD_LINE, ST_KEYWORD_LINE_PARAM, ST_DATA_LINE };

  state = ST_KEYWORD_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_KEYWORD_LINE) {
      if (read_conductivity_keyword(&token)) return -1;
      if (token == ',') {
        state = ST_KEYWORD_LINE_PARAM;
      } else if (token == HECMW_ABLEX_NL) {
        state = ST_DATA_LINE;
      } else {
        HECMW_assert(0);
      }
    } else if (state == ST_KEYWORD_LINE_PARAM) {
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_K_DEPENDENCIES) {
        if (read_conductivity_param_dependencies(&dependencies)) return -1;
        flag_dependencies = 1;
      } else if (token == HECMW_ABLEX_K_TYPE) {
        if (read_conductivity_param_type(&type)) return -1;
        flag_type = 1;
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E2500, "Unknown parameter");
        return -1;
      }
      token = HECMW_ablex_next_token();
      if (token != ',' && token != HECMW_ABLEX_NL) {
        set_err_token(token, HECMW_IO_ABAQUS_E2500, "Unknown parameter");
        return -1;
      }
      if (token == HECMW_ABLEX_NL) {
        state = ST_DATA_LINE;
      }
    } else if (state == ST_DATA_LINE) {
      struct hecmw_io_matitem *item;
      if (read_conductivity_data(type, dependencies, &item)) return -1;
      if (add_mat_data(HECMW_ABLEX_H_CONDUCTIVITY, item)) return -1;
      state = ST_FINISHED;
    } else {
      HECMW_assert(0);
    }
  }
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_density_keyword(int *last_token) {
  int token;

  /* *DENSITY */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_H_DENSITY) {
    set_err_token(token, HECMW_IO_ABAQUS_E2200, "*DENSITY required");
    return -1;
  }

  token = HECMW_ablex_next_token();
  if (token != ',' && token != HECMW_ABLEX_NL) {
    set_err_token(token, HECMW_IO_ABAQUS_E2200,
                  "',' or NL required after *DENSITY");
    return -1;
  }
  *last_token = token;

  return 0;
}

static int read_density_param_dependencies(int *dependencies) {
  return read_param_dependencies(dependencies, HECMW_IO_ABAQUS_E2200);
}

static int read_density_data(int dependencies,
                             struct hecmw_io_matitem **matitem) {
  return read_mat_data_common(2 + dependencies, 8, 2, matitem,
                              HECMW_IO_ABAQUS_E2200);
}

static int read_density(void) {
  int token, state;
  int flag_dependencies = 0; /* flag for DEPENDENCIES */
  int dependencies      = 0;
  enum { ST_FINISHED, ST_KEYWORD_LINE, ST_KEYWORD_LINE_PARAM, ST_DATA_LINE };

  state = ST_KEYWORD_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_KEYWORD_LINE) {
      if (read_density_keyword(&token)) return -1;
      if (token == ',') {
        state = ST_KEYWORD_LINE_PARAM;
      } else if (token == HECMW_ABLEX_NL) {
        state = ST_DATA_LINE;
      } else {
        HECMW_assert(0);
      }
    } else if (state == ST_KEYWORD_LINE_PARAM) {
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_K_DEPENDENCIES) {
        if (read_density_param_dependencies(&dependencies)) return -1;
        flag_dependencies = 1;
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E2200, "Unknown parameter");
        return -1;
      }
      token = HECMW_ablex_next_token();
      if (token != HECMW_ABLEX_NL) {
        set_err_token(token, HECMW_IO_ABAQUS_E2200, "NL required");
        return -1;
      }
      state = ST_DATA_LINE;
    } else if (state == ST_DATA_LINE) {
      struct hecmw_io_matitem *item;
      if (read_density_data(dependencies, &item)) return -1;
      if (add_mat_data(HECMW_ABLEX_H_DENSITY, item)) return -1;
      state = ST_FINISHED;
    } else {
      HECMW_assert(0);
    }
  }
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_elastic_keyword(int *last_token) {
  int token;

  /* *ELASTIC */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_H_ELASTIC) {
    set_err_token(token, HECMW_IO_ABAQUS_E2300, "*ELASTIC required");
    return -1;
  }

  token = HECMW_ablex_next_token();
  if (token != ',' && token != HECMW_ABLEX_NL) {
    set_err_token(token, HECMW_IO_ABAQUS_E2300,
                  "',' or NL required after *ELASTIC");
    return -1;
  }
  *last_token = token;

  return 0;
}

static int read_elastic_param_dependencies(int *dependencies) {
  return read_param_dependencies(dependencies, HECMW_IO_ABAQUS_E2300);
}

static int read_elastic_param_type(int *type) {
  int token;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E2300, "'=' required after TYPE");
    return -1;
  }
  token = HECMW_ablex_next_token();
  switch (token) {
    case HECMW_ABLEX_K_ISOTROPIC:             /* fall through */
    case HECMW_ABLEX_K_ENGINEERING_CONSTANTS: /* fall through */
    case HECMW_ABLEX_K_LAMINA:                /* fall through */
    case HECMW_ABLEX_K_ORTHOTROPIC:           /* fall through */
    case HECMW_ABLEX_K_ANISOTROPIC:
      break;
    default:
      set_err_token(token, HECMW_IO_ABAQUS_E2300, "Invalid TYPE");
      return -1;
  }
  *type = token;

  return 0;
}

static int read_elastic_param_moduli(int *moduli) {
  int token;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E2300, "'=' required after TYPE");
    return -1;
  }
  token = HECMW_ablex_next_token();
  switch (token) {
    case HECMW_ABLEX_K_INSTANTANEOUS:
      break;
    default:
      set_err_token(token, HECMW_IO_ABAQUS_E2300, "Invalid TYPE");
      return -1;
  }
  *moduli = token;

  return 0;
}

static int read_elastic_data(int type, int dependencies,
                             struct hecmw_io_matitem **matitem) {
  int n, temp_col;
  n        = 0;
  temp_col = 0;

  switch (type) {
    case HECMW_ABLEX_K_ISOTROPIC:
      n = temp_col = 3;
      break;
    case HECMW_ABLEX_K_ENGINEERING_CONSTANTS: /* fall through */
    case HECMW_ABLEX_K_ORTHOTROPIC:
      n = temp_col = 10;
      break;
    case HECMW_ABLEX_K_LAMINA:
      n = temp_col = 7;
      break;
    case HECMW_ABLEX_K_ANISOTROPIC:
      n = temp_col = 22;
      break;
    default:
      HECMW_assert(0);
  }
  return read_mat_data_common(n + dependencies, 8, temp_col, matitem,
                              HECMW_IO_ABAQUS_E2300);
}

static int read_elastic(void) {
  int token, state;
  int flag_dependencies = 0; /* flag for DEPENDENCIES */
  int flag_type         = 0; /* flag for TYPE */
  int flag_moduli       = 0; /* flag for MODULI */
  int dependencies      = 0;
  int type              = HECMW_ABLEX_K_ISOTROPIC;
  int moduli            = 0;
  enum { ST_FINISHED, ST_KEYWORD_LINE, ST_KEYWORD_LINE_PARAM, ST_DATA_LINE };

  state = ST_KEYWORD_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_KEYWORD_LINE) {
      if (read_elastic_keyword(&token)) return -1;
      if (token == ',') {
        state = ST_KEYWORD_LINE_PARAM;
      } else if (token == HECMW_ABLEX_NL) {
        state = ST_DATA_LINE;
      } else {
        HECMW_assert(0);
      }
    } else if (state == ST_KEYWORD_LINE_PARAM) {
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_K_DEPENDENCIES) {
        if (read_elastic_param_dependencies(&dependencies)) return -1;
        flag_dependencies = 1;
      } else if (token == HECMW_ABLEX_K_TYPE) {
        if (read_elastic_param_type(&type)) return -1;
        flag_type = 1;
      } else if (token == HECMW_ABLEX_K_MODULI) {
        if (read_elastic_param_moduli(&moduli)) return -1;
        flag_moduli = 1;
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E2300, "Unknown parameter");
        return -1;
      }
      token = HECMW_ablex_next_token();
      if (token != ',' && token != HECMW_ABLEX_NL) {
        set_err_token(token, HECMW_IO_ABAQUS_E2300, "Unknown parameter");
        return -1;
      }
      if (token == HECMW_ABLEX_NL) {
        state = ST_DATA_LINE;
      }
    } else if (state == ST_DATA_LINE) {
      struct hecmw_io_matitem *item;
      if (read_elastic_data(type, dependencies, &item)) return -1;
      if (add_mat_data(HECMW_ABLEX_H_ELASTIC, item)) return -1;
      state = ST_FINISHED;
    } else {
      HECMW_assert(0);
    }
  }
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_specific_keyword(int *last_token) {
  int token;

  /* *SPECIFIC HEAT */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_H_SPECIFIC_HEAT) {
    set_err_token(token, HECMW_IO_ABAQUS_E2400, "*SPECIFIC HEAT required");
    return -1;
  }

  token = HECMW_ablex_next_token();
  if (token != ',' && token != HECMW_ABLEX_NL) {
    set_err_token(token, HECMW_IO_ABAQUS_E2400,
                  "',' or NL required after *SPECIFIC HEAT");
    return -1;
  }
  *last_token = token;

  return 0;
}

static int read_specific_param_dependencies(int *dependencies) {
  return read_param_dependencies(dependencies, HECMW_IO_ABAQUS_E2400);
}

static int read_specific_data(int dependencies,
                              struct hecmw_io_matitem **matitem) {
  return read_mat_data_common(2 + dependencies, 8, 2, matitem,
                              HECMW_IO_ABAQUS_E2400);
}

static int read_specific_heat(void) {
  int token, state;
  int flag_dependencies = 0; /* flag for DEPENDENCIES */
  int dependencies      = 0;
  enum { ST_FINISHED, ST_KEYWORD_LINE, ST_KEYWORD_LINE_PARAM, ST_DATA_LINE };

  state = ST_KEYWORD_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_KEYWORD_LINE) {
      if (read_specific_keyword(&token)) return -1;
      if (token == ',') {
        state = ST_KEYWORD_LINE_PARAM;
      } else if (token == HECMW_ABLEX_NL) {
        state = ST_DATA_LINE;
      } else {
        HECMW_assert(0);
      }
    } else if (state == ST_KEYWORD_LINE_PARAM) {
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_K_DEPENDENCIES) {
        if (read_specific_param_dependencies(&dependencies)) return -1;
        flag_dependencies = 1;
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E2400, "Unknown parameter");
        return -1;
      }
      token = HECMW_ablex_next_token();
      if (token != HECMW_ABLEX_NL) {
        set_err_token(token, HECMW_IO_ABAQUS_E2400, "NL required");
        return -1;
      }
      state = ST_DATA_LINE;
    } else if (state == ST_DATA_LINE) {
      struct hecmw_io_matitem *item;
      if (read_specific_data(dependencies, &item)) return -1;
      if (add_mat_data(HECMW_ABLEX_H_SPECIFIC_HEAT, item)) return -1;
      state = ST_FINISHED;
    } else {
      HECMW_assert(0);
    }
  }
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_material_param_name(char *name) {
  int token;
  char *p;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E1100, "'=' required after NAME");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_NAME) {
    set_err_token(token, HECMW_IO_ABAQUS_E1100,
                  "NAME must begin with a letter");
    return -1;
  }
  p = HECMW_ablex_get_text();
  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }
  strcpy(name, p);
  HECMW_toupper(name);
  if (HECMW_io_is_reserved_name(name)) {
    set_err(HECMW_IO_E0003, "");
    return -1;
  }
  if (HECMW_io_get_mat(name)) {
    set_err(HECMW_IO_ABAQUS_E1102, "%s already exists", name);
    return -1;
  }
  return 0;
}

static int read_material(void) {
  int token;
  char name[HECMW_NAME_LEN + 1] = "";

  /* *MATERIAL */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_H_MATERIAL) {
    set_err_token(token, HECMW_IO_ABAQUS_E1100, "*MATERIAL required");
    return -1;
  }

  token = HECMW_ablex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_ABAQUS_E1101, "");
    return -1;
  }

  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_K_NAME) {
    set_err_token(token, HECMW_IO_ABAQUS_E1100, "Unknown parameter");
    return -1;
  }

  if (read_material_param_name(name)) return -1;

  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_NL) {
    set_err_token(token, HECMW_IO_ABAQUS_E1100, "NL required");
    return -1;
  }

  strcpy(matname, name);

  return 0;
}

/*----------------------------------------------------------------------------*/
#if 0
static int
read_ncopy(void)
{
	fprintf(stderr, "*NCOPY has not implemented yet\n");
	HECMW_abort(HECMW_comm_get_comm());
	return 0;
}


/*----------------------------------------------------------------------------*/

static int
read_nfill(void)
{
	fprintf(stderr, "*NFILL has not implemented yet\n");
	HECMW_abort(HECMW_comm_get_comm());
	return 0;
}


/*----------------------------------------------------------------------------*/

static int
read_ngen(void)
{
	fprintf(stderr, "*NGEN has not implemented yet\n");
	HECMW_abort(HECMW_comm_get_comm());
	return 0;
}
#endif

/*----------------------------------------------------------------------------*/

static int read_nset_keyword(void) {
  int token;

  /* *NSET */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_H_NSET) {
    set_err_token(token, HECMW_IO_ABAQUS_E1500, "*NSET required");
    return -1;
  }

  token = HECMW_ablex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_ABAQUS_E1500, "',' required after *NSET");
    return -1;
  }
  return 0;
}

static int read_nset_param_nset(char *nset, int *isAll) {
  int token;
  char *p;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E1500, "'=' required after NSET");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_NAME) {
    set_err_token(token, HECMW_IO_ABAQUS_E1500,
                  "NSET must begin with a letter");
    return -1;
  }
  p = HECMW_ablex_get_text();
  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }
  strcpy(nset, p);
  HECMW_toupper(nset);
  if (HECMW_io_is_reserved_name(nset)) {
    set_err(HECMW_IO_E0003, "");
    return -1;
  }
  if (strcmp(nset, "EQUATION_BLOCK") == 0) {
    set_err(HECMW_IO_E0003, "Reserved name: %s", nset);
    return -1;
  }
  if (strcmp(nset, "ALL") == 0) {
    HECMW_print_msg(HECMW_LOG_WARN, HECMW_IO_W1030, "");
    strcpy(nset, "ABAQUS_NSET_ALL");
    *isAll = 1;
    // return -1;
  }
  return 0;
}

static int read_nset_param_instance() {
  int token;
  char *p;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E1500, "'=' required after INSTANCE");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_NAME) {
    set_err_token(token, HECMW_IO_ABAQUS_E1500,
                  "NSET must begin with a letter");
    return -1;
  }
  return 0;
}

static int read_nset_data(int *nnode, int **node_array) {
  int i, n, *node, token;
  struct hecmw_io_id *head, *prev, *p, *q;

  n    = 0;
  head = NULL;
  node = NULL;
  while (1) {
    struct hecmw_io_id *id;

    token = HECMW_ablex_next_token();
    if (n != 0 && token == HECMW_ABLEX_NL) break;

    id = HECMW_malloc(sizeof(*id));
    if (id == NULL) {
      set_err(errno, "");
      goto error;
    }

    /* nodX */
    if (token != HECMW_ABLEX_INT && token != ',') {
      set_err_token(token, HECMW_IO_ABAQUS_E1500, "Node ID required");
      goto error;
    }
    if (token == ',') {
      id->id = 0;
      HECMW_ablex_unput_token();
    } else {
      id->id = HECMW_ablex_get_number();
    }
    id->next = NULL;
    if (head == NULL) {
      head = id;
    } else {
      prev->next = id;
    }
    prev = id;
    n++;

    /* ',' or NL */
    token = HECMW_ablex_next_token();
    if (token != ',' && token != HECMW_ABLEX_NL) {
      set_err_token(token, HECMW_IO_ABAQUS_E1500,
                    "',' or NL required after node ID");
      goto error;
    }
    if (token == HECMW_ABLEX_NL) break;
  }

  HECMW_assert(head);
  HECMW_assert(n > 0);

  node = HECMW_malloc(sizeof(*node) * n);
  if (node == NULL) {
    set_err(errno, "");
    goto error;
  }
  i = 0;
  for (p = head; p; p = q) {
    q         = p->next;
    node[i++] = p->id;
    HECMW_free(p);
  }
  head = NULL;

  *nnode      = n;
  *node_array = node;
  return 0;
error:
  for (p = head; p; p = q) {
    q = p->next;
    HECMW_free(p);
  }
  HECMW_free(node);
  return -1;
}

static int read_nset_data_generate(int *nnode, int **node_array) {
  int token, i, n, id, *node, nod1, nod2, nod3;

  /* nod1 */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_INT) {
    set_err_token(token, HECMW_IO_ABAQUS_E1500, "nod1 required");
    return -1;
  }
  nod1 = HECMW_ablex_get_number();
  if (nod1 <= 0) {
    set_err(HECMW_IO_ABAQUS_E1502, "");
    return -1;
  }

  /* ',' */
  token = HECMW_ablex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_ABAQUS_E1500, "',' required after nod1");
    return -1;
  }

  /* nod2 */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_INT) {
    set_err_token(token, HECMW_IO_ABAQUS_E1500, "nod2 required");
    return -1;
  }
  nod2 = HECMW_ablex_get_number();
  if (nod2 <= 0) {
    set_err(HECMW_IO_ABAQUS_E1502, "");
    return -1;
  }

  /* ',' or NL */
  token = HECMW_ablex_next_token();
  if (token == ',') {
    /* nod3 */
    token = HECMW_ablex_next_token();
    if (token != HECMW_ABLEX_INT) {
      set_err_token(token, HECMW_IO_ABAQUS_E1500, "Increment required");
      return -1;
    }
    nod3 = HECMW_ablex_get_number();
    if (nod3 <= 0) {
      set_err(HECMW_IO_ABAQUS_E1502, "");
      return -1;
    }

    /* NL */
    token = HECMW_ablex_next_token();
    if (token != HECMW_ABLEX_NL) {
      set_err_token(token, HECMW_IO_ABAQUS_E1500,
                    "NL required after increment");
      return -1;
    }
  } else if (token == HECMW_ABLEX_NL) {
    nod3 = 1;
  } else {
    set_err_token(token, HECMW_IO_ABAQUS_E1500,
                  "',' or NL required after nod2");
    return -1;
  }
  HECMW_assert(token == HECMW_ABLEX_NL);

  /* make node */
  if (nod1 > nod2) {
    set_err(HECMW_IO_ABAQUS_E1503,
            "Cannot generate between %d and %d with an increment of %d", nod1,
            nod2, nod3);
    return -1;
  }
  if ((nod2 - nod1) % nod3) {
    set_err(HECMW_IO_ABAQUS_E1503,
            "Cannot generate between %d and %d with an increment of %d", nod1,
            nod2, nod3);
    return -1;
  }

  n    = (nod2 - nod1) / nod3 + 1;
  node = HECMW_malloc(sizeof(int) * n);
  if (node == NULL) {
    set_err(errno, "");
    return -1;
  }

  i = 0;
  for (id = nod1; id <= nod2; id += nod3) {
    node[i++] = id;
  }
  HECMW_assert(i == n);

  *nnode      = n;
  *node_array = node;

  return 0;
}

static int read_nset(void) {
  int token, state;
  int flag_nset                 = 0; /* flag for NSET */
  int flag_generate             = 0; /* flag for GENERATE */
  int flag_unsorted             = 0; /* flag for UNSORTED */
  int isAll                     = 0;
  char nset[HECMW_NAME_LEN + 1] = "";
  enum {
    ST_FINISHED,
    ST_KEYWORD_LINE,
    ST_KEYWORD_LINE_PARAM,
    ST_DATA_LINE,
    ST_DATA_LINE_GENERATE
  };

  state = ST_KEYWORD_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_KEYWORD_LINE) {
      if (read_nset_keyword()) return -1;
      state = ST_KEYWORD_LINE_PARAM;
    } else if (state == ST_KEYWORD_LINE_PARAM) {
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_K_NSET) {
        /* must */
        if (read_nset_param_nset(nset, &isAll)) return -1;
        flag_nset = 1;
      } else if (token == HECMW_ABLEX_K_GENERATE) {
        /* oprtional */
        flag_generate = 1;
      } else if (token == HECMW_ABLEX_K_UNSORTED) {
        /* oprtional */
        log_warn(HECMW_IO_ABAQUS_W0097, "UNSORTED is not suppotred. Ignored.");
        flag_unsorted = 0; /* always ignore in this version */
      } else if (token == HECMW_ABLEX_K_INSTANCE) {
        read_nset_param_instance();
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E1500, "Unknown parameter");
        return -1;
      }

      /* check next parameter */
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_NL) {
        /* check */
        if (!flag_nset) {
          set_err(HECMW_IO_ABAQUS_E1501, "");
          return -1;
        }
        if (flag_generate) {
          state = ST_DATA_LINE_GENERATE;
        } else {
          state = ST_DATA_LINE;
        }
      } else if (token == ',') {
        ; /* continue this state */
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E1500, "Unknown parameter");
        return -1;
      }
    } else if (state == ST_DATA_LINE) {
      int n, *node;

      HECMW_assert(flag_nset);
      if (read_nset_data(&n, &node)) return -1;

      /* add node to group */
      if (HECMW_io_add_ngrp(nset, n, node) < 0) return -1;
      HECMW_free(node);

      /* check next state */
      token = HECMW_ablex_next_token();
      if (token != HECMW_ABLEX_INT) {
        state = ST_FINISHED;
      } else {
        state = ST_DATA_LINE;
      }
      HECMW_ablex_unput_token();
    } else if (state == ST_DATA_LINE_GENERATE) {
      int n, *node;

      HECMW_assert(flag_generate);
      HECMW_assert(flag_nset);

      if (read_nset_data_generate(&n, &node)) return -1;

      /* add node to group */
      if (HECMW_io_add_ngrp(nset, n, node) < 0) return -1;
      HECMW_free(node);

      /* check next state */
      token = HECMW_ablex_next_token();
      if (token != HECMW_ABLEX_INT) {
        state = ST_FINISHED;
      }
      HECMW_ablex_unput_token();
    } else {
      HECMW_assert(0);
    }
  }
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_node_keyword(int *token) {
  /* *NODE */
  *token = HECMW_ablex_next_token();
  if (*token != HECMW_ABLEX_H_NODE) {
    set_err_token(*token, HECMW_IO_ABAQUS_E1600, "*NODE required");
    return -1;
  }

  *token = HECMW_ablex_next_token();
  if (*token != HECMW_ABLEX_NL && *token != ',') {
    set_err_token(*token, HECMW_IO_ABAQUS_E1600,
                  "',' or NL required after *NODE");
    return -1;
  }
  return 0;
}

static int read_node_param_system(int *system) {
  int token;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E1600, "'=' required after SYSTEM");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != 'C' && token != 'R') {
    set_err_token(token, HECMW_IO_ABAQUS_E1600, "Invalid SYSTEM");
    return -1;
  }
  *system = token;
  return 0;
}

static int read_node_param_nset(char *nset, int *isAll) {
  char *p;
  int token;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E1600, "'=' required after NSET");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_NAME) {
    set_err_token(token, HECMW_IO_ABAQUS_E1600,
                  "NSET must begin with a letter");
    return -1;
  }
  p = HECMW_ablex_get_text();
  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }
  strcpy(nset, p);
  HECMW_toupper(nset);
  if (HECMW_io_is_reserved_name(nset)) {
    set_err(HECMW_IO_E0003, "");
    return -1;
  }
  if (strcmp(nset, "ALL") == 0) {
    HECMW_print_msg(HECMW_LOG_WARN, HECMW_IO_W1030, "");
    strcpy(nset, "ABAQUS_ESET_ALL");
    *isAll = 1;
    // return -1;
  }
  return 0;
}

static int read_node_data(int *id, double *x, double *y, double *z) {
  int token;

  /* node ID */
  *id   = 0;
  token = HECMW_ablex_next_token();
  if (token == ',') {
    HECMW_ablex_unput_token();
  } else if (token == HECMW_ABLEX_INT) {
    *id = HECMW_ablex_get_number();
  } else {
    set_err(HECMW_IO_ABAQUS_E1600, "");
    return -1;
  }
  if (*id <= 0) {
    set_err(HECMW_IO_ABAQUS_E1601, "");
    return -1;
  }

  /* ',' */
  token = HECMW_ablex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_ABAQUS_E1600, "',' required after nood ID");
    return -1;
  }

  *x = *y = *z = 0.0;
  while (1) {
    /* X */
    token = HECMW_ablex_next_token();
    if (token == HECMW_ABLEX_NL) break;
    if (token == ',') {
      HECMW_ablex_unput_token();
    } else if (token == HECMW_ABLEX_DOUBLE || token == HECMW_ABLEX_INT) {
      *x = HECMW_ablex_get_number();
    } else {
      set_err_token(token, HECMW_IO_ABAQUS_E1600, "X required");
      return -1;
    }

    /* ',' */
    token = HECMW_ablex_next_token();
    if (token == HECMW_ABLEX_NL) break;
    if (token != ',') {
      set_err_token(token, HECMW_IO_ABAQUS_E1600, "',' required after X");
      return -1;
    }

    /* Y */
    token = HECMW_ablex_next_token();
    if (token == HECMW_ABLEX_NL) break;
    if (token == ',') {
      HECMW_ablex_unput_token();
    } else if (token == HECMW_ABLEX_DOUBLE || token == HECMW_ABLEX_INT) {
      *y = HECMW_ablex_get_number();
    } else {
      set_err_token(token, HECMW_IO_ABAQUS_E1600, "Y required");
      return -1;
    }

    /* ',' */
    token = HECMW_ablex_next_token();
    if (token == HECMW_ABLEX_NL) break;
    if (token != ',') {
      set_err_token(token, HECMW_IO_ABAQUS_E1600, "',' required after Y");
      return -1;
    }

    /* Z */
    token = HECMW_ablex_next_token();
    if (token == HECMW_ABLEX_NL) break;
    if (token == HECMW_ABLEX_DOUBLE || token == HECMW_ABLEX_INT) {
      *z = HECMW_ablex_get_number();
    } else {
      set_err_token(token, HECMW_IO_ABAQUS_E1600, "Z required");
      return -1;
    }

    /* ',' or NL */
    token = HECMW_ablex_next_token();
    if (token == HECMW_ABLEX_NL) break;
    if (token == ',') {
      token = HECMW_ablex_next_token();
      if (token != HECMW_ABLEX_NL) {
        set_err_token(token, HECMW_IO_ABAQUS_E1600, "NL required after Z");
        return -1;
      }
    }

    break;
  }
  return 0;
}

static int read_node_data_system(int system, double *x, double *y, double *z) {
  struct hecmw_coord coord, result;

  /* prepare */
  coord.x = *x;
  coord.y = *y;
  coord.z = *z;

  /* reflect parameter SYSTEM */
  if (system == 'C') {
    coord.y = HECMW_degree_to_radian(coord.y);
    if (HECMW_cylindrical_to_cartesian(&coord, &result)) {
      HECMW_assert(0);
    }
    coord = result;
  }

  /* reflect *SYSTEM */
  if (HECMW_system(HECMW_io_get_system(), &coord, &result)) {
    HECMW_assert(0);
  }
  *x = result.x;
  *y = result.y;
  *z = result.z;

  return 0;
}

static int read_node(void) {
  int token, state;
  int system = 'R'; /* C:cylindrical coordinates, R:cartesian coordinates */
  int flag_system               = 0; /* flag for SYSTEM */
  int flag_nset                 = 0; /* flag for NSET */
  int flag_input                = 0; /* flag for INPUT */
  int isAll                     = 0;
  char nset[HECMW_NAME_LEN + 1] = "";
  enum {
    ST_FINISHED,
    ST_KEYWORD_LINE,
    ST_KEYWORD_LINE_PARAM,
    ST_DATA_INCLUDE,
    ST_DATA_LINE
  };

  state = ST_KEYWORD_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_KEYWORD_LINE) {
      if (read_node_keyword(&token)) return -1;
      if (token == HECMW_ABLEX_NL) {
        state = ST_DATA_LINE;
      } else if (token == ',') {
        state = ST_KEYWORD_LINE_PARAM;
      } else {
        HECMW_assert(0);
      }
    } else if (state == ST_KEYWORD_LINE_PARAM) {
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_K_SYSTEM) {
        /* optional */
        if (read_node_param_system(&system)) return -1;
        flag_system = 1;
      } else if (token == HECMW_ABLEX_K_NSET) {
        /* optional */
        if (read_node_param_nset(nset, &isAll)) return -1;
        if (isAll == 0) {
          flag_nset = 1;
        }
      } else if (token == HECMW_ABLEX_K_INPUT) {
        /* optional */
        if (read_input(HECMW_IO_ABAQUS_E1600)) return -1;
        flag_input = 1;
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E1600, "Unknown parameter");
        return -1;
      }

      /* check next parameter */
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_NL) {
        state = flag_input ? ST_DATA_INCLUDE : ST_DATA_LINE;
      } else if (token == ',') {
        ; /* continue this state */
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E1600, "Unknown parameter");
        return -1;
      }
    } else if (state == ST_DATA_INCLUDE) {
      HECMW_assert(flag_input);
      if (HECMW_ablex_switch_to_include(include_filename)) return -1;
      state = ST_DATA_LINE;
    } else if (state == ST_DATA_LINE) {
      int id;
      double x, y, z;

      if (read_node_data(&id, &x, &y, &z)) return -1;

      /* check next state */
      token = HECMW_ablex_next_token();
      if (token != HECMW_ABLEX_INT) {
        state = ST_FINISHED;
      } else {
        state = ST_DATA_LINE;
      }
      HECMW_ablex_unput_token();

      /* reflect SYSTEM, *SYSTEM */
      if (read_node_data_system(system, &x, &y, &z)) return -1;

      /* add node */
      if (HECMW_io_add_node(id, x, y, z) == NULL) return -1;

      /* add node to group */
      if (HECMW_io_add_ngrp("ALL", 1, &id) < 0) return -1;
      if (flag_nset) {
        if (HECMW_io_add_ngrp(nset, 1, &id) < 0) return -1;
      }
    } else {
      HECMW_assert(0);
    }
  }
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_shellsect_keyword(void) {
  int token;

  /* *SHELL SECTION */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_H_SHELL_SECTION) {
    set_err_token(token, HECMW_IO_ABAQUS_E1700, "*SHELL SECTION required");
    return -1;
  }

  token = HECMW_ablex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_ABAQUS_E1700,
                  "',' required after *SHELL SECTION");
    return -1;
  }
  return 0;
}

static int read_shellsect_param_elset(char *elset) {
  int token;
  char *p;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E1700, "'=' reuqired after ELSET");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_NAME) {
    set_err_token(token, HECMW_IO_ABAQUS_E1700,
                  "ELSET must begin with a letter");
    return -1;
  }
  p = HECMW_ablex_get_text();
  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }
  strcpy(elset, p);
  HECMW_toupper(elset);
  if (HECMW_io_is_reserved_name(elset)) {
    set_err(HECMW_IO_E0003, "");
    return -1;
  }
  return 0;
}

static int read_shellsect_param_material(char *material) {
  int token;
  char *p;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E1700, "'=' reuqired after MATERIAL");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_NAME) {
    set_err_token(token, HECMW_IO_ABAQUS_E1700,
                  "MATERIAL must begin with a letter");
    return -1;
  }
  p = HECMW_ablex_get_text();
  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }
  strcpy(material, p);
  HECMW_toupper(material);
  if (HECMW_io_is_reserved_name(material)) {
    set_err(HECMW_IO_E0003, "");
    return -1;
  }
  return 0;
}

static int read_shellsect_data(double *thickness, int *integpoints) {
  int token;

  /* THICKNESS */
  *thickness = 0.0;
  token      = HECMW_ablex_next_token();
  if (token == ',') {
    HECMW_ablex_unput_token();
  } else if (token == HECMW_ABLEX_DOUBLE || token == HECMW_ABLEX_INT) {
    *thickness = HECMW_ablex_get_number();
  } else {
    set_err_token(token, HECMW_IO_ABAQUS_E1700, "THICKNESS reuiqred");
    return -1;
  }
  if (*thickness <= 0.0) {
    set_err(HECMW_IO_ABAQUS_E1705, "");
    return -1;
  }

  /* ',' */
  token = HECMW_ablex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_ABAQUS_E1700, "',' required after THICKNESS");
    return -1;
  }

  /* INTEGPOINTS */
  *integpoints = 0;
  token        = HECMW_ablex_next_token();
  if (token == HECMW_ABLEX_NL) {
    HECMW_ablex_unput_token();
  } else if (token == HECMW_ABLEX_INT) {
    *integpoints = HECMW_ablex_get_number();
  } else {
    set_err_token(token, HECMW_IO_ABAQUS_E1700, "INTEGPOINTS required");
    return -1;
  }
  if (*integpoints <= 0) {
    set_err(HECMW_IO_ABAQUS_E1706, "");
    return -1;
  }

  /* NL */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_NL) {
    set_err_token(token, HECMW_IO_ABAQUS_E1700,
                  "NL required after INTEGPOINTS");
    return -1;
  }

  return 0;
}

static int read_shell_section(void) {
  int token, state;
  int composite                     = -1;
  int flag_elset                    = 0; /* flag for ELSET */
  int flag_material                 = 0; /* flag for MATERIAL */
  int flag_composite                = 0; /* flag for COMPOSITE */
  char elset[HECMW_NAME_LEN + 1]    = "";
  char material[HECMW_NAME_LEN + 1] = "ALL";
  enum { ST_FINISHED, ST_KEYWORD_LINE, ST_KEYWORD_LINE_PARAM, ST_DATA_LINE };

  state = ST_KEYWORD_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_KEYWORD_LINE) {
      if (read_shellsect_keyword()) return -1;
      state = ST_KEYWORD_LINE_PARAM;
    } else if (state == ST_KEYWORD_LINE_PARAM) {
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_K_ELSET) {
        /* must */
        if (read_shellsect_param_elset(elset)) return -1;
        flag_elset = 1;
      } else if (token == HECMW_ABLEX_K_MATERIAL) {
        /* must */
        if (flag_composite) {
          set_err(HECMW_IO_ABAQUS_E1703, "");
          return -1;
        }
        if (read_shellsect_param_material(material)) return -1;
        flag_material = 1;
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E1700, "Unknown parameter");
        return -1;
      }

      /* check next parameter */
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_NL) {
        /* check */
        if (!flag_elset) {
          set_err(HECMW_IO_ABAQUS_E1702, "");
          return -1;
        }
        if (!flag_material) {
          set_err(HECMW_IO_ABAQUS_E1707, "");
          return -1;
        }
        /* set next state */
        state = ST_DATA_LINE;
      } else if (token == ',') {
        ; /* continue this state */
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E1700, "Unknown parameter");
        return -1;
      }
    } else if (state == ST_DATA_LINE) {
      double thickness;
      int integpoints;
      struct hecmw_io_section sect;
      union hecmw_io_section_item sect_item;

      HECMW_assert(flag_elset);

      if (read_shellsect_data(&thickness, &integpoints)) return -1;

      /* set */
      sect_item.shell.thickness   = thickness;
      sect_item.shell.integpoints = integpoints;
      strcpy(sect.egrp, elset);
      strcpy(sect.material, material);
      sect.composite = composite;
      sect.secopt    = 0;
      sect.type      = HECMW_SECT_TYPE_SHELL;
      sect.sect      = sect_item;
      sect.next      = NULL;

      /* add */
      if (HECMW_io_add_sect(&sect) == NULL) return -1;

      /* set next state */
      state = ST_FINISHED;
    } else {
      HECMW_assert(0);
    }
  }
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_solidsect_keyword(void) {
  int token;

  /* *SOLID SECTION */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_H_SOLID_SECTION) {
    set_err_token(token, HECMW_IO_ABAQUS_E2100, "*SOLID SECTION required");
    return -1;
  }

  token = HECMW_ablex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_ABAQUS_E2100,
                  "',' required after *SOLID SECTION");
    return -1;
  }
  return 0;
}

static int read_solidsect_param_elset(char *elset) {
  int token;
  char *p;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E2100, "'=' reuqired after ELSET");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_NAME) {
    set_err_token(token, HECMW_IO_ABAQUS_E2100,
                  "ELSET must begin with a letter");
    return -1;
  }
  p = HECMW_ablex_get_text();
  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }
  strcpy(elset, p);
  HECMW_toupper(elset);
  if (HECMW_io_is_reserved_name(elset)) {
    set_err(HECMW_IO_E0003, "");
    return -1;
  }
  return 0;
}

static int read_solidsect_param_material(char *material) {
  int token;
  char *p;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E2100, "'=' reuqired after MATERIAL");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_NAME) {
    set_err_token(token, HECMW_IO_ABAQUS_E2100,
                  "MATERIAL must begin with a letter");
    return -1;
  }
  p = HECMW_ablex_get_text();
  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }
  strcpy(material, p);
  HECMW_toupper(material);
  if (HECMW_io_is_reserved_name(material)) {
    set_err(HECMW_IO_E0003, "");
    return -1;
  }
  return 0;
}

static int read_solidsect_param_orientation(char *orientation) {
  int token;
  char *p;

  token = HECMW_ablex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_ABAQUS_E2100,
                  "'=' reuqired after ORIENTATION");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_NAME) {
    set_err_token(token, HECMW_IO_ABAQUS_E2100,
                  "ORIENTATION must begin with a letter");
    return -1;
  }
  p = HECMW_ablex_get_text();
  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }
  strcpy(orientation, p);
  HECMW_toupper(orientation);
  if (HECMW_io_is_reserved_name(orientation)) {
    set_err(HECMW_IO_E0003, "");
    return -1;
  }
  return 0;
}

static int read_solidsect_data(double *thickness) {
  int token;

  /* THICKNESS */
  *thickness = 1.0;
  token      = HECMW_ablex_next_token();
  if (token == HECMW_ABLEX_DOUBLE || token == HECMW_ABLEX_INT) {
    *thickness = HECMW_ablex_get_number();
    /* ',' */
    token = HECMW_ablex_next_token();
    if (token != ',') {
      HECMW_ablex_unput_token();
    }
    /* NL */
    token = HECMW_ablex_next_token();
    if (token != HECMW_ABLEX_NL) {
      set_err_token(token, HECMW_IO_ABAQUS_E2100, "NL required");
      return -1;
    }
  } else {
    HECMW_ablex_unput_token();
  }

  if (*thickness <= 0.0) {
    set_err(HECMW_IO_ABAQUS_E2105, "");
    return -1;
  }

  return 0;
}

static int read_solid_section(void) {
  int token, state;
  int composite                        = -1;
  int flag_elset                       = 0; /* flag for ELSET */
  int flag_material                    = 0; /* flag for MATERIAL */
  int flag_composite                   = 0; /* flag for COMPOSITE */
  char elset[HECMW_NAME_LEN + 1]       = "";
  char material[HECMW_NAME_LEN + 1]    = "ALL";
  char orientation[HECMW_NAME_LEN + 1] = "";
  enum { ST_FINISHED, ST_KEYWORD_LINE, ST_KEYWORD_LINE_PARAM, ST_DATA_LINE };

  state = ST_KEYWORD_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_KEYWORD_LINE) {
      if (read_solidsect_keyword()) return -1;
      state = ST_KEYWORD_LINE_PARAM;
    } else if (state == ST_KEYWORD_LINE_PARAM) {
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_K_ELSET) {
        /* must */
        if (read_solidsect_param_elset(elset)) return -1;
        flag_elset = 1;
      } else if (token == HECMW_ABLEX_K_MATERIAL) {
        /* must */
        if (flag_composite) {
          set_err(HECMW_IO_ABAQUS_E2103, "");
          return -1;
        }
        if (read_solidsect_param_material(material)) return -1;
        flag_material = 1;
      } else if (token == HECMW_ABLEX_K_ORIENTATION) {
        if (read_solidsect_param_orientation(orientation)) return -1;
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E2100, "Unknown parameter");
        return -1;
      }

      /* check next parameter */
      token = HECMW_ablex_next_token();
      if (token == HECMW_ABLEX_NL) {
        /* check */
        if (!flag_elset) {
          set_err(HECMW_IO_ABAQUS_E2102, "");
          return -1;
        }
        if (!flag_material) {
          set_err(HECMW_IO_ABAQUS_E2108, "");
          return -1;
        }
        /* set next state */
        state = ST_DATA_LINE;
      } else if (token == ',') {
        ; /* continue this state */
      } else {
        set_err_token(token, HECMW_IO_ABAQUS_E2100, "Unknown parameter");
        return -1;
      }
    } else if (state == ST_DATA_LINE) {
      int secopt;
      double thickness;
      struct hecmw_io_section sect;
      union hecmw_io_section_item sect_item;

      HECMW_assert(flag_elset);

      if (read_solidsect_data(&thickness)) return -1;

      secopt = get_secopt(elset, HECMW_IO_ABAQUS_E2107);
      if (secopt < 0) {
        set_err_noloc(HECMW_IO_E1026, "Two or more secopt found in %s", elset);
        return -1;
      }

      /* set */
      sect_item.solid.thickness = thickness;
      strcpy(sect.egrp, elset);
      strcpy(sect.material, material);
      sect.composite = composite;
      sect.secopt    = secopt;
      sect.type      = HECMW_SECT_TYPE_SOLID;
      sect.sect      = sect_item;
      sect.next      = NULL;

      /* add */
      if (HECMW_io_add_sect(&sect) == NULL) return -1;

      /* set next state */
      state = ST_FINISHED;
    } else {
      HECMW_assert(0);
    }
  }
  return 0;
}

/*----------------------------------------------------------------------------*/
#if 0
static int
read_system_keyword(void)
{
	int token;

	/* *SYSTEM */
	token = HECMW_ablex_next_token();
	if(token != HECMW_ABLEX_H_SYSTEM) {
		set_err_token(token, HECMW_IO_ABAQUS_E1900, "*SYSTEM required");
		return -1;
	}

	/* NL */
	token = HECMW_ablex_next_token();
	if(token != HECMW_ABLEX_NL) {
		set_err_token(token, HECMW_IO_ABAQUS_E1900, "NL required after *SYSTEM");
		return -1;
	}
	return 0;
}


static int
read_system_data_line1a(struct hecmw_system_param *system, int *last_token)
{
	int token;

	/* Xa */
	token = HECMW_ablex_next_token();
	if(token == HECMW_ABLEX_DOUBLE || token == HECMW_ABLEX_INT) {
		system->xa = HECMW_ablex_get_number();
	} else if(token == ',') {
		system->xa = 0.0;
		HECMW_ablex_unput_token();
	} else {
		set_err_token(token, HECMW_IO_ABAQUS_E1900, "Xa required");
		return -1;
	}

	/* ',' */
	token = HECMW_ablex_next_token();
	if(token != ',') {
		set_err_token(token, HECMW_IO_ABAQUS_E1900, "',' required after Xa");
		return -1;
	}

	/* Ya */
	token = HECMW_ablex_next_token();
	if(token == HECMW_ABLEX_DOUBLE || token == HECMW_ABLEX_INT) {
		system->ya = HECMW_ablex_get_number();
	} else if(token == ',') {
		system->ya = 0.0;
		HECMW_ablex_unput_token();
	} else {
		set_err_token(token, HECMW_IO_ABAQUS_E1900, "Ya required");
		return -1;
	}

	/* ',' */
	token = HECMW_ablex_next_token();
	if(token != ',') {
		set_err_token(token, HECMW_IO_ABAQUS_E1900, "',' required after Ya");
		return -1;
	}

	/* Za */
	token = HECMW_ablex_next_token();
	if(token == HECMW_ABLEX_DOUBLE || token == HECMW_ABLEX_INT) {
		system->za = HECMW_ablex_get_number();
	} else if(token == ',' || token == HECMW_ABLEX_NL) {
		system->za = 0.0;
		HECMW_ablex_unput_token();
	} else {
		set_err_token(token, HECMW_IO_ABAQUS_E1900, "Za required");
		return -1;
	}

	token = HECMW_ablex_next_token();
	if(token != ',' && token != HECMW_ABLEX_NL) {
		set_err_token(token, HECMW_IO_ABAQUS_E1900, "',' or NL required after Za");
		return -1;
	}

	*last_token = token;

	return 0;
}


static int
read_system_data_line1b(struct hecmw_system_param *system)
{
	int token;
	/* Xb */
	token = HECMW_ablex_next_token();
	if(token == HECMW_ABLEX_DOUBLE || token == HECMW_ABLEX_INT) {
		system->xb = HECMW_ablex_get_number();
	} else if(token == ',') {
		system->xb = 0.0;
		HECMW_ablex_unput_token();
	} else {
		set_err_token(token, HECMW_IO_ABAQUS_E1900, "Xb required");
		return -1;
	}

	/* ',' */
	token = HECMW_ablex_next_token();
	if(token != ',') {
		set_err_token(token, HECMW_IO_ABAQUS_E1900, "',' required after Xb");
		return -1;
	}

	/* Yb */
	token = HECMW_ablex_next_token();
	if(token == HECMW_ABLEX_DOUBLE || token == HECMW_ABLEX_INT) {
		system->yb = HECMW_ablex_get_number();
	} else if(token == ',') {
		system->yb = 0.0;
		HECMW_ablex_unput_token();
	} else {
		set_err_token(token, HECMW_IO_ABAQUS_E1900, "Yb required");
		return -1;
	}

	/* ',' */
	token = HECMW_ablex_next_token();
	if(token != ',') {
		set_err_token(token, HECMW_IO_ABAQUS_E1900, "',' required after Yb");
		return -1;
	}

	/* Zb */
	token = HECMW_ablex_next_token();
	if(token == HECMW_ABLEX_DOUBLE || token == HECMW_ABLEX_INT) {
		system->zb = HECMW_ablex_get_number();
	} else if(token == HECMW_ABLEX_NL) {
		system->zb = 0.0;
		HECMW_ablex_unput_token();
	} else {
		set_err_token(token, HECMW_IO_ABAQUS_E1900, "Zb required");
		return -1;
	}

	/*NL */
	token = HECMW_ablex_next_token();
	if(token != HECMW_ABLEX_NL) {
		set_err_token(token, HECMW_IO_ABAQUS_E1900, "NL required after Zb");
		return -1;
	}

	return 0;
}


static int
read_system_data_line2(struct hecmw_system_param *system)
{
	int token;

	/* Xc */
	token = HECMW_ablex_next_token();
	if(token == HECMW_ABLEX_DOUBLE || token == HECMW_ABLEX_INT) {
		system->xc = HECMW_ablex_get_number();
	} else if(token == ',') {
		system->xc = 0.0;
		HECMW_ablex_unput_token();
	} else {
		set_err_token(token, HECMW_IO_ABAQUS_E1900, "Xc required");
		return -1;
	}

	/* ',' */
	token = HECMW_ablex_next_token();
	if(token != ',') {
		set_err_token(token, HECMW_IO_ABAQUS_E1900, "',' required after Xc");
		return -1;
	}

	/* Yc */
	token = HECMW_ablex_next_token();
	if(token == HECMW_ABLEX_DOUBLE || token == HECMW_ABLEX_INT) {
		system->yc = HECMW_ablex_get_number();
	} else if(token == ',') {
		system->yc = 0.0;
		HECMW_ablex_unput_token();
	} else {
		set_err_token(token, HECMW_IO_ABAQUS_E1900, "Yc required");
		return -1;
	}

	/* ',' */
	token = HECMW_ablex_next_token();
	if(token != ',') {
		set_err_token(token, HECMW_IO_ABAQUS_E1900, "',' required after Yc");
		return -1;
	}

	/* Zc */
	token = HECMW_ablex_next_token();
	if(token == HECMW_ABLEX_DOUBLE || token == HECMW_ABLEX_INT) {
		system->zc = HECMW_ablex_get_number();
	} else if(token == HECMW_ABLEX_NL) {
		system->zc = 0.0;
		HECMW_ablex_unput_token();
	} else {
		set_err_token(token, HECMW_IO_ABAQUS_E1900, "Zc required");
		return -1;
	}

	/* NL */
	token = HECMW_ablex_next_token();
	if(token != HECMW_ABLEX_NL) {
		set_err_token(token, HECMW_IO_ABAQUS_E1900, "NL required after Zc");
		return -1;
	}

	return 0;
}


static int
read_system(void)
{
	int token,state;
	struct hecmw_system_param *system = NULL;
	enum {
		ST_FINISHED,
		ST_KEYWORD_LINE,
		ST_DATA_LINE1,
		ST_DATA_LINE2
	};

	fprintf(stderr, "*SYSTEM has not implemented yet\n");
	HECMW_abort(HECMW_comm_get_comm());

	system = HECMW_malloc(sizeof(*system));
	if(system == NULL) {
		set_err(errno, "");
		return -1;
	}

	/* default values */
	system->xa = 0.0;
	system->ya = 0.0;
	system->za = 0.0;
	system->xb = 0.0;
	system->yb = 0.0;
	system->zb = 0.0;
	system->xc = 0.0;
	system->yc = 0.0;
	system->zc = 0.0;

	state = ST_KEYWORD_LINE;
	while(state != ST_FINISHED) {
		if(state == ST_KEYWORD_LINE) {
			if(read_system_keyword()) return -1;
			/* check next state */
			token = HECMW_ablex_next_token();
			if(token != HECMW_ABLEX_DOUBLE && token != HECMW_ABLEX_INT && token != ',') {
				/* clear *SYSTEM */
				HECMW_free(system);
				system = NULL;
				state = ST_FINISHED;
			} else {
				state = ST_DATA_LINE1;
			}
			HECMW_ablex_unput_token();
		} else if(state == ST_DATA_LINE1) {
			if(read_system_data_line1a(system, &token)) return -1;
			if(token == HECMW_ABLEX_NL) {
				state = ST_FINISHED;
				continue;
			}
			HECMW_assert(token == ',');

			if(read_system_data_line1b(system)) return -1;
			token = HECMW_ablex_next_token();
			if(token != HECMW_ABLEX_DOUBLE && token != HECMW_ABLEX_INT && token != ',') {
				state = ST_FINISHED;
			} else {
				state = ST_DATA_LINE2;
			}
			HECMW_ablex_unput_token();
		} else if(state == ST_DATA_LINE2) {
			if(read_system_data_line2(system)) return -1;
			state = ST_FINISHED;
		} else {
			HECMW_assert(0);
		}
	}

	/* set */
	HECMW_io_set_system(system);

	return 0;
}
#endif

/*----------------------------------------------------------------------------*/
static int read_boundary_keyword(void) {
  int token;
  /* static int isFirst = 0; */

  /* *BOUNDARY */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_H_BOUNDARY) {
    set_err_token(token, HECMW_IO_ABAQUS_E1500, "*BOUNDARY required");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_NL) {
    set_err_token(token, HECMW_IO_ABAQUS_E1700,
                  "',' is not required after *BOUNDARY SECTION");
    return -1;
  }
  /* if(isFirst == 0) { */
  fprintf(stderr,
          "Auto-generated cards should be added in !BOUNDARY section of *.cnt "
          "file \n");
  /*	isFirst =1; */
  /*} */
  return 0;
}

static int read_boundary_data(int *nnode, int **node_array) {
  int i, n, *node, token;
  int isFirst, isSuggest, isNode;
  struct hecmw_io_id *head, *prev, *p, *q;

  n         = 0;
  isNode    = 0;
  isSuggest = 0;
  head      = NULL;
  node      = NULL;
  while (1) {
    struct hecmw_io_id *id;

    token = HECMW_ablex_next_token();
    if (token == HECMW_ABLEX_NL) break;

    id = HECMW_malloc(sizeof(*id));
    if (id == NULL) {
      set_err(errno, "");
      goto error;
    }

    /* node or NGRP */
    if (token != HECMW_ABLEX_INT && token != HECMW_ABLEX_DOUBLE &&
        token != ',') {
      isNode = 1;
    }

    if (token == HECMW_ABLEX_NAME) {
      isSuggest = 1;
      char *c   = HECMW_ablex_get_text();
      fprintf(stderr, "%s", c);
    }
    if (token == HECMW_ABLEX_INT && isSuggest == 1) {
      int i = HECMW_ablex_get_number();
      fprintf(stderr, ", %d", i);
    }
    if (token == HECMW_ABLEX_DOUBLE && isSuggest == 1) {
      double a = HECMW_ablex_get_number();
      fprintf(stderr, ", %f", a);
    }

    if (token == ',') {
      id->id = 0;
      HECMW_ablex_unput_token();
    } else if (token == HECMW_ABLEX_INT) {
      id->id = HECMW_ablex_get_number();
    }

    if (isNode == 0) {
      id->next = NULL;
      if (head == NULL) {
        head = id;
      } else {
        prev->next = id;
      }
      prev = id;
      n++;
      isNode = 1;

      HECMW_assert(head);
      HECMW_assert(n > 0);
    }

    /* ',' or NL */
    token = HECMW_ablex_next_token();
    if (token != ',' && token != HECMW_ABLEX_NL) {
      set_err_token(token, HECMW_IO_ABAQUS_E1500,
                    "',' or NL required after node ID");
      goto error;
    }
    if (token == HECMW_ABLEX_NL) break;
  }

  node = HECMW_malloc(sizeof(*node) * n);
  if (node == NULL) {
    set_err(errno, "");
    goto error;
  }
  i = 0;
  for (p = head; p; p = q) {
    q         = p->next;
    node[i++] = p->id;
    HECMW_free(p);
  }
  head   = NULL;
  isNode = 0;

  *nnode      = n;
  *node_array = node;

  if (isSuggest == 1) {
    fprintf(stderr, "\n");
    isSuggest = 0;
  }

  return 0;
error:
  for (p = head; p; p = q) {
    q = p->next;
    HECMW_free(p);
  }
  HECMW_free(node);
  return -1;
}

static int read_boundary(void) {
  int token, state;
  int flag_boundary             = 0; /* flag for BOUNDARY */
  int isNodeInput               = 0;
  static int nbound             = 0;
  char nset[HECMW_NAME_LEN + 1] = "";
  enum { ST_FINISHED, ST_KEYWORD_LINE, ST_KEYWORD_LINE_PARAM, ST_DATA_LINE };

  state = ST_KEYWORD_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_KEYWORD_LINE) {
      if (read_boundary_keyword()) return -1;
      state         = ST_DATA_LINE;
      flag_boundary = 1;
    } else if (state == ST_KEYWORD_LINE_PARAM) {
      ;
    } else if (state == ST_DATA_LINE) {
      int n = 0;
      int *node;

      HECMW_assert(flag_boundary);
      if (read_boundary_data(&n, &node)) return -1;

      if (n > 0) {
        isNodeInput = 1;

        /* add node to group */
        sprintf(nset, "BND%d", nbound);
        if (HECMW_io_add_ngrp(nset, n, node) < 0) return -1;
        HECMW_free(node);
      }

      /* check next state */
      token = HECMW_ablex_next_token();
      if (token != HECMW_ABLEX_INT && token != HECMW_ABLEX_NAME) {
        state = ST_FINISHED;
      } else {
        state = ST_DATA_LINE;
      }
      HECMW_ablex_unput_token();
    } else {
      HECMW_assert(0);
    }
  }

  if (isNodeInput == 1) {
    fprintf(stderr, "NGRP=BND%d\n", nbound);
    nbound++;
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
static int read_cload_keyword(void) {
  int token;
  static int isFirst = 0;

  /* *BOUNDARY */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_H_CLOAD) {
    set_err_token(token, HECMW_IO_ABAQUS_E1500, "*CLOAD required");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_NL) {
    set_err_token(token, HECMW_IO_ABAQUS_E1700,
                  "',' is not required after *CLOAD SECTION");
    return -1;
  }
  /*if(isFirst == 0) {*/
  fprintf(stderr,
          "Auto-generated cards should be added in !CLOAD section of *.cnt "
          "file \n");
  isFirst = 1;
  /*}*/
  return 0;
}

static int read_cload(void) {
  int token, state;
  int flag_cload                = 0; /* flag for CLOAD */
  int isCload                   = 0;
  static int ncload             = 0;
  char nset[HECMW_NAME_LEN + 1] = "";
  enum { ST_FINISHED, ST_KEYWORD_LINE, ST_KEYWORD_LINE_PARAM, ST_DATA_LINE };

  state = ST_KEYWORD_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_KEYWORD_LINE) {
      if (read_cload_keyword()) return -1;
      state      = ST_DATA_LINE;
      flag_cload = 1;
    } else if (state == ST_KEYWORD_LINE_PARAM) {
      ;
    } else if (state == ST_DATA_LINE) {
      int n, *node;

      HECMW_assert(flag_cload);
      if (read_boundary_data(&n, &node)) return -1;

      /* add node to group */
      if (n != 0) {
        isCload = 1;

        sprintf(nset, "CLOAD%d", ncload);
        if (HECMW_io_add_ngrp(nset, n, node) < 0) return -1;
        HECMW_free(node);
      }

      /* check next state */
      token = HECMW_ablex_next_token();
      if (token != HECMW_ABLEX_INT) {
        state = ST_FINISHED;
      } else {
        state = ST_DATA_LINE;
      }
      HECMW_ablex_unput_token();
    } else {
      HECMW_assert(0);
    }
  }

  if (isCload != 0) {
    fprintf(stderr, "NGRP=CLOAD%d\n", ncload);
    ncload++;
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
static int read_dload_keyword(void) {
  int token;
  static int isFirst = 0;

  /* *DLOAD */
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_H_DLOAD) {
    set_err_token(token, HECMW_IO_ABAQUS_E1500, "*DLOAD required");
    return -1;
  }
  token = HECMW_ablex_next_token();
  if (token != HECMW_ABLEX_NL) {
    set_err_token(token, HECMW_IO_ABAQUS_E1700,
                  "',' is not required after *DLOAD SECTION");
    return -1;
  }
  /*if(isFirst == 0) {*/
  fprintf(stderr,
          "Auto-generated cards should be added in !DLOAD section of *.cnt "
          "file \n");
  isFirst = 1;
  /*}*/
  return 0;
}

static int read_dload(void) {
  int token, state;
  int flag_dload                 = 0; /* flag for DLOAD */
  int isDload                    = 0;
  static int ndload              = 0;
  char elset[HECMW_NAME_LEN + 1] = "";
  enum { ST_FINISHED, ST_KEYWORD_LINE, ST_KEYWORD_LINE_PARAM, ST_DATA_LINE };

  state = ST_KEYWORD_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_KEYWORD_LINE) {
      if (read_dload_keyword()) return -1;
      state      = ST_DATA_LINE;
      flag_dload = 1;
    } else if (state == ST_KEYWORD_LINE_PARAM) {
      ;
    } else if (state == ST_DATA_LINE) {
      int n, *elem;

      HECMW_assert(flag_dload);
      if (read_boundary_data(&n, &elem)) return -1;

      if (n != 0) {
        isDload = 1;
        /* add node to group */
        sprintf(elset, "DLOAD%d", ndload);
        if (HECMW_io_add_egrp(elset, n, elem) < 0) return -1;
        HECMW_free(elem);
      }

      /* check next state */
      token = HECMW_ablex_next_token();
      if (token != HECMW_ABLEX_INT) {
        state = ST_FINISHED;
      } else {
        state = ST_DATA_LINE;
      }
      HECMW_ablex_unput_token();
    } else {
      HECMW_assert(0);
    }
  }

  if (isDload != 0) {
    fprintf(stderr, "NGRP=DLOAD%d is automatically generated\n", ndload);
    ndload++;
  }

  return 0;
}

/*------------------------------------------------------------------------------
  ReadFunc table
*/

typedef int (*ReadFunc)(void);

static struct read_func_table {
  int token;
  ReadFunc func;
} read_func_table[] = {
    {HECMW_ABLEX_H_AMPLITUDE, read_amplitude},
    {HECMW_ABLEX_H_BOUNDARY, read_boundary},
    {HECMW_ABLEX_H_CONDUCTIVITY, read_conductivity},
    {HECMW_ABLEX_H_CLOAD, read_cload},
    {HECMW_ABLEX_H_DLOAD, read_dload},
    {HECMW_ABLEX_H_DENSITY, read_density},
    /*	{ HECMW_ABLEX_H_ECOPY,     read_ecopy     }, */
    /*	{ HECMW_ABLEX_H_EGEN,      read_egen      }, */
    {HECMW_ABLEX_H_ELASTIC, read_elastic},
    {HECMW_ABLEX_H_ELEMENT, read_element},
    {HECMW_ABLEX_H_ELSET, read_elset},
    {HECMW_ABLEX_H_EQUATION, read_equation},
    {HECMW_ABLEX_H_HEADING, read_heading},
    {HECMW_ABLEX_H_INCLUDE, read_include},
    {HECMW_ABLEX_H_INITIAL, read_initial},
    {HECMW_ABLEX_H_MATERIAL, read_material},
    /*	{ HECMW_ABLEX_H_NCOPY,     read_ncopy     }, */
    /*	{ HECMW_ABLEX_H_NFILL,     read_nfill     }, */
    /*	{ HECMW_ABLEX_H_NGEN,      read_ngen      }, */
    {HECMW_ABLEX_H_NODE, read_node},
    {HECMW_ABLEX_H_NSET, read_nset},
    {HECMW_ABLEX_H_SHELL_SECTION, read_shell_section},
    {HECMW_ABLEX_H_SOLID_SECTION, read_solid_section},
    {HECMW_ABLEX_H_SPECIFIC_HEAT, read_specific_heat},
    /*	{ HECMW_ABLEX_H_SYSTEM,    read_system    }, */
};

#define N_READ_FUNC (sizeof(read_func_table) / sizeof(read_func_table[0]))

/* static int (* get_read_func(int token))(void) */
static ReadFunc get_read_func(int token) {
  int i;

  for (i = 0; i < N_READ_FUNC; i++) {
    if (token == read_func_table[i].token) {
      return read_func_table[i].func;
    }
  }
  return NULL;
}

static int parse(void) {
  int token, head;
  ReadFunc func;

  while ((token = HECMW_ablex_next_token())) {
    if (token == HECMW_ABLEX_NL) continue;
    head = token;
    func = get_read_func(token);
    if (func == NULL) {
      char *p = HECMW_ablex_get_text();
      if (p[0] != '*') {
        set_err(HECMW_IO_ABAQUS_E0098, "");
        return -1;
      }
      /* skip unsupported keyword */
      token = HECMW_ablex_next_token();
      p     = token ? HECMW_ablex_get_text() : "";
      HECMW_print_msg(HECMW_LOG_WARN, HECMW_IO_ABAQUS_W0099, "*%s", p);
      if (!token) break;
      while ((token = HECMW_ablex_next_token())) {
        p = HECMW_ablex_get_text();
        if (p[0] == '*') break;
      }
      if (!token) break;
      HECMW_ablex_unput_token(); /* unput *XXXX */
      continue;
    }
    if (is_material_keyword(token) && !is_material_zone()) {
      set_err(HECMW_IO_ABAQUS_E0096, "keyword: %s", HECMW_ablex_get_text());
      return -1;
    }
    if (!is_material_keyword(token) && is_material_zone()) {
      if (regist_material()) return -1;
      set_material_zone(0);
    }
    HECMW_ablex_unput_token(); /* unput *XXXX */
    if ((*func)()) return -1;
    if (head == HECMW_ABLEX_H_MATERIAL) {
      set_material_zone(1);
    }
  }
  if (is_material_zone()) {
    if (regist_material()) return -1;
    set_material_zone(0);
  }
  return 0;
}

/*----------------------------------------------------------------------------*/

static int post_abaqus(void) {
  HECMW_map_int_finalize(elem_secopt);
  HECMW_free(elem_secopt);
  return 0;
}

/* read only. Not make hecmwST_local_mesh */
int HECMW_read_abaqus_mesh(const char *filename) {
  FILE *fp;

  HECMW_log(HECMW_LOG_DEBUG, "Start to read ABAQUS mesh");

  if (filename == NULL) {
    set_err_noloc(HECMW_IO_E0001,
                  "Not specified filename for ABAQUS mesh input routine");
    return -1;
  }
  HECMW_log(HECMW_LOG_DEBUG, "ABAQUS mesh file is '%s'", filename);

  if (strlen(filename) > HECMW_FILENAME_LEN) {
    set_err_noloc(HECMW_IO_E0002, "");
    return -1;
  }

  strcpy(grid_filename, filename);
  HECMW_io_set_gridfile(grid_filename);

  if ((fp = fopen(filename, "r")) == NULL) {
    set_err_noloc(HECMW_IO_ABAQUS_E0001, "File: %s, %s", filename,
                  strerror(errno));
    return -1;
  }

  if (HECMW_ablex_set_input(fp)) return -1;

  HECMW_log(HECMW_LOG_DEBUG, "Parsing...");
  if (parse()) {
    return -1;
  }

  if (fclose(fp)) {
    set_err_noloc(HECMW_IO_ABAQUS_E0002, "File: %s, %s", filename,
                  strerror(errno));
    return -1;
  }

  if (post_abaqus()) {
    return -1;
  }

  strcpy(grid_filename, "Unknown");

  return 0;
}

struct hecmwST_local_mesh *HECMW_get_abaqus_mesh(const char *filename) {
  struct hecmwST_local_mesh *local_mesh;

  if (HECMW_io_init()) return NULL;
  if (HECMW_io_pre_process()) return NULL;
  if (HECMW_read_abaqus_mesh(filename)) return NULL;
  if (HECMW_io_post_process()) return NULL;
  local_mesh = HECMW_io_make_local_mesh();
  if (local_mesh == NULL) return NULL;
  if (HECMW_io_finalize()) return NULL;

  strcpy(grid_filename, "Unknown");

  return local_mesh;
}
