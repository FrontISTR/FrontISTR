/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "hecmw_util.h"
#include "hecmw_control.h"
#include "hecmw_ctrllex.h"
#include "hecmw_path.h"

static char ctrl_filename[HECMW_FILENAME_LEN + 1];

static int table_entire_mesh[] = {
    HECMW_CTRL_FTYPE_HECMW_ENTIRE, HECMW_CTRL_FTYPE_GEOFEM,
    HECMW_CTRL_FTYPE_ABAQUS,       HECMW_CTRL_FTYPE_NASTRAN,
    HECMW_CTRL_FTYPE_FEMAP,
};

struct mesh_entry {
  char *name_ID;
  int type;
  int io;
  int refine;
  char *filename;
  struct mesh_entry *next;
};

struct mesh_grp_entry {
  char *name_ID;
  int n_mesh;
  struct mesh_entry **mesh;
  struct mesh_grp_entry *next;
};

struct restart_entry {
  char *name_ID;
  int io;
  /*#define HECMW_CTRL_FILE_IO_IN 1*/
  /*#define HECMW_CTRL_FILE_IO_OUT 2*/
  /* #define HECMW_CTRL_FILE_IO_INOUT 4 */
  char *filename;
  struct restart_entry *next;
};

struct result_entry {
  char *name_ID;
  int io;
  /*#define HECMW_CTRL_FILE_IO_IN 1*/
  /*#define HECMW_CTRL_FILE_IO_OUT 2*/
  int fg_text; /* 1:text(default), 0:binary */
  char *filename;
  struct result_entry *next;
};

struct ctrl_entry {
  char *name_ID;
  char *filename;
  struct ctrl_entry *next;
};

static int subdir_on = 0;

static int nlimit;

static struct ctrl_entry *ctrl_ent;

static struct mesh_entry *mesh_ent;

static struct mesh_grp_entry *mesh_grp_ent;

static struct restart_entry *restart_ent;

static struct result_entry *result_ent;

/*----------------------------------------------------------------------------*/

static int is_entire_mesh(int type) {
  int i, n;
  n = sizeof(table_entire_mesh) / sizeof(table_entire_mesh[0]);

  for (i = 0; i < n; i++) {
    if (table_entire_mesh[i] == type) return 1;
  }

  return 0;
}

/* return dir + subdir + prefix + file + suffix + .rank */
static char *make_filename(char *dir, char *subdir, char *prefix, char *file,
                           char *suffix, int myrank, int flag_rank) {
  static char filename[HECMW_FILENAME_LEN + 1];
  char rank[10];
  char separator[10];
  strcpy(filename, "");

  if (dir && strlen(dir) > 0) {
    sprintf(separator, "%c", HECMW_get_path_separator());

    if ((strlen(dir) + strlen(separator)) > HECMW_FILENAME_LEN) return NULL;

    sprintf(filename, "%s%s", dir, separator);
  }

  if (subdir && strlen(subdir) > 0) {
    sprintf(separator, "%c", HECMW_get_path_separator());

    if ((strlen(filename) + strlen(subdir) + strlen(separator)) >
        HECMW_FILENAME_LEN)
      return NULL;

    strcat(filename, subdir);
    strcat(filename, separator);
  }

  if (prefix && strlen(prefix) > 0) {
    sprintf(separator, "%c", HECMW_get_path_separator());

    if ((strlen(filename) + strlen(prefix) + strlen(separator)) >
        HECMW_FILENAME_LEN)
      return NULL;

    strcat(filename, prefix);
    strcat(filename, separator);
  }

  if ((strlen(filename) + strlen(file)) > HECMW_FILENAME_LEN) return NULL;

  strcat(filename, file);

  if (suffix) {
    if ((strlen(filename) + strlen(suffix)) > HECMW_FILENAME_LEN) return NULL;

    strcat(filename, suffix);
  }

  if (flag_rank) {
    sprintf(rank, ".%d", myrank);

    if ((strlen(filename) + strlen(rank)) > HECMW_FILENAME_LEN) return NULL;

    strcat(filename, rank);
  }

  return filename;
}

/* return dir + subdir + prefix + file + suffix + .rank (thread-safe version) */
static char *make_filename_r(char *dir, char *subdir, char *prefix, char *file,
                             char *suffix, int myrank, int flag_rank,
                             char *filename, int filename_len) {
  char rank[10];
  char separator[10];
  strcpy(filename, "");

  HECMW_assert( filename );

  if (dir && strlen(dir) > 0) {
    sprintf(separator, "%c", HECMW_get_path_separator());

    if ((strlen(dir) + strlen(separator)) > filename_len) return NULL;

    sprintf(filename, "%s%s", dir, separator);
  }

  if (subdir && strlen(subdir) > 0) {
    sprintf(separator, "%c", HECMW_get_path_separator());

    if ((strlen(filename) + strlen(subdir) + strlen(separator)) >
        filename_len)
      return NULL;

    strcat(filename, subdir);
    strcat(filename, separator);
  }

  if (prefix && strlen(prefix) > 0) {
    sprintf(separator, "%c", HECMW_get_path_separator());

    if ((strlen(filename) + strlen(prefix) + strlen(separator)) >
        filename_len)
      return NULL;

    strcat(filename, prefix);
    strcat(filename, separator);
  }

  if ((strlen(filename) + strlen(file)) > filename_len) return NULL;

  strcat(filename, file);

  if (suffix) {
    if ((strlen(filename) + strlen(suffix)) > filename_len) return NULL;

    strcat(filename, suffix);
  }

  if (flag_rank) {
    sprintf(rank, ".%d", myrank);

    if ((strlen(filename) + strlen(rank)) > filename_len) return NULL;

    strcat(filename, rank);
  }

  return filename;
}

/*----------------------------------------------------------------------------*/

static void free_mesh_entry(void) {
  struct mesh_entry *p, *q;

  for (p = mesh_ent; p; p = q) {
    q = p->next;
    HECMW_free(p->name_ID);
    HECMW_free(p->filename);
    HECMW_free(p);
  }

  mesh_ent = NULL;
}

static void free_mesh_grp_entry(void) {
  struct mesh_grp_entry *p, *q;

  for (p = mesh_grp_ent; p; p = q) {
    q = p->next;
    HECMW_free(p->name_ID);
    HECMW_free(p->mesh); /* free only mesh array */
    HECMW_free(p);
  }

  mesh_grp_ent = NULL;
}

static void free_restart_entry(void) {
  struct restart_entry *p, *q;

  for (p = restart_ent; p; p = q) {
    q = p->next;
    HECMW_free(p->name_ID);
    HECMW_free(p->filename);
    HECMW_free(p);
  }

  restart_ent = NULL;
}

static void free_result_entry(void) {
  struct result_entry *p, *q;

  for (p = result_ent; p; p = q) {
    q = p->next;
    HECMW_free(p->name_ID);
    HECMW_free(p->filename);
    HECMW_free(p);
  }

  result_ent = NULL;
}

static void free_ctrl_entry(void) {
  struct ctrl_entry *p, *q;

  for (p = ctrl_ent; p; p = q) {
    q = p->next;
    HECMW_free(p->name_ID);
    HECMW_free(p->filename);
    HECMW_free(p);
  }

  ctrl_ent = NULL;
}

/*----------------------------------------------------------------------------*/

static struct mesh_entry *get_mesh_entry(char *name_ID) {
  struct mesh_entry *p;

  if (name_ID == NULL) return NULL;

  for (p = mesh_ent; p; p = p->next) {
    if (strcmp(p->name_ID, name_ID) == 0) return p;
  }

  return NULL;
}

static struct mesh_grp_entry *get_mesh_grp_entry(char *name_ID) {
  struct mesh_grp_entry *p;

  if (name_ID == NULL) return NULL;

  for (p = mesh_grp_ent; p; p = p->next) {
    if (strcmp(p->name_ID, name_ID) == 0) return p;
  }

  return NULL;
}

static struct restart_entry *get_restart_entry(char *name_ID) {
  struct restart_entry *p;

  if (name_ID == NULL) return NULL;

  for (p = restart_ent; p; p = p->next) {
    if (strcmp(p->name_ID, name_ID) == 0) return p;
  }

  return NULL;
}

static struct restart_entry *get_restart_entry_by_io(int io) {
  struct restart_entry *p;

  for (p = restart_ent; p; p = p->next) {
    if (p->io & io) return p; /* attention!! arg io is bitmap */
  }

  return NULL;
}

static struct result_entry *get_result_entry(char *name_ID) {
  struct result_entry *p;

  if (name_ID == NULL) return NULL;

  for (p = result_ent; p; p = p->next) {
    if (strcmp(p->name_ID, name_ID) == 0) return p;
  }

  return NULL;
}

static struct ctrl_entry *get_ctrl_entry(char *name_ID) {
  struct ctrl_entry *p;

  if (name_ID == NULL) return NULL;

  for (p = ctrl_ent; p; p = p->next) {
    if (strcmp(p->name_ID, name_ID) == 0) return p;
  }

  return NULL;
}

/*----------------------------------------------------------------------------*/

static struct mesh_entry *make_mesh_entry(char *name_ID, int type, int io,
                                          int refine, char *filename) {
  char *p;
  struct mesh_entry *mesh = NULL;
  mesh                    = HECMW_calloc(1, sizeof(*mesh));

  if (mesh == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  mesh->type   = type;
  mesh->io     = io;
  mesh->refine = refine;
  mesh->next   = NULL;
  p            = HECMW_strdup(name_ID);

  if (p == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  mesh->name_ID = p;
  p             = HECMW_strdup(filename);

  if (p == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  mesh->filename = p;
  return mesh;
error:

  if (mesh) {
    HECMW_free(mesh->name_ID);
    HECMW_free(mesh->filename);
    HECMW_free(mesh);
  }

  return NULL;
}

static int add_mesh_entry(struct mesh_entry *mesh) {
  struct mesh_entry *p, *q;
  q = NULL;

  for (p = mesh_ent; p; p = (q = p)->next)
    ;

  if (q == NULL) {
    mesh_ent = mesh;

  } else {
    q->next = mesh;
  }

  return 0;
}

static struct mesh_grp_entry *make_mesh_group_entry(char *name_ID, int n_mesh,
                                                    char **mesh) {
  int i;
  struct mesh_grp_entry *meshgrp = NULL;
  meshgrp                        = HECMW_calloc(1, sizeof(*meshgrp));

  if (meshgrp == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  meshgrp->name_ID = HECMW_strdup(name_ID);

  if (meshgrp->name_ID == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  meshgrp->n_mesh = n_mesh;
  meshgrp->mesh   = HECMW_calloc(n_mesh, sizeof(*meshgrp->mesh));

  if (meshgrp->mesh == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  for (i = 0; i < n_mesh; i++) {
    meshgrp->mesh[i] = get_mesh_entry(mesh[i]);

    if (meshgrp->mesh[i] == NULL) goto error;
  }

  meshgrp->next = NULL;
  return meshgrp;
error:

  if (meshgrp) {
    HECMW_free(meshgrp->name_ID);
    HECMW_free(meshgrp->mesh);
    HECMW_free(meshgrp);
  }

  return NULL;
}

static int add_mesh_group_entry(struct mesh_grp_entry *mesh) {
  struct mesh_grp_entry *p, *q;
  q = NULL;

  for (p = mesh_grp_ent; p; p = (q = p)->next)
    ;

  if (q == NULL) {
    mesh_grp_ent = mesh;

  } else {
    q->next = mesh;
  }

  return 0;
}

static struct result_entry *make_result_entry(char *name_ID, int io,
                                              int fg_text, char *filename) {
  char *p;
  struct result_entry *result = NULL;
  result                      = HECMW_calloc(1, sizeof(*result));

  if (result == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  result->io      = io;
  result->fg_text = fg_text;
  result->next    = NULL;
  p               = HECMW_strdup(name_ID);

  if (p == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  result->name_ID = p;
  p               = HECMW_strdup(filename);

  if (p == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  result->filename = p;
  return result;
error:

  if (result) {
    HECMW_free(result->name_ID);
    HECMW_free(result->filename);
    HECMW_free(result);
  }

  return NULL;
}

static int add_result_entry(struct result_entry *result) {
  struct result_entry *p, *q;
  q = NULL;

  for (p = result_ent; p; p = (q = p)->next)
    ;

  if (q == NULL) {
    result_ent = result;

  } else {
    q->next = result;
  }

  return 0;
}

static struct restart_entry *make_restart_entry(char *name_ID, int io,
                                                char *filename) {
  char *p;
  struct restart_entry *restart = NULL;
  restart                       = HECMW_calloc(1, sizeof(*restart));

  if (restart == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  restart->io   = io;
  restart->next = NULL;
  p             = HECMW_strdup(name_ID);

  if (p == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  restart->name_ID = p;
  p                = HECMW_strdup(filename);

  if (p == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  restart->filename = p;
  return restart;
error:

  if (restart) {
    HECMW_free(restart->name_ID);
    HECMW_free(restart->filename);
    HECMW_free(restart);
  }

  return NULL;
}

static int add_restart_entry(struct restart_entry *restart) {
  struct restart_entry *p, *q;
  q = NULL;

  for (p = restart_ent; p; p = (q = p)->next)
    ;

  if (q == NULL) {
    restart_ent = restart;

  } else {
    q->next = restart;
  }

  return 0;
}

static struct ctrl_entry *make_ctrl_entry(char *name_ID, char *filename) {
  char *p;
  struct ctrl_entry *ctrl = NULL;
  ctrl                    = HECMW_calloc(1, sizeof(*ctrl));

  if (ctrl == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  ctrl->next = NULL;
  p          = HECMW_strdup(name_ID);

  if (p == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  ctrl->name_ID = p;
  p             = HECMW_strdup(filename);

  if (p == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  ctrl->filename = p;
  return ctrl;
error:

  if (ctrl) {
    HECMW_free(ctrl->name_ID);
    HECMW_free(ctrl->filename);
    HECMW_free(ctrl);
  }

  return NULL;
}

static int add_ctrl_entry(struct ctrl_entry *ctrl) {
  struct ctrl_entry *p, *q;
  q = NULL;

  for (p = ctrl_ent; p; p = (q = p)->next)
    ;

  if (q == NULL) {
    ctrl_ent = ctrl;

  } else {
    q->next = ctrl;
  }

  return 0;
}

/*----------------------------------------------------------------------------*/

static void do_logging(int loglv, int msgno, int add_location, const char *fmt,
                       va_list ap) {
  char line[100] = "";
  char msg[HECMW_MSG_LEN + 1];
  HECMW_vsnprintf(msg, sizeof(msg), fmt, ap);

  if (add_location) {
    char *s = "";

    if (strlen(msg) > 0) s = ": ";

    HECMW_snprintf(line, sizeof(line), "%s:%d%s", ctrl_filename,
                   HECMW_ctrllex_get_lineno(), s);
  }

  HECMW_set_error(msgno, "%s%s", line, msg);
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
    msg_no = HECMW_UTIL_E0003;

  } else {
    msg_no = msgno;
  }

  va_start(ap, fmt);
  do_logging(HECMW_LOG_ERROR, msg_no, 1, fmt, ap);
  va_end(ap);
}

/*----------------------------------------------------------------------------*/

static int read_mesh_header(void) {
  int token;
  /* !MESH */
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_H_MESH) {
    set_err(HECMW_UTIL_E0010, "!MESH required");
    return -1;
  }

  /* ',' */
  token = HECMW_ctrllex_next_token();

  if (token != ',') {
    set_err_token(token, HECMW_UTIL_E0010, "',' required after !MESH");
    return -1;
  }

  return 0;
}

static int read_mesh_head_param_name(char *name) {
  int token;
  char *p;
  token = HECMW_ctrllex_next_token();

  if (token != '=') {
    set_err_token(token, HECMW_UTIL_E0010, "'=' required after NAME");
    return -1;
  }

  /* NAME value */
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_NAME) {
    set_err_token(token, HECMW_UTIL_E0010,
                  "NAME must begin with a letter or '_'");
    return -1;
  }

  p = HECMW_ctrllex_get_text();

  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }

  strcpy(name, p);

  /* check */
  if (get_mesh_entry(name) || get_mesh_grp_entry(name)) {
    set_err(HECMW_UTIL_E0013, "");
    return -1;
  }

  return 0;
}

static int read_mesh_head_param_type(int *type) {
  int token;
  token = HECMW_ctrllex_next_token();

  if (token != '=') {
    set_err_token(token, HECMW_UTIL_E0010, "'=' required after TYPE");
    return -1;
  }

  /* TYPE value */
  token = HECMW_ctrllex_next_token();

  if (token == HECMW_CTRLLEX_K_HECMW_DIST) {
    *type = HECMW_CTRL_FTYPE_HECMW_DIST;

  } else if (token == HECMW_CTRLLEX_K_HECMW_ENTIRE) {
    *type = HECMW_CTRL_FTYPE_HECMW_ENTIRE;

  } else if (token == HECMW_CTRLLEX_K_GEOFEM) {
    *type = HECMW_CTRL_FTYPE_GEOFEM;

  } else if (token == HECMW_CTRLLEX_K_ABAQUS) {
    *type = HECMW_CTRL_FTYPE_ABAQUS;

  } else if (token == HECMW_CTRLLEX_K_NASTRAN) {
    *type = HECMW_CTRL_FTYPE_NASTRAN;
#if 0

  } else  if (token == HECMW_CTRLLEX_K_FEMAP) {
    *type = HECMW_CTRL_FTYPE_FEMAP;
#endif

  } else {
    set_err_token(token, HECMW_UTIL_E0010, "Invalid TYPE");
    return -1;
  }

  return 0;
}

static int read_mesh_head_param_io(int *io) {
  int token;
  token = HECMW_ctrllex_next_token();

  if (token != '=') {
    set_err_token(token, HECMW_UTIL_E0010, "'=' required after IO");
    return -1;
  }

  /* IO value */
  token = HECMW_ctrllex_next_token();

  if (token == HECMW_CTRLLEX_K_IN) {
    *io = HECMW_CTRL_FILE_IO_IN;

  } else if (token == HECMW_CTRLLEX_K_OUT) {
    *io = HECMW_CTRL_FILE_IO_OUT;

  } else {
    set_err_token(token, HECMW_UTIL_E0010, "Invalid IO");
    return -1;
  }

  return 0;
}

static int read_mesh_head_param_refine(int *refine) {
  int token;
  token = HECMW_ctrllex_next_token();

  if (token != '=') {
    set_err_token(token, HECMW_UTIL_E0010, "'=' required after REFINE");
    return -1;
  }

  /* REFINE value */
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_INT) {
    set_err_token(token, HECMW_UTIL_E0010, "Invalid REFINE");
    return -1;

  } else {
    *refine = HECMW_ctrllex_get_number();
  }

  return 0;
}

static int read_mesh_data(char *name, int type, int io, int refine) {
  int token;
  char *p;
  struct mesh_entry *mesh;
  /* filename */
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_NAME && token != HECMW_CTRLLEX_FILENAME) {
    set_err_token(token, HECMW_UTIL_E0010, "Invalid filename");
    return -1;
  }

  p = HECMW_ctrllex_get_text();

  if (strlen(p) > HECMW_FILENAME_LEN) {
    set_err(HECMW_IO_E0002, "");
    return -1;
  }

  /* create */
  mesh = make_mesh_entry(name, type, io, refine, p);

  if (mesh == NULL) return -1;

  /* add */
  if (add_mesh_entry(mesh)) return -1;

  /* NL*/
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_NL) {
    set_err_token(token, HECMW_UTIL_E0010, "NL required after filename");
    return -1;
  }

  return 0;
}

static int read_mesh(void) {
  int state;
  int token                     = -1;
  int flag_name                 = 0; /* flag for NAME */
  int flag_type                 = 0; /* flag for TYPE */
  int flag_io                   = 0; /* flag for IO */
  int flag_refine               = 0; /* flag for REFINE */
  int type                      = -1;
  int io                        = HECMW_CTRL_FILE_IO_IN;
  int refine                    = 0;
  char name[HECMW_NAME_LEN + 1] = "";
  enum {
    ST_FINISHED,
    ST_HEADER_LINE,
    ST_HEADER_LINE_PARAM,
    ST_DATA_LINE,
  };
  state = ST_HEADER_LINE;

  while (state != ST_FINISHED) {
    if (state == ST_HEADER_LINE) {
      if (read_mesh_header()) return -1;

      /* set next state */
      state = ST_HEADER_LINE_PARAM;

    } else if (state == ST_HEADER_LINE_PARAM) {
      token = HECMW_ctrllex_next_token();

      if (token == HECMW_CTRLLEX_K_NAME) {
        /* must */
        if (read_mesh_head_param_name(name)) return -1;

        flag_name = 1;

      } else if (token == HECMW_CTRLLEX_K_TYPE) {
        /* must */
        if (read_mesh_head_param_type(&type)) return -1;

        flag_type = 1;

      } else if (token == HECMW_CTRLLEX_K_IO) {
        /* optional */
        if (read_mesh_head_param_io(&io)) return -1;

        flag_io = 1;

      } else if (token == HECMW_CTRLLEX_K_REFINE) {
        /* optional */
        if (read_mesh_head_param_refine(&refine)) return -1;

        flag_refine = 1;

      } else {
        set_err_token(token, HECMW_UTIL_E0010, "Unknown parameter");
        return -1;
      }

      /* check next parameter */
      token = HECMW_ctrllex_next_token();

      if (token == HECMW_CTRLLEX_NL) {
        /* check NAME */
        if (!flag_name) {
          set_err(HECMW_UTIL_E0011, "");
          return -1;
        }

        /* check TYPE */
        if (!flag_type) {
          set_err(HECMW_UTIL_E0012, "");
          return -1;
        }

        state = ST_DATA_LINE;

      } else if (token == ',') {
        ; /* continue this state */

      } else {
        set_err_token(token, HECMW_UTIL_E0010, "Unknown parameter");
        return -1;
      }

    } else if (state == ST_DATA_LINE) {
      HECMW_assert(flag_name);
      HECMW_assert(flag_type);

      if (read_mesh_data(name, type, io, refine)) return -1;

      state = ST_FINISHED;

    } else {
      HECMW_assert(0);
    }
  }

  /* check */
  if (!strcmp(name, "fstrMSH") && type == HECMW_CTRL_FTYPE_HECMW_ENTIRE &&
      HECMW_comm_get_size() > 1) {
    set_err_token(token, HECMW_UTIL_E0010, "Invalid TYPE");
    return -1;
  }

  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_meshgrp_header(void) {
  int token;
  /* !MESH GROUP */
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_H_MESH_GROUP) {
    set_err(HECMW_UTIL_E0050, "!MESH GROUP required");
    return -1;
  }

  /* ',' */
  token = HECMW_ctrllex_next_token();

  if (token != ',') {
    set_err_token(token, HECMW_UTIL_E0050, "',' required after !MESH GROUP");
    return -1;
  }

  return 0;
}

static int read_meshgrp_head_param_name(char *name) {
  int token;
  char *p;
  token = HECMW_ctrllex_next_token();

  if (token != '=') {
    set_err_token(token, HECMW_UTIL_E0050, "'=' required after NAME");
    return -1;
  }

  /* NAME value */
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_NAME) {
    set_err_token(token, HECMW_UTIL_E0050,
                  "NAME must begin with a letter or '_'");
    return -1;
  }

  p = HECMW_ctrllex_get_text();

  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }

  strcpy(name, p);

  /* check */
  if (get_mesh_entry(name) || get_mesh_grp_entry(name)) {
    set_err(HECMW_UTIL_E0053, "");
    return -1;
  }

  return 0;
}

static int read_meshgrp_data(char *name) {
  int i, token, n_mesh, n_mesh_max;
  char *p, **q;
  char **mesh                    = NULL;
  struct mesh_grp_entry *meshgrp = NULL;
  struct mesh_entry *ment;
  n_mesh_max = 10; /* default */
  mesh       = HECMW_malloc(sizeof(*mesh) * n_mesh_max);

  if (mesh == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  n_mesh = 0;

  while (1) {
    /* filename */
    token = HECMW_ctrllex_next_token();

    if (token != HECMW_CTRLLEX_NAME && token != HECMW_CTRLLEX_FILENAME) {
      if (n_mesh == 0) {
        set_err_token(token, HECMW_UTIL_E0050, "name_ID required");
        return -1;
      }

      HECMW_ctrllex_unput_token();
      break;
    }

    p = HECMW_ctrllex_get_text();

    if (strlen(p) > HECMW_FILENAME_LEN) {
      set_err(HECMW_IO_E0002, "");
      return -1;
    }

    if ((ment = get_mesh_entry(p)) == NULL) {
      set_err_token(token, HECMW_UTIL_E0052, "name_ID: %s", p);
      return -1;
    }

    if (!is_entire_mesh(ment->type)) {
      set_err_token(token, HECMW_UTIL_E0055, "name_ID: %s", p);
      return -1;
    }

    if (n_mesh == n_mesh_max) {
      n_mesh_max *= 2;
      q = HECMW_realloc(mesh, sizeof(*mesh) * n_mesh_max);

      if (q == NULL) {
        HECMW_set_error(errno, "");
        return -1;
      }

      mesh = q;
    }

    mesh[n_mesh] = HECMW_strdup(p);

    if (mesh[n_mesh] == NULL) {
      HECMW_set_error(errno, "");
      return -1;
    }

    n_mesh++;
    /* ',' or NL*/
    token = HECMW_ctrllex_next_token();

    if (token != HECMW_CTRLLEX_NL && token != ',') {
      set_err_token(token, HECMW_UTIL_E0050,
                    "','  or NL required after name_ID");
      return -1;
    }
  }

  /* create */
  meshgrp = make_mesh_group_entry(name, n_mesh, mesh);

  if (meshgrp == NULL) return -1;

  /* add */
  if (add_mesh_group_entry(meshgrp)) return -1;

  for (i = 0; i < n_mesh; i++) {
    HECMW_free(mesh[i]);
  }

  HECMW_free(mesh);
  return 0;
}

static int read_meshgrp(void) {
  int token, state;
  int flag_name                 = 0; /* flag for NAME */
  char name[HECMW_NAME_LEN + 1] = "";
  enum {
    ST_FINISHED,
    ST_HEADER_LINE,
    ST_HEADER_LINE_PARAM,
    ST_DATA_LINE,
  };
  state = ST_HEADER_LINE;

  while (state != ST_FINISHED) {
    if (state == ST_HEADER_LINE) {
      if (read_meshgrp_header()) return -1;

      /* set next state */
      state = ST_HEADER_LINE_PARAM;

    } else if (state == ST_HEADER_LINE_PARAM) {
      token = HECMW_ctrllex_next_token();

      if (token == HECMW_CTRLLEX_K_NAME) {
        /* must */
        if (read_meshgrp_head_param_name(name)) return -1;

        flag_name = 1;

      } else {
        set_err_token(token, HECMW_UTIL_E0050, "Unknown parameter");
        return -1;
      }

      /* check next parameter */
      token = HECMW_ctrllex_next_token();

      if (token == HECMW_CTRLLEX_NL) {
        /* check NAME */
        if (!flag_name) {
          set_err(HECMW_UTIL_E0051, "");
          return -1;
        }

        state = ST_DATA_LINE;

      } else if (token == ',') {
        ; /* continue this state */

      } else {
        set_err_token(token, HECMW_UTIL_E0050, "Unknown parameter");
        return -1;
      }

    } else if (state == ST_DATA_LINE) {
      HECMW_assert(flag_name);

      if (read_meshgrp_data(name)) return -1;

      state = ST_FINISHED;

    } else {
      HECMW_assert(0);
    }
  }

  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_result_head(void) {
  int token;
  /* !RESULT */
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_H_RESULT) {
    set_err(HECMW_UTIL_E0020, "!RESULT required");
    return -1;
  }

  /* ',' */
  token = HECMW_ctrllex_next_token();

  if (token != ',') {
    set_err_token(token, HECMW_UTIL_E0020, "',' required after !RESULT");
    return -1;
  }

  return 0;
}

static int read_result_param_name(char *name) {
  int token;
  char *p;
  token = HECMW_ctrllex_next_token();

  if (token != '=') {
    set_err_token(token, HECMW_UTIL_E0020, "'=' required after NAME");
    return -1;
  }

  /* NAME value */
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_NAME) {
    set_err_token(token, HECMW_UTIL_E0020,
                  "NAME must begin with a letter or '_'");
    return -1;
  }

  p = HECMW_ctrllex_get_text();

  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }

  strcpy(name, p);

  /* check */
  if (get_result_entry(name)) {
    set_err(HECMW_UTIL_E0023, "");
    return -1;
  }

  return 0;
}

static int read_result_param_io(int *io) {
  int token;
  token = HECMW_ctrllex_next_token();

  if (token != '=') {
    set_err_token(token, HECMW_UTIL_E0020, "'=' required after IO");
    return -1;
  }

  /* IO value */
  token = HECMW_ctrllex_next_token();

  if (token == HECMW_CTRLLEX_K_IN) {
    *io = HECMW_CTRL_FILE_IO_IN;

  } else if (token == HECMW_CTRLLEX_K_OUT) {
    *io = HECMW_CTRL_FILE_IO_OUT;

  } else {
    set_err_token(token, HECMW_UTIL_E0020, "Invalid IO");
    return -1;
  }

  return 0;
}

static int read_result_param_type(int *fg_text) {
  int token;
  char s[HECMW_NAME_LEN + 1];
  char *p;
  char *sp;
  token = HECMW_ctrllex_next_token();

  if (token != '=') {
    set_err_token(token, HECMW_UTIL_E0020, "'=' required after TYPE");
    return -1;
  }

  /* TYPE value */
  token = HECMW_ctrllex_next_token();
  p     = HECMW_ctrllex_get_text();

  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_UTIL_E0020, "");
    return -1;
  }

  sp = s;

  while (*p) {
    *sp = (char)toupper(*p);
    p++;
    sp++;
  }

  *sp = 0;

  if (strcmp(s, "TEXT") == 0) {
    *fg_text = 1;

  } else if (strcmp(s, "BINARY") == 0) {
    *fg_text = 0;

  } else {
    set_err(HECMW_UTIL_E0020, "TEXT or BINARY required");
    return -1;
  }

  return 0;
}

static int read_result_data(char *name, int io, int fg_text) {
  int token;
  char *p;
  struct result_entry *result;
  /* filename */
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_NAME && token != HECMW_CTRLLEX_FILENAME) {
    set_err_token(token, HECMW_UTIL_E0020, "Invalid filename");
  }

  p = HECMW_ctrllex_get_text();

  if (strlen(p) > HECMW_FILENAME_LEN) {
    set_err(HECMW_IO_E0002, "NL required after filename");
    return -1;
  }

  /* create */
  result = make_result_entry(name, io, fg_text, p);

  if (result == NULL) return -1;

  /* add */
  if (add_result_entry(result)) return -1;

  /* NL*/
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_NL) {
    set_err_token(token, HECMW_UTIL_E0020, "NL required after filename");
    return -1;
  }

  return 0;
}

static int read_result(void) {
  int token, state;
  int flag_name = 0; /* flag for NAME */
  int flag_io   = 0; /* flag for IO */
  int io;
  int fg_text; /* default : text */
  char name[HECMW_NAME_LEN + 1] = "";
  enum {
    ST_FINISHED,
    ST_HEADER_LINE,
    ST_HEADER_LINE_PARAM,
    ST_DATA_LINE,
  };
  fg_text = 1;
  state   = ST_HEADER_LINE;

  while (state != ST_FINISHED) {
    if (state == ST_HEADER_LINE) {
      if (read_result_head()) return -1;

      state = ST_HEADER_LINE_PARAM;

    } else if (state == ST_HEADER_LINE_PARAM) {
      token = HECMW_ctrllex_next_token();

      if (token == HECMW_CTRLLEX_K_NAME) {
        /* must */
        if (read_result_param_name(name)) return -1;

        flag_name = 1;

      } else if (token == HECMW_CTRLLEX_K_IO) {
        /* must */
        if (read_result_param_io(&io)) return -1;

        flag_io = 1;

      } else if (token == HECMW_CTRLLEX_K_TYPE) {
        /* option */
        if (read_result_param_type(&fg_text)) return -1;

      } else {
        set_err_token(token, HECMW_UTIL_E0020, "Unknown parameter");
        return -1;
      }

      /* check next parameter */
      token = HECMW_ctrllex_next_token();

      if (token == HECMW_CTRLLEX_NL) {
        /* check NAME */
        if (!flag_name) {
          set_err(HECMW_UTIL_E0021, "");
          return -1;
        }

        /* check IO */
        if (!flag_io) {
          set_err(HECMW_UTIL_E0022, "");
          return -1;
        }

        state = ST_DATA_LINE;

      } else if (token == ',') {
        ; /* continue this state */

      } else {
        set_err_token(token, HECMW_UTIL_E0020, "Unknown parameter");
        return -1;
      }

    } else if (state == ST_DATA_LINE) {
      HECMW_assert(flag_name);

      if (read_result_data(name, io, fg_text)) return -1;

      state = ST_FINISHED;

    } else {
      HECMW_assert(0);
    }
  }

  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_restart_head(void) {
  int token;
  /* !RESULT */
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_H_RESTART) {
    set_err(HECMW_UTIL_E0040, "!RESTART required");
    return -1;
  }

  /* ',' */
  token = HECMW_ctrllex_next_token();

  if (token != ',') {
    set_err_token(token, HECMW_UTIL_E0040, "',' required after !RESTART");
    return -1;
  }

  return 0;
}

static int read_restart_param_name(char *name) {
  int token;
  char *p;
  token = HECMW_ctrllex_next_token();

  if (token != '=') {
    set_err_token(token, HECMW_UTIL_E0040, "'=' required after NAME");
    return -1;
  }

  /* NAME value */
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_NAME) {
    set_err_token(token, HECMW_UTIL_E0040,
                  "NAME must begin with a letter or '_'");
    return -1;
  }

  p = HECMW_ctrllex_get_text();

  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }

  strcpy(name, p);

  /* check */
  if (get_restart_entry(name)) {
    set_err(HECMW_UTIL_E0043, "");
    return -1;
  }

  return 0;
}

static int read_restart_param_io(int *io) {
  int token;
  token = HECMW_ctrllex_next_token();

  if (token != '=') {
    set_err_token(token, HECMW_UTIL_E0040, "'=' required after IO");
    return -1;
  }

  /* IO value */
  token = HECMW_ctrllex_next_token();

  if (token == HECMW_CTRLLEX_K_IN) {
    *io = HECMW_CTRL_FILE_IO_IN;

  } else if (token == HECMW_CTRLLEX_K_OUT) {
    *io = HECMW_CTRL_FILE_IO_OUT;

  } else if (token == HECMW_CTRLLEX_K_INOUT) {
    *io = HECMW_CTRL_FILE_IO_INOUT;

  } else {
    set_err_token(token, HECMW_UTIL_E0040, "Invalid IO");
    return -1;
  }

  return 0;
}

static int read_restart_data(char *name, int io) {
  int token;
  char *p;
  struct restart_entry *restart;
  /* filename */
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_NAME && token != HECMW_CTRLLEX_FILENAME) {
    set_err_token(token, HECMW_UTIL_E0040, "Invalid filename");
  }

  p = HECMW_ctrllex_get_text();

  if (strlen(p) > HECMW_FILENAME_LEN) {
    set_err(HECMW_IO_E0002, "NL required after filename");
    return -1;
  }

  /* create */
  restart = make_restart_entry(name, io, p);

  if (restart == NULL) return -1;

  /* add */
  if (add_restart_entry(restart)) return -1;

  /* NL*/
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_NL) {
    set_err_token(token, HECMW_UTIL_E0040, "NL required after filename");
    return -1;
  }

  return 0;
}

static int read_restart(void) {
  int token, state;
  int flag_name = 0; /* flag for NAME */
  int flag_io   = 0; /* flag for IO */
  int io;
  char name[HECMW_NAME_LEN + 1] = "";
  enum {
    ST_FINISHED,
    ST_HEADER_LINE,
    ST_HEADER_LINE_PARAM,
    ST_DATA_LINE,
  };
  state = ST_HEADER_LINE;

  while (state != ST_FINISHED) {
    if (state == ST_HEADER_LINE) {
      if (read_restart_head()) return -1;

      /* set next state */
      state = ST_HEADER_LINE_PARAM;

    } else if (state == ST_HEADER_LINE_PARAM) {
      token = HECMW_ctrllex_next_token();

      if (token == HECMW_CTRLLEX_K_NAME) {
        /* must */
        if (read_restart_param_name(name)) return -1;

        flag_name = 1;

      } else if (token == HECMW_CTRLLEX_K_IO) {
        /* must */
        if (read_restart_param_io(&io)) return -1;

        flag_io = 1;

      } else {
        set_err_token(token, HECMW_UTIL_E0040, "Unknown parameter");
        return -1;
      }

      /* check next parameter */
      token = HECMW_ctrllex_next_token();

      if (token == HECMW_CTRLLEX_NL) {
        /* check NAME */
        if (!flag_name) {
          set_err(HECMW_UTIL_E0041, "");
          return -1;
        }

        /* check IO */
        if (!flag_io) {
          set_err(HECMW_UTIL_E0042, "");
          return -1;
        }

        state = ST_DATA_LINE;

      } else if (token == ',') {
        ; /* continue this state */

      } else {
        set_err_token(token, HECMW_UTIL_E0040, "Unknown parameter");
        return -1;
      }

    } else if (state == ST_DATA_LINE) {
      HECMW_assert(flag_name);
      HECMW_assert(flag_io);

      if (read_restart_data(name, io)) return -1;

      state = ST_FINISHED;

    } else {
      HECMW_assert(0);
    }
  }

  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_control_head(void) {
  int token;
  /* !CONTROL */
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_H_CONTROL) {
    set_err(HECMW_UTIL_E0030, "!CONTROL required");
    return -1;
  }

  /* ',' */
  token = HECMW_ctrllex_next_token();

  if (token != ',') {
    set_err_token(token, HECMW_UTIL_E0030, "',' required after !CONTROL");
    return -1;
  }

  return 0;
}

static int read_control_head_param_name(char *name) {
  int token;
  char *p;
  token = HECMW_ctrllex_next_token();

  if (token != '=') {
    set_err_token(token, HECMW_UTIL_E0030, "'=' required after NAME");
    return -1;
  }

  /* NAME value */
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_NAME) {
    set_err_token(token, HECMW_UTIL_E0030,
                  "NAME must begin with a letter or '_'");
    return -1;
  }

  p = HECMW_ctrllex_get_text();

  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }

  strcpy(name, p);

  /* check */
  if (get_ctrl_entry(name)) {
    set_err(HECMW_UTIL_E0032, "");
    return -1;
  }

  return 0;
}

static int read_control_data(char *name) {
  int token;
  char *p;
  struct ctrl_entry *control;
  /* filename */
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_NAME && token != HECMW_CTRLLEX_FILENAME) {
    set_err_token(token, HECMW_UTIL_E0030, "Invalid filename");
    return -1;
  }

  p = HECMW_ctrllex_get_text();

  if (strlen(p) > HECMW_FILENAME_LEN) {
    set_err(HECMW_IO_E0002, "NL required after filename");
    return -1;
  }

  /* create */
  control = make_ctrl_entry(name, p);

  if (control == NULL) {
    return -1;
  }

  /* add */
  if (add_ctrl_entry(control)) {
    return -1;
  }

  /* NL*/
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_NL) {
    set_err_token(token, HECMW_UTIL_E0030, "NL required after filename");
    return -1;
  }

  return 0;
}

static int read_control(void) {
  int token, state;
  int flag_name                 = 0; /* flag for NAME */
  char name[HECMW_NAME_LEN + 1] = "";
  enum {
    ST_FINISHED,
    ST_HEADER_LINE,
    ST_HEADER_LINE_PARAM,
    ST_DATA_LINE,
  };
  state = ST_HEADER_LINE;

  while (state != ST_FINISHED) {
    if (state == ST_HEADER_LINE) {
      if (read_control_head()) return -1;

      /* set next state */
      state = ST_HEADER_LINE_PARAM;

    } else if (state == ST_HEADER_LINE_PARAM) {
      token = HECMW_ctrllex_next_token();

      if (token == HECMW_CTRLLEX_K_NAME) {
        /* must */
        if (read_control_head_param_name(name)) return -1;

        flag_name = 1;

      } else {
        set_err_token(token, HECMW_UTIL_E0030, "Unknown parameter");
        return -1;
      }

      /* check next parameter */
      token = HECMW_ctrllex_next_token();

      if (token == HECMW_CTRLLEX_NL) {
        /* check NAME */
        if (!flag_name) {
          set_err(HECMW_UTIL_E0031, "");
          return -1;
        }

        state = ST_DATA_LINE;

      } else if (token == ',') {
        ; /* continue this state */

      } else {
        set_err_token(token, HECMW_UTIL_E0030, "Unknown parameter");
        return -1;
      }

    } else if (state == ST_DATA_LINE) {
      HECMW_assert(flag_name);

      if (read_control_data(name)) return -1;

      state = ST_FINISHED;

    } else {
      HECMW_assert(0);
    }
  }

  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_subdir_head(void) {
  int token;
  /* !SUBDIR */
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_H_SUBDIR) {
    set_err(HECMW_UTIL_E0060, "!SUBDIR required");
    return -1;
  }

  /* ',' */
  token = HECMW_ctrllex_next_token();

  if (token != ',') {
    set_err_token(token, HECMW_UTIL_E0060, "',' required after !SUBDIR");
    return -1;
  }

  return 0;
}

static int read_subdir_head_param_limit(void) {
  int token;
  char *p;
  token = HECMW_ctrllex_next_token();

  if (token != '=') {
    set_err_token(token, HECMW_UTIL_E0060, "'=' required after LIMIT");
    return -1;
  }

  /* LIMIT value */
  token = HECMW_ctrllex_next_token();

  if (token != HECMW_CTRLLEX_INT) {
    set_err_token(token, HECMW_UTIL_E0060, "Invalid LIMIT");
    return -1;

  } else {
    nlimit = HECMW_ctrllex_get_number();
  }

  return 0;
}

static int read_subdir(void) {
  int token, state;
  int flag_name                 = 0; /* flag for NAME */
  char name[HECMW_NAME_LEN + 1] = "";
  enum {
    ST_FINISHED,
    ST_HEADER_LINE,
    ST_HEADER_LINE_PARAM,
  };
  nlimit = 5000;
  state  = ST_HEADER_LINE;

  while (state != ST_FINISHED) {
    if (state == ST_HEADER_LINE) {
      if (read_subdir_head()) return -1;

      /* set next state */
      state = ST_HEADER_LINE_PARAM;

    } else if (state == ST_HEADER_LINE_PARAM) {
      token = HECMW_ctrllex_next_token();

      if (token == HECMW_CTRLLEX_K_ON) {
        /* must */
        subdir_on = 1;
        flag_name = 1;

      } else if (token == HECMW_CTRLLEX_K_OFF) {
        /* must */
        subdir_on = 0;
        flag_name = 1;

      } else if (token == HECMW_CTRLLEX_K_LIMIT) {
        if (read_subdir_head_param_limit()) return -1;

      } else {
        set_err_token(token, HECMW_UTIL_E0060, "Unknown parameter");
        return -1;
      }

      /* check next parameter */
      token = HECMW_ctrllex_next_token();

      if (token == HECMW_CTRLLEX_NL) {
        if (!flag_name) {
          set_err(HECMW_UTIL_E0061, "");
          return -1;
        }

        state = ST_FINISHED;

      } else if (token == ',') {
        ; /* continue this state */

      } else {
        set_err_token(token, HECMW_UTIL_E0060, "Unknown parameter");
        return -1;
      }

    } else {
      HECMW_assert(0);
    }
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
    {HECMW_CTRLLEX_H_CONTROL, read_control},
    {HECMW_CTRLLEX_H_MESH, read_mesh},
    {HECMW_CTRLLEX_H_MESH_GROUP, read_meshgrp},
    {HECMW_CTRLLEX_H_RESULT, read_result},
    {HECMW_CTRLLEX_H_RESTART, read_restart},
    {HECMW_CTRLLEX_H_SUBDIR, read_subdir},
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

/*----------------------------------------------------------------------------*/

static int parse(void) {
  int token;
  ReadFunc func;

  while ((token = HECMW_ctrllex_next_token())) {
    if (token == HECMW_CTRLLEX_NL) continue;

    func = get_read_func(token);

    if (func == NULL) {
      char *p = HECMW_ctrllex_get_text();

      if (p[0] == '!') {
        set_err(HECMW_UTIL_E0004, "");

      } else {
        set_err(HECMW_UTIL_E0005, "");
      }

      return -1;
    }

    HECMW_ctrllex_unput_token(); /* unput !XXXX */

    if ((*func)()) return -1;
  }

  return 0;
}

int HECMW_ctrl_init_ex(const char *ctrlfile) {
  FILE *fp;
  HECMW_log(HECMW_LOG_DEBUG, "Getting control data");

  if (ctrlfile == NULL) {
    HECMW_set_error(HECMW_ALL_E0101, "Not specified control file name");
    return -1;
  }

  strcpy(ctrl_filename, ctrlfile);
  HECMW_log(HECMW_LOG_DEBUG, "Control file is '%s'", ctrl_filename);

  if ((fp = fopen(ctrl_filename, "r")) == NULL) {
    HECMW_set_error(HECMW_UTIL_E0001, "File: %s, %s", ctrl_filename,
                    strerror(errno));
    return -1;
  }

  if (HECMW_ctrllex_set_input(fp)) return -1;

  if (parse()) {
    return -1;
  }

  if (fclose(fp)) {
    HECMW_set_error(HECMW_UTIL_E0002, "File: %s, %s", ctrl_filename,
                    strerror(errno));
    return -1;
  }

  return 0;
}

int HECMW_ctrl_init(void) { return HECMW_ctrl_init_ex(HECMW_CTRL_FILE); }

int HECMW_ctrl_finalize(void) {
  HECMW_log(HECMW_LOG_DEBUG, "Finalizing control data");
  free_ctrl_entry();
  free_mesh_entry();
  free_mesh_grp_entry();
  free_result_entry();
  free_restart_entry();
  return 0;
}

void HECMW_ctrl_free_meshfiles(struct hecmw_ctrl_meshfiles *meshfiles) {
  int i;

  for (i = 0; i < meshfiles->n_mesh; i++) {
    HECMW_free(meshfiles->meshfiles[i].filename);
  }

  HECMW_free(meshfiles->meshfiles);
  HECMW_free(meshfiles);
}

static struct hecmw_ctrl_meshfiles *make_meshfiles_struct(
    int n_mesh, struct mesh_entry **mesh, int n_rank, int i_rank,
    int flag_rank_none) {
  int i, flag_rank, nrank, myrank, irank;
  char *fname;
  char *retval;
  struct hecmw_ctrl_meshfiles *files;
  char prefix[10];
  files = HECMW_malloc(sizeof(*files));

  if (files == NULL) {
    HECMW_set_error(errno, "");
    return NULL;
  }

  if (n_rank == 0) {
    nrank  = HECMW_comm_get_size();
    myrank = HECMW_comm_get_rank();

  } else {
    nrank  = n_rank;
    myrank = i_rank;
  }

  files->n_mesh    = n_mesh;
  files->meshfiles = HECMW_malloc(sizeof(*files->meshfiles) * n_mesh);

  if (files->meshfiles == NULL) {
    HECMW_set_error(errno, "");
    return NULL;
  }

  for (i = 0; i < n_mesh; i++) {
    struct hecmw_ctrl_meshfile *file = &files->meshfiles[i];
    struct mesh_entry *ment          = mesh[i];
    file->type                       = ment->type;
    file->io                         = ment->io;
    file->refine                     = ment->refine;

    if (flag_rank_none) {
      flag_rank = 0;

    } else {
      if (is_entire_mesh(file->type)) {
        flag_rank = 0;

      } else {
        flag_rank = 1;
      }
    }

    fname = HECMW_malloc(sizeof(char) * (HECMW_FILENAME_LEN + 1));

    if (fname == NULL) {
      HECMW_set_error(errno, "");
      return NULL;
    }

    if (ment->type == HECMW_CTRL_FTYPE_HECMW_ENTIRE) {
      retval = make_filename_r(NULL, NULL, NULL, ment->filename, "", myrank,
                               flag_rank, fname, HECMW_FILENAME_LEN);

    } else {
      if (subdir_on && nrank > nlimit) {
        irank = myrank / nlimit;
        sprintf(prefix, "TRUNK%d", irank);
        retval = make_filename_r("MESH", NULL, prefix, ment->filename, "", myrank,
                                 flag_rank, fname, HECMW_FILENAME_LEN);

      } else if (subdir_on && ment->type == HECMW_CTRL_FTYPE_HECMW_DIST) {
        retval = make_filename_r("MESH", NULL, NULL, ment->filename, "", myrank,
                                 flag_rank, fname, HECMW_FILENAME_LEN);

      } else if (subdir_on && nrank > 1) {
        retval = make_filename_r("MESH", NULL, NULL, ment->filename, "", myrank,
                                 flag_rank, fname, HECMW_FILENAME_LEN);

      } else {
        retval = make_filename_r(NULL, NULL, NULL, ment->filename, "", myrank,
                                 flag_rank, fname, HECMW_FILENAME_LEN);
      }
    }

    if (retval == NULL) {
      HECMW_set_error(HECMW_IO_E0002, "Cannot create mesh filename");
      HECMW_free(fname);
      return NULL;
    }

    file->filename = fname;
  }

  return files;
}

static struct hecmw_ctrl_meshfiles *get_meshfiles(char *name_ID, int n_rank,
                                                  int i_rank,
                                                  int flag_rank_none) {
  struct mesh_entry *mesh;
  struct mesh_grp_entry *mesh_grp;
  struct hecmw_ctrl_meshfiles *files = NULL;
  mesh_grp                           = get_mesh_grp_entry(name_ID);

  if (mesh_grp) {
    files = make_meshfiles_struct(mesh_grp->n_mesh, mesh_grp->mesh, n_rank,
                                  i_rank, flag_rank_none);

    if (files == NULL) return NULL;
  }

  if (files == NULL) {
    mesh = get_mesh_entry(name_ID);

    if (mesh) {
      files = make_meshfiles_struct(1, &mesh, n_rank, i_rank, flag_rank_none);

      if (files == NULL) return NULL;
    }
  }

  if (files == NULL) {
    HECMW_set_error(HECMW_UTIL_E0014, "NAME: %s", name_ID);
    return NULL;
  }

  return files;
}

struct hecmw_ctrl_meshfiles *HECMW_ctrl_get_meshfiles(char *name_ID) {
  return get_meshfiles(name_ID, 0, 0, 0);
}

struct hecmw_ctrl_meshfiles *HECMW_ctrl_get_meshfiles_header(char *name_ID) {
  return get_meshfiles(name_ID, 0, 0, 1);
}

struct hecmw_ctrl_meshfiles *HECMW_ctrl_get_meshfiles_sub(char *name_ID,
                                                          int n_rank,
                                                          int i_rank) {
  return get_meshfiles(name_ID, n_rank, i_rank, 0);
}

struct hecmw_ctrl_meshfiles *HECMW_ctrl_get_meshfiles_header_sub(char *name_ID,
                                                                 int n_rank,
                                                                 int i_rank) {
  return get_meshfiles(name_ID, n_rank, i_rank, 1);
}

static char *get_result_file(char *name_ID, int istep, int n_rank,
                             int i_rank, int *fg_text, int flag_rank_none) {
  int nrank, myrank, irank;
  char *fname, *retfname;
  struct result_entry *result;
  char subname[10], prefix[10];
  result = get_result_entry(name_ID);

  if (result == NULL) {
    HECMW_set_error(HECMW_UTIL_E0024, "NAME: %s",
                    name_ID ? name_ID : "Not specified");
    return NULL;
  }

  if (n_rank == 0) {
    nrank  = HECMW_comm_get_size();
    myrank = HECMW_comm_get_rank();

  } else {
    nrank  = n_rank;
    myrank = i_rank;
  }

  if (subdir_on && !strcmp(name_ID, "vis_out")) {
    fname = make_filename(name_ID, NULL, NULL, result->filename, "", myrank,
                          flag_rank_none);

  } else if (subdir_on && nrank > nlimit) {
    sprintf(subname, "STEP%d", istep);
    irank = myrank / nlimit;
    sprintf(prefix, "TRUNK%d", irank);
    fname = make_filename(name_ID, subname, prefix, result->filename, "",
                          myrank, flag_rank_none);

  } else if (subdir_on) {
    sprintf(subname, "STEP%d", istep);
    fname = make_filename(name_ID, subname, NULL, result->filename, "", myrank,
                          flag_rank_none);

  } else {
    fname = make_filename(NULL, NULL, NULL, result->filename, "", myrank,
                          flag_rank_none);
  }

  if (fname == NULL) {
    HECMW_set_error(HECMW_IO_E0002, "Cannot create result filename");
    return NULL;
  }

  retfname = HECMW_strdup(fname);

  if (retfname == NULL) {
    HECMW_set_error(errno, "");
    return NULL;
  }

  *fg_text = result->fg_text;
  return retfname;
}

static char *get_filename_body(char *file) {
  static char filename[HECMW_FILENAME_LEN+1];
  strcpy(filename, "");
  strcat(filename, file);
  return filename;
}

static char *get_result_filebody(char *name_ID) {
  struct result_entry *result;
  char *fname, *retfname;

  result = get_result_entry(name_ID);

  if (result == NULL) {
    HECMW_set_error(HECMW_UTIL_E0024, "NAME: %s",
                    name_ID ? name_ID : "Not specified");
    return NULL;
  }

  fname = get_filename_body(result->filename);

  retfname = HECMW_strdup(fname);

  if (retfname == NULL) {
    HECMW_set_error(errno, "");
    return NULL;
  }

  return retfname;
}

char *HECMW_ctrl_get_result_file(char *name_ID, int istep,
                                 int *fg_text) {
  return get_result_file(name_ID, istep, 0, 0, fg_text, 1);
}

char *HECMW_ctrl_get_result_fileheader(char *name_ID, int istep,
                                       int *fg_text) {
  return get_result_file(name_ID, istep, 0, 0, fg_text, 0);
}

char *HECMW_ctrl_get_result_file_sub(char *name_ID, int istep,
                                     int n_rank, int i_rank, int *fg_text) {
  return get_result_file(name_ID, istep, n_rank, i_rank, fg_text, 1);
}

char *HECMW_ctrl_get_result_fileheader_sub(char *name_ID, int istep,
                                           int n_rank, int i_rank,
                                           int *fg_text) {
  return get_result_file(name_ID, istep, n_rank, i_rank, fg_text, 0);
}

char *HECMW_ctrl_get_result_filebody(char *name_ID) {
  return get_result_filebody(name_ID);
}

char *HECMW_ctrl_get_restart_file(char *name_ID) {
  int nrank, myrank, irank;
  char *fname, *retfname;
  struct restart_entry *restart;
  char prefix[10];
  restart = get_restart_entry(name_ID);

  if (restart == NULL) {
    HECMW_set_error(HECMW_UTIL_E0044, "NAME: %s",
                    name_ID ? name_ID : "Not specified");
    return NULL;
  }

  nrank  = HECMW_comm_get_size();
  myrank = HECMW_comm_get_rank();

  if (subdir_on && nrank > nlimit) {
    irank = myrank / nlimit;
    sprintf(prefix, "TRUNK%d", irank);
    fname =
        make_filename(name_ID, NULL, prefix, restart->filename, "", myrank, 1);

  } else if (subdir_on) {
    fname =
        make_filename(name_ID, NULL, NULL, restart->filename, "", myrank, 1);

  } else {
    fname = make_filename(NULL, NULL, NULL, restart->filename, "", myrank, 1);
  }

  if (fname == NULL) {
    HECMW_set_error(HECMW_IO_E0002, "Cannot create restart filename");
    return NULL;
  }

  retfname = HECMW_strdup(fname);

  if (retfname == NULL) {
    HECMW_set_error(errno, "");
    return NULL;
  }

  return retfname;
}

char *HECMW_ctrl_get_restart_file_by_io(int io) {
  int nrank, myrank, irank;
  char *fname, *retfname;
  struct restart_entry *p;
  char prefix[10];
  p = get_restart_entry_by_io(io);

  if (p == NULL) {
    HECMW_set_error(HECMW_UTIL_E0045, "");
    return NULL;
  }

  nrank  = HECMW_comm_get_size();
  myrank = HECMW_comm_get_rank();

  if (subdir_on && nrank > nlimit) {
    irank = myrank / nlimit;
    sprintf(prefix, "TRUNK%d", irank);
    fname = make_filename(p->name_ID, NULL, prefix, p->filename, "", myrank, 1);

  } else if (subdir_on) {
    fname = make_filename(p->name_ID, NULL, NULL, p->filename, "", myrank, 1);

  } else {
    fname = make_filename(NULL, NULL, NULL, p->filename, "", myrank, 1);
  }

  if (fname == NULL) {
    HECMW_set_error(HECMW_IO_E0002, "Cannot create restart filename");
    return NULL;
  }

  retfname = HECMW_strdup(fname);

  if (retfname == NULL) {
    HECMW_set_error(errno, "");
    return NULL;
  }

  return retfname;
}

char *HECMW_ctrl_get_control_file(char *name_ID) {
  char *fname;
  struct ctrl_entry *ctrl;
  ctrl = get_ctrl_entry(name_ID);

  if (ctrl == NULL) {
    HECMW_set_error(HECMW_UTIL_E0033, "NAME: %s",
                    name_ID ? name_ID : "Not specified");
    return NULL;
  }

  fname = HECMW_strdup(ctrl->filename);
  return fname;
}

int HECMW_ctrl_is_exists_control(char *name_ID) {
  return get_ctrl_entry(name_ID) ? 1 : 0;
}

int HECMW_ctrl_make_subdir(char *filename) {
  char fname[HECMW_FILENAME_LEN + 1];
  char dirname[HECMW_FILENAME_LEN + 1];
  char *token;
  char separator[10];
  char *saveptr;
  mode_t mode;
  DIR *dp;
#ifndef _WINDOWS
  mode = S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
#endif
  strcpy(fname, filename);
  sprintf(separator, "%c", HECMW_get_path_separator());
#if defined(__WIN32__) || defined(__WIN64__)
  /* strtok is thread-safe on Windows */
  token = strtok(fname, separator);
  sprintf(dirname, "%s", token);
  token = strtok(NULL, separator);
#else
  token = strtok_r(fname, separator, &saveptr);
  sprintf(dirname, "%s", token);
  token = strtok_r(NULL, separator, &saveptr);
#endif

  while (token) {
    if ((dp = opendir(dirname)) == NULL) {
#ifndef _WINDOWS

      if (mkdir(dirname, mode) != 0) {
#else

      if (mkdir(dirname) != 0) {
#endif

        if (errno != EEXIST) return -1;
      }

    } else {
      closedir(dp);
    }

    strcat(dirname, separator);
    strcat(dirname, token);
#if defined(__WIN32__) || defined(__WIN64__)
    token = strtok(NULL, separator);
#else
    token = strtok_r(NULL, separator, &saveptr);
#endif
  }

  return 0;
}

int HECMW_ctrl_is_subdir(void) { return subdir_on; }

/*---------------------------------------------------------------------------*/

void hecmw_ctrl_init_if(int *err) {
  *err = 1;

  if (HECMW_ctrl_init()) return;

  *err = 0;
}

void hecmw_ctrl_init_if_(int *err) { hecmw_ctrl_init_if(err); }

void hecmw_ctrl_init_if__(int *err) { hecmw_ctrl_init_if(err); }

void HECMW_CTRL_INIT_IF(int *err) { hecmw_ctrl_init_if(err); }

/*---------------------------------------------------------------------------*/

void hecmw_ctrl_init_ex_if(char *ctrlfile, int *err, int len) {
  char c_ctrlfile[HECMW_FILENAME_LEN + 1];
  *err = 1;

  if (HECMW_strcpy_f2c_r(ctrlfile, len, c_ctrlfile, sizeof(c_ctrlfile)) == NULL)
    return;

  if (HECMW_ctrl_init_ex(c_ctrlfile)) return;

  *err = 0;
}

void hecmw_ctrl_init_ex_if_(char *ctrlfile, int *err, int len) {
  hecmw_ctrl_init_ex_if(ctrlfile, err, len);
}

void hecmw_ctrl_init_ex_if__(char *ctrlfile, int *err, int len) {
  hecmw_ctrl_init_ex_if(ctrlfile, err, len);
}

void HECMW_CTRL_INIT_EX_IF(char *ctrlfile, int *err, int len) {
  hecmw_ctrl_init_ex_if(ctrlfile, err, len);
}

/*---------------------------------------------------------------------------*/

void hecmw_ctrl_finalize_if(void) { HECMW_ctrl_finalize(); }

void hecmw_ctrl_finalize_if_(void) { hecmw_ctrl_finalize_if(); }

void hecmw_ctrl_finalize_if__(void) { hecmw_ctrl_finalize_if(); }

void HECMW_CTRL_FINALIZE_IF(void) { hecmw_ctrl_finalize_if(); }

/*---------------------------------------------------------------------------*/

void hecmw_ctrl_get_control_file_if(char *name_ID, char *buf, int *err,
                                    int name_len, int buf_len) {
  char c_name_ID[HECMW_NAME_LEN + 1];
  char *c_buf;
  int ret;
  *err = 1;

  if (HECMW_strcpy_f2c_r(name_ID, name_len, c_name_ID, sizeof(c_name_ID)) ==
      NULL)
    return;

  if ((c_buf = HECMW_ctrl_get_control_file(c_name_ID)) == NULL) return;

  ret = HECMW_strcpy_c2f(c_buf, buf, buf_len);
  HECMW_free(c_buf);

  if (ret == 0) return;

  *err = 0;
}

void hecmw_ctrl_get_control_file_if_(char *name_ID, char *buf, int *err,
                                     int name_len, int buf_len) {
  hecmw_ctrl_get_control_file_if(name_ID, buf, err, name_len, buf_len);
}

void hecmw_ctrl_get_control_file_if__(char *name_ID, char *buf, int *err,
                                      int name_len, int buf_len) {
  hecmw_ctrl_get_control_file_if(name_ID, buf, err, name_len, buf_len);
}

void HECMW_CTRL_GET_CONTROL_FILE_IF(char *name_ID, char *buf, int *err,
                                    int name_len, int buf_len) {
  hecmw_ctrl_get_control_file_if(name_ID, buf, err, name_len, buf_len);
}

/*---------------------------------------------------------------------------*/

void hecmw_ctrl_make_subdir(char *filename, int *err, int len) {
  char fname[HECMW_FILENAME_LEN + 1];
  *err = 1;

  if (HECMW_strcpy_f2c_r(filename, len, fname, sizeof(fname)) == NULL) return;

  if (HECMW_ctrl_make_subdir(fname) != 0) return;

  *err = 0;
}

void hecmw_ctrl_make_subdir_(char *filename, int *err, int len) {
  hecmw_ctrl_make_subdir(filename, err, len);
}

void hecmw_ctrl_make_subdir__(char *filename, int *err, int len) {
  hecmw_ctrl_make_subdir(filename, err, len);
}

void HECMW_CTRL_MAKE_SUBDIR(char *filename, int *err, int len) {
  hecmw_ctrl_make_subdir(filename, err, len);
}

/*---------------------------------------------------------------------------*/

void hecmw_ctrl_is_subdir(int *flag, int *limit) {
  *flag  = subdir_on;
  *limit = nlimit;
}

void hecmw_ctrl_is_subdir_(int *flag, int *limit) {
  hecmw_ctrl_is_subdir(flag, limit);
}

void hecmw_ctrl_is_subdir__(int *flag, int *limit) {
  hecmw_ctrl_is_subdir(flag, limit);
}

void HECMW_CTRL_IS_SUBDIR(int *flag, int *limit) {
  hecmw_ctrl_is_subdir(flag, limit);
}
