/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdbool.h>
#include "hecmw_util.h"
#include "hecmw_heclex.h"
#include "hecmw_io_hec.h"
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

static char grid_filename[HECMW_FILENAME_LEN + 1]    = "Unknown";
static char include_filename[HECMW_FILENAME_LEN + 1] = "Unknown";

static int connectivity_type = HECMW_CONNTYPE_HECMW;

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
    p = HECMW_heclex_is_including() ? include_filename : grid_filename;
    HECMW_snprintf(line, sizeof(line), "%s:%d%s", p, HECMW_heclex_get_lineno(),
                   s);
  }
  if (loglv == HECMW_LOG_ERROR) {
    HECMW_set_error(msgno, "%s%s", line, msg);
  } else {
    HECMW_print_msg(loglv, msgno, "%s%s", line, msg);
  }
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
    msg_no = HECMW_IO_HEC_E0003;
  } else {
    msg_no = msgno;
  }
  va_start(ap, fmt);
  do_logging(HECMW_LOG_ERROR, msg_no, 1, fmt, ap);
  va_end(ap);
}

/*-----------------------------------------------------------------------------
  ReadFunc
*/

static int read_input(int msgno_invalid_token) {
  int token;
  char *p;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, msgno_invalid_token, "'=' required after INPUT");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_FILENAME && token != HECMW_HECLEX_NAME) {
    set_err_token(token, msgno_invalid_token, "Invalid filename for INPUT");
    return -1;
  }
  p = HECMW_heclex_get_text();
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

static int read_amp_head(void) {
  int token;

  /* !AMPLITUDE */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_H_AMPLITUDE) {
    set_err_token(token, HECMW_IO_HEC_E0100, "!AMPLITUDE required");
    return -1;
  }

  token = HECMW_heclex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_HEC_E0100, "',' required after !AMPLITUDE");
    return -1;
  }

  return 0;
}

static int read_amp_param_name(char *name) {
  int token;
  char *p;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E0100, "'=' required after NAME");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_NAME) {
    set_err_token(token, HECMW_IO_HEC_E0100,
                  "NAME must begin with a letter or '_'");
    return -1;
  }
  p = HECMW_heclex_get_text();
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

static int read_amp_param_type(int *type) {
  int token;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E0100, "'=' required after TYPE");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_K_TIMEVALUE) {
    set_err_token(token, HECMW_IO_HEC_E0100, "Invalid TYPE");
    return -1;
  }
  *type = HECMW_HECLEX_K_TIMEVALUE;
  return 0;
}

static int read_amp_param_definition(int *definition) {
  int token;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E0100, "'=' required after DEFINITION");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_K_TABULAR) {
    set_err_token(token, HECMW_IO_HEC_E0100, "Invalid DEFINITION");
    return -1;
  }
  *definition = HECMW_AMP_TYPEDEF_TABULAR;
  return 0;
}

static int read_amp_param_time(int *time) {
  int token;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E0100, "'=' after TIME required");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_K_STEP_TIME) {
    set_err_token(token, HECMW_IO_HEC_E0100, "Invalid TIME");
    return -1;
  }
  *time = HECMW_AMP_TYPETIME_STEP;
  return 0;
}

static int read_amp_param_value(int *value) {
  int token;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E0100, "'=' required after VALUE");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token == HECMW_HECLEX_K_RELATIVE) {
    *value = HECMW_AMP_TYPEVAL_RELATIVE;
  } else if (token == HECMW_HECLEX_K_ABSOLUTE) {
    *value = HECMW_AMP_TYPEVAL_ABSOLUTE;
  } else {
    set_err_token(token, HECMW_IO_HEC_E0100, "Invalid VALUE");
    return -1;
  }
  return 0;
}

static int read_amp_data(char *name, int type, int definition, int time,
                         int value) {
  int i, token;
  const int NITEM = 4;

  i = 0;
  while (1) {
    double val, t, tmp;

    token = HECMW_heclex_next_token();
    if (i != 0 && token == HECMW_HECLEX_NL) break;
    /* VAL */
    if (token != HECMW_HECLEX_DOUBLE && token != HECMW_HECLEX_INT) {
      set_err_token(token, HECMW_IO_HEC_E0100, "VAL required");
      return -1;
    }
    val = HECMW_heclex_get_number();

    /* ',' */
    token = HECMW_heclex_next_token();
    if (token != ',') {
      set_err_token(token, HECMW_IO_HEC_E0100, "',' required after VAL");
      return -1;
    }

    /* T */
    token = HECMW_heclex_next_token();
    if (token != HECMW_HECLEX_DOUBLE && token != HECMW_HECLEX_INT) {
      set_err_token(token, HECMW_IO_HEC_E0100, "T required");
      return -1;
    }
    t = HECMW_heclex_get_number();

    /* type ABAQUS*/
    if (type == HECMW_HECLEX_K_TIMEVALUE) {
      tmp = val;
      val = t;
      t   = tmp;
    }

    /* add */
    if (HECMW_io_add_amp(name, definition, time, value, val, t) == NULL) {
      return -1;
    }

    i++;

    /* ',' or NL */
    token = HECMW_heclex_next_token();
    if (token != ',' && token != HECMW_HECLEX_NL) {
      set_err_token(token, HECMW_IO_HEC_E0100, "',' or NL required");
      return -1;
    }
    if (token == ',' && i == NITEM) {
      token = HECMW_heclex_next_token();
      if (token != HECMW_HECLEX_NL) {
        set_err_token(token, HECMW_IO_HEC_E0100, "Only %d items allow per line",
                      NITEM);
        return -1;
      }
      break;
    }
    if (token == HECMW_HECLEX_NL) break;
  }
  return 0;
}

static int read_amplitude(void) {
  int token, state;
  int type                      = -1;
  int definition                = HECMW_AMP_TYPEDEF_TABULAR;
  int time                      = HECMW_AMP_TYPETIME_STEP;
  int value                     = HECMW_AMP_TYPEVAL_RELATIVE;
  int flag_name                 = 0; /* flag for NAME */
  int flag_type                 = 0; /* flag for TYPE */
  int flag_definition           = 0; /* flag for DEFINITION */
  int flag_time                 = 0; /* flag for TIME */
  int flag_value                = 0; /* flag for VALUE */
  int flag_input                = 0; /* flag for INPUT */
  char name[HECMW_NAME_LEN + 1] = "";
  enum {
    ST_FINISHED,
    ST_HEADER_LINE,
    ST_HEADER_LINE_PARAM,
    ST_DATA_INCLUDE,
    ST_DATA_LINE,
    ST_FINALIZE
  };

  state = ST_HEADER_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_HEADER_LINE) {
      if (read_amp_head()) return -1;
      state = ST_HEADER_LINE_PARAM;
    } else if (state == ST_HEADER_LINE_PARAM) {
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_K_NAME) {
        /* must */
        if (read_amp_param_name(name)) return -1;
        flag_name = 1;
      } else if (token == HECMW_HECLEX_K_TYPE) {
        /* optional */
        if (read_amp_param_type(&type)) return -1;
        flag_type = 1;
      } else if (token == HECMW_HECLEX_K_DEFINITION) {
        /* optional */
        if (read_amp_param_definition(&definition)) return -1;
        flag_definition = 1;
      } else if (token == HECMW_HECLEX_K_TIME) {
        /* optional */
        if (read_amp_param_time(&time)) return -1;
        flag_time = 1;
      } else if (token == HECMW_HECLEX_K_VALUE) {
        /* optional */
        if (read_amp_param_value(&value)) return -1;
        flag_value = 1;
      } else if (token == HECMW_HECLEX_K_INPUT) {
        /* optional */
        if (read_input(HECMW_IO_HEC_E0100)) return -1;
        flag_input = 1;
      } else {
        set_err_token(token, HECMW_IO_HEC_E0100, "Unknown parameter");
        return -1;
      }

      /* check next state */
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_NL) {
        /* check NAME */
        if (!flag_name) {
          set_err(HECMW_IO_HEC_E0101, "");
          return -1;
        }
        state = flag_input ? ST_DATA_INCLUDE : ST_DATA_LINE;
      } else if (token == ',') {
        ; /* continue this state */
      } else {
        set_err_token(token, HECMW_IO_HEC_E0100, "Unknown parameter");
        return -1;
      }
    } else if (state == ST_DATA_INCLUDE) {
      HECMW_assert(flag_input);
      HECMW_assert(flag_name);

      if (HECMW_heclex_switch_to_include(include_filename)) return -1;
      state = ST_DATA_LINE;
    } else if (state == ST_DATA_LINE) {
      HECMW_assert(flag_name);
      if (read_amp_data(name, type, definition, time, value)) return -1;

      /* check next state */
      token = HECMW_heclex_next_token();
      if (token != HECMW_HECLEX_DOUBLE && token != HECMW_HECLEX_INT) {
        state = ST_FINISHED;
      }
      HECMW_heclex_unput_token();
    } else {
      HECMW_assert(0);
    }
  }
  HECMW_log(HECMW_LOG_DEBUG, "read_amplitude done");
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_contact_or_insert_pair_head(int *last_token, int *type) {
  int token;

  /* !CONTACT */
  token = HECMW_heclex_next_token();
  if ( token == HECMW_HECLEX_H_CONTACT_PAIR ){
    *type = HECMW_CONTACT_TYPE_NODE_SURF;
  } else if ( token == HECMW_HECLEX_H_EMBED_PAIR ){
    *type = HECMW_CONTACT_TYPE_NODE_ELEM;
  } else {
    set_err_token(token, HECMW_IO_HEC_E2100, "!CONTACT or !EMBED required");
    return -1;
  }

  token = HECMW_heclex_next_token();
  if (token != ',' && token != HECMW_HECLEX_NL) {
    set_err_token(token, HECMW_IO_HEC_E2100,
                  "',' or NL required after !CONTACT");
    return -1;
  }
  *last_token = token;

  return 0;
}

static int read_contact_pair_param_name(char *name) {
  int token;
  char *p;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E2100, "'=' required after NAME");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_NAME) {
    set_err_token(token, HECMW_IO_HEC_E2100,
                  "NAME must begin with a letter or '_'");
    return -1;
  }
  p = HECMW_heclex_get_text();
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

static int read_contact_pair_param_type(int *type) {
  int token;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E2100, "'=' required after TYPE");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token == HECMW_HECLEX_K_NODE_SURF) {
    *type = HECMW_CONTACT_TYPE_NODE_SURF;
  } else if (token == HECMW_HECLEX_K_SURF_SURF) {
    *type = HECMW_CONTACT_TYPE_SURF_SURF;
  } else if (token == HECMW_HECLEX_K_NODE_ELEM) {
    *type = HECMW_CONTACT_TYPE_NODE_ELEM;
  } else {
    set_err_token(token, HECMW_IO_HEC_E2100, "Invalid  TYPE");
    return -1;
  }
  return 0;
}

static int read_contact_pair_data(char *name, int type) {
  int token;
  char *slave_grp, *master_grp;

  slave_grp  = NULL;
  master_grp = NULL;

  /* SLAVE NODE/SURF GROUP */
  token = HECMW_heclex_next_token();
  if (token == HECMW_HECLEX_NAME) {
    slave_grp = HECMW_heclex_get_text();
    if (strlen(slave_grp) > HECMW_NAME_LEN) {
      set_err(HECMW_IO_E0001, "");
      return -1;
    }
    HECMW_toupper(slave_grp);
    slave_grp = HECMW_strdup(slave_grp);
    if (slave_grp == NULL) {
      HECMW_set_error(errno, "");
      return -1;
    }
  } else {
    if (type == HECMW_CONTACT_TYPE_NODE_SURF) {
      set_err_token(token, HECMW_IO_HEC_E2100, "NGROUP name required");
    } else if (type == HECMW_CONTACT_TYPE_SURF_SURF) {
      set_err_token(token, HECMW_IO_HEC_E2100, "SGROUP name required");
    } else if (type == HECMW_CONTACT_TYPE_NODE_ELEM) {
      set_err_token(token, HECMW_IO_HEC_E2100, "NGROUP name required");
    } else {
      HECMW_assert(0);
    }
    return -1;
  }

  /* ',' */
  token = HECMW_heclex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_HEC_E2100, "',' required after SGROUP");
    return -1;
  }

  /* MASTER SURF GROUP */
  token = HECMW_heclex_next_token();
  if (token == HECMW_HECLEX_NAME) {
    master_grp = HECMW_heclex_get_text();
    if (strlen(master_grp) > HECMW_NAME_LEN) {
      set_err(HECMW_IO_E0001, "");
      return -1;
    }
    HECMW_toupper(master_grp);
    master_grp = HECMW_strdup(master_grp);
    if (master_grp == NULL) {
      HECMW_set_error(errno, "");
      return -1;
    }
  } else {
    set_err_token(token, HECMW_IO_HEC_E2100, "SGROUP name required");
    return -1;
  }

  /* NL */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_NL) {
    set_err_token(token, HECMW_IO_HEC_E1000, "NL required after NGROUP");
    return -1;
  }

  /* add */
  if (HECMW_io_add_contact(name, type, slave_grp, master_grp) == NULL) {
    return -1;
  };
  HECMW_free(slave_grp);
  HECMW_free(master_grp);

  return 0;
}

static int read_contact_pair(void) {
  int token, state;
  int type                      = HECMW_CONTACT_TYPE_NODE_SURF;
  int flag_name                 = 0; /* flag for NAME */
  int flag_input                = 0; /* flag for INPUT */
  char name[HECMW_NAME_LEN + 1] = "";
  enum {
    ST_FINISHED,
    ST_HEADER_LINE,
    ST_HEADER_LINE_PARAM,
    ST_DATA_INCLUDE,
    ST_DATA_LINE,
  };

  state = ST_HEADER_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_HEADER_LINE) {
      if (read_contact_or_insert_pair_head(&token,&type)) return -1;
      if (token == HECMW_HECLEX_NL) {
        state = ST_DATA_LINE;
      } else if (token == ',') {
        state = ST_HEADER_LINE_PARAM;
      } else {
        HECMW_assert(0);
      }
    } else if (state == ST_HEADER_LINE_PARAM) {
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_K_NAME) {
        /* must */
        if (read_contact_pair_param_name(name)) return -1;
        flag_name = 1;
      } else if (token == HECMW_HECLEX_K_TYPE) {
        /* optional */
        if (read_contact_pair_param_type(&type)) return -1;
      } else if (token == HECMW_HECLEX_K_INPUT) {
        /* oprtional */
        if (read_input(HECMW_IO_HEC_E2100)) return -1;
        flag_input = 1;
      } else {
        set_err_token(token, HECMW_IO_HEC_E2100, "Unknown parameter");
        return -1;
      }

      /* check next parameter */
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_NL) {
        /* check NAME */
        if (!flag_name) {
          set_err(HECMW_IO_HEC_E2101, "");
          return -1;
        }
        if (flag_input) {
          state = ST_DATA_INCLUDE;
        } else {
          state = ST_DATA_LINE;
        }
      } else if (token == ',') {
        ; /* continue this state */
      } else {
        set_err_token(token, HECMW_IO_HEC_E2100, "Unknown parameter");
        return -1;
      }
    } else if (state == ST_DATA_INCLUDE) {
      HECMW_assert(flag_input);

      if (HECMW_heclex_switch_to_include(include_filename)) return -1;
      state = ST_DATA_LINE;
    } else if (state == ST_DATA_LINE) {
      if (read_contact_pair_data(name, type)) return -1;

      /* check next state */
      token = HECMW_heclex_next_token();
      if (token != HECMW_HECLEX_NAME) {
        state = ST_FINISHED;
      }
      HECMW_heclex_unput_token();
    } else {
      HECMW_assert(0);
    }
  }
  HECMW_log(HECMW_LOG_DEBUG, "read_contact done");
  return 0;
}

/*----------------------------------------------------------------------------*/
#if 0
static int
read_ecopy(void)
{
	fprintf(stderr, "!ECOPY has not implemented yet\n");
	HECMW_abort(HECMW_comm_get_comm());
	return 0;
}


static int
read_egen(void)
{
	fprintf(stderr, "!EGEN has not implemented yet\n");
	HECMW_abort(HECMW_comm_get_comm());
	return 0;
}
#endif

/*----------------------------------------------------------------------------*/

static int read_egrp_head(void) {
  int token;

  /* !EGROUP */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_H_EGROUP) {
    set_err_token(token, HECMW_IO_HEC_E0500, "!EGROUP required");
    return -1;
  }

  token = HECMW_heclex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_HEC_E0500, "',' required after !EGROUP");
    return -1;
  }

  return 0;
}

static int read_egrp_param_egrp(char *egrp) {
  int token;
  char *p;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E0500, "'=' required after EGRP");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_NAME) {
    set_err_token(token, HECMW_IO_HEC_E0500,
                  "EGRP must begin with a letter or '_'");
    return -1;
  }
  p = HECMW_heclex_get_text();
  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }
  strcpy(egrp, p);
  HECMW_toupper(egrp);
  if (HECMW_io_is_reserved_name(egrp)) {
    set_err(HECMW_IO_E0003, "");
    return -1;
  }
  if (strcmp(egrp, "ALL") == 0) {
    set_err(HECMW_IO_E0003, "Reserved name: %s", egrp);
    return -1;
  }
  return 0;
}

static int read_egrp_data(char *egrp) {
  int i, n, *elem, token;
  struct hecmw_io_id *head, *prev, *p, *q;

  n    = 0;
  prev = NULL;
  head = NULL;
  while (1) {
    struct hecmw_io_id *id;

    token = HECMW_heclex_next_token();
    if (n != 0 && token == HECMW_HECLEX_NL) break;

    id = HECMW_malloc(sizeof(*id));
    if (id == NULL) {
      HECMW_set_error(errno, "");
      return -1;
    }

    /* elemX */
    if (token != HECMW_HECLEX_INT) {
      set_err_token(token, HECMW_IO_HEC_E0500, "Element ID required");
      return -1;
    }
    id->id   = HECMW_heclex_get_number();
    id->next = NULL;
    if (head == NULL) {
      head = id;
    } else {
      prev->next = id;
    }
    prev = id;
    n++;

    /* ',' or NL */
    token = HECMW_heclex_next_token();
    if (token != ',' && token != HECMW_HECLEX_NL) {
      set_err_token(token, HECMW_IO_HEC_E0500,
                    "',' or NL required after element ID");
      return -1;
    }
    if (token == HECMW_HECLEX_NL) break;
  }
  HECMW_assert(head);
  HECMW_assert(n > 0);

  /* add elem to group */
  elem = HECMW_malloc(sizeof(*elem) * n);
  if (elem == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }
  i = 0;
  p = head;
  while (p) {
    elem[i++] = p->id;
    q         = p;
    p         = p->next;
    HECMW_free(q);
  }
  if (HECMW_io_add_egrp(egrp, n, elem) < 0) {
    return -1;
  }
  HECMW_free(elem);

  return 0;
}

static int read_egrp_generate(char *egrp) {
  int i, n, id, *elem, token;
  int elem1, elem2, elem3;

  while (1) {
    /* elem1 */
    token = HECMW_heclex_next_token();
    if (token != HECMW_HECLEX_INT) {
      set_err_token(token, HECMW_IO_HEC_E0500, "elem1 required");
      return -1;
    }
    elem1 = HECMW_heclex_get_number();
    if (elem1 <= 0) {
      set_err(HECMW_IO_HEC_E0502, "");
      return -1;
    }

    /* ',' */
    token = HECMW_heclex_next_token();
    if (token != ',') {
      set_err_token(token, HECMW_IO_HEC_E0500, "',' required after elem1");
      return -1;
    }

    /* elem2 */
    token = HECMW_heclex_next_token();
    if (token != HECMW_HECLEX_INT) {
      set_err_token(token, HECMW_IO_HEC_E0500, "elem2 required");
      return -1;
    }
    elem2 = HECMW_heclex_get_number();
    if (elem2 <= 0) {
      set_err(HECMW_IO_HEC_E0502, "");
      return -1;
    }

    /* ',' or NL */
    token = HECMW_heclex_next_token();
    if (token == ',') {
      /* elem3 */
      token = HECMW_heclex_next_token();
      if (token != HECMW_HECLEX_INT) {
        set_err_token(token, HECMW_IO_HEC_E0500, "Increment required");
        return -1;
      }
      elem3 = HECMW_heclex_get_number();
      if (elem3 <= 0) {
        set_err(HECMW_IO_HEC_E0502, "");
        return -1;
      }

      /* NL */
      token = HECMW_heclex_next_token();
      if (token != HECMW_HECLEX_NL) {
        set_err_token(token, HECMW_IO_HEC_E0500, "NL required after increment");
        return -1;
      }
    } else if (token == HECMW_HECLEX_NL) {
      elem3 = 1;
    } else {
      set_err_token(token, HECMW_IO_HEC_E0500,
                    "',' or NL required after elem2");
      return -1;
    }
    HECMW_assert(token == HECMW_HECLEX_NL);

    /* make element */
    if (elem1 > elem2) {
      set_err(HECMW_IO_HEC_E0503,
              "Cannot generate between %d and %d with an increment of %d",
              elem1, elem2, elem3);
      return -1;
    }
    if ((elem2 - elem1) % elem3) {
      set_err(HECMW_IO_HEC_E0503,
              "Cannot generate between %d and %d with an increment of %d",
              elem1, elem2, elem3);
      return -1;
    }

    n    = (elem2 - elem1) / elem3 + 1;
    elem = HECMW_malloc(sizeof(*elem) * n);
    if (elem == NULL) {
      HECMW_set_error(errno, "");
      return -1;
    }

    i = 0;
    for (id = elem1; id <= elem2; id += elem3) {
      elem[i++] = id;
    }
    HECMW_assert(i == n);
    if (HECMW_io_add_egrp(egrp, n, elem) < 0) return -1;
    HECMW_free(elem);

    /* check next state */
    token = HECMW_heclex_next_token();
    if (token != HECMW_HECLEX_INT) {
      HECMW_heclex_unput_token();
      break;
    }
    HECMW_heclex_unput_token();
  }
  return 0;
}

static int read_egroup(void) {
  int token, state;
  int flag_egrp                 = 0; /* flag for EGRP */
  int flag_generate             = 0; /* flag for GENERATE */
  int flag_input                = 0; /* flag for INPUT */
  char egrp[HECMW_NAME_LEN + 1] = "";
  enum {
    ST_FINISHED,
    ST_HEADER_LINE,
    ST_HEADER_LINE_PARAM,
    ST_DATA_INCLUDE,
    ST_DATA_LINE,
    ST_DATA_LINE_GENERATE
  };

  state = ST_HEADER_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_HEADER_LINE) {
      if (read_egrp_head()) return -1;
      state = ST_HEADER_LINE_PARAM;
    } else if (state == ST_HEADER_LINE_PARAM) {
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_K_EGRP) {
        /* must */
        if (read_egrp_param_egrp(egrp)) return -1;
        flag_egrp = 1;
      } else if (token == HECMW_HECLEX_K_GENERATE) {
        /* oprtional */
        flag_generate = 1;
      } else if (token == HECMW_HECLEX_K_INPUT) {
        /* oprtional */
        if (read_input(HECMW_IO_HEC_E0500)) return -1;
        flag_input = 1;
      } else {
        set_err_token(token, HECMW_IO_HEC_E0500, "Unknown parameter");
        return -1;
      }

      /* check next parameter */
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_NL) {
        if (flag_input) {
          state = ST_DATA_INCLUDE;
        } else if (flag_generate) {
          state = ST_DATA_LINE_GENERATE;
        } else {
          state = ST_DATA_LINE;
        }
        /* check */
        if (!flag_egrp) {
          set_err(HECMW_IO_HEC_E0501, "");
          return -1;
        }
      } else if (token == ',') {
        ; /* continue this state */
      } else {
        set_err_token(token, HECMW_IO_HEC_E0500, "Unknown parameter");
        return -1;
      }
    } else if (state == ST_DATA_INCLUDE) {
      HECMW_assert(flag_input);
      HECMW_assert(flag_egrp);

      if (HECMW_heclex_switch_to_include(include_filename)) return -1;
      state = flag_generate ? ST_DATA_LINE_GENERATE : ST_DATA_LINE;
    } else if (state == ST_DATA_LINE) {
      HECMW_assert(!flag_generate);
      HECMW_assert(flag_egrp);

      if (read_egrp_data(egrp)) return -1;

      /* check next state */
      token = HECMW_heclex_next_token();
      if (token != HECMW_HECLEX_INT) {
        state = ST_FINISHED;
      }
      HECMW_heclex_unput_token();
    } else if (state == ST_DATA_LINE_GENERATE) {
      HECMW_assert(flag_generate);
      HECMW_assert(flag_egrp);

      if (read_egrp_generate(egrp)) return -1;
      state = ST_FINISHED;
    } else {
      HECMW_assert(0);
    }
  }
  HECMW_log(HECMW_LOG_DEBUG, "read_egroup done");
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_elem_head(void) {
  int token;

  /* !ELEMENT */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_H_ELEMENT) {
    set_err_token(token, HECMW_IO_HEC_E0600, "!ELEMENT required");
    return -1;
  }

  token = HECMW_heclex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_HEC_E0600, "',' required after !ELEMENT");
    return -1;
  }

  return 0;
}

static int read_elem_param_type(int *type) {
  int token;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E0600, "'=' required after TYPE");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_INT) {
    set_err_token(token, HECMW_IO_HEC_E0600, "Invalid TYPE");
    return -1;
  }
  *type = HECMW_heclex_get_number();
  if (HECMW_get_max_node(*type) == -1) {
    set_err(HECMW_IO_HEC_E0601, "Invalid type: %d", *type);
    return -1;
  }
  return 0;
}

static int read_elem_param_egrp(char *egrp) {
  int token;
  char *p;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E0600, "'=' required after EGRP");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_NAME) {
    set_err_token(token, HECMW_IO_HEC_E0600,
                  "EGRP must begin with a letter or '_'");
    return -1;
  }
  p = HECMW_heclex_get_text();
  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }
  strcpy(egrp, p);
  HECMW_toupper(egrp);
  if (HECMW_io_is_reserved_name(egrp)) {
    set_err(HECMW_IO_E0003, "");
    return -1;
  }
  if (strcmp(egrp, "ALL") == 0) {
    set_err(HECMW_IO_E0003, "Reserved name: %s", egrp);
    return -1;
  }
  return 0;
}

static int read_elem_param_nmatitem(int *nmatitem) {
  int token;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E0600, "'=' required after MATITEM");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_INT) {
    set_err_token(token, HECMW_IO_HEC_E0602, "");
    return -1;
  }
  *nmatitem = HECMW_heclex_get_number();
  if (*nmatitem < 0) {
    set_err_token(token, HECMW_IO_HEC_E0602, "");
    return -1;
  }
  return 0;
}

static int read_elem_data_conn(int *id, int nnode, int *node) {
  int token, i;

  /* element ID */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_INT) {
    set_err_token(token, HECMW_IO_HEC_E0603, "");
    return -1;
  }
  *id = HECMW_heclex_get_number();
  if (*id <= 0) {
    set_err_token(token, HECMW_IO_HEC_E0603, "");
    return -1;
  }

  /* ',' */
  token = HECMW_heclex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_HEC_E0600, "',' required after element ID");
    return -1;
  }

  /* connectivity */
  i = 0;
  while (1) {
    token = HECMW_heclex_next_token();
    if (i != 0 && token == HECMW_HECLEX_NL) continue;
    if (token != HECMW_HECLEX_INT) {
      set_err(HECMW_IO_HEC_E0604, "");
      return -1;
    }
    node[i] = HECMW_heclex_get_number();
    if (node[i] <= 0) {
      set_err(HECMW_IO_HEC_E0604, "");
      return -1;
    }

    if (i == nnode - 1) break;

    /* ',' or NL */
    token = HECMW_heclex_next_token();
    if (token != ',' && token != HECMW_HECLEX_NL) {
      set_err_token(token, HECMW_IO_HEC_E0600,
                    "',' or NL required after connectivity");
      return -1;
    }

    i++;
  }
  return 0;
}

static int read_elem_data_mat(int nmatitem, double *matitem) {
  int token, i;

  /* default value */
  for (i = 0; i < nmatitem; i++) {
    matitem[i] = 0.0;
  }

  /* MATITEM */
  for (i = 0; i < nmatitem; i++) {
    token = HECMW_heclex_next_token();
    if (token != HECMW_HECLEX_DOUBLE && token != HECMW_HECLEX_INT) {
      set_err_token(token, HECMW_IO_HEC_E0600, "required MATITEM");
      return -1;
    }
    matitem[i] = HECMW_heclex_get_number();

    /* ',' or NL */
    token = HECMW_heclex_next_token();
    if (token != ',' && token != HECMW_HECLEX_NL) {
      set_err_token(token, HECMW_IO_HEC_E0600, "',' or NL required after MAT");
      return -1;
    }
    if (i == nmatitem - 1) {
      if (token != HECMW_HECLEX_NL) {
        set_err_token(token, HECMW_IO_HEC_E0600, "NL required after MAT");
        return -1;
      }
    } else {
      if (token != ',') {
        set_err_token(token, HECMW_IO_HEC_E0600, "',' required after MAT");
        return -1;
      }
    }
  }
  return 0;
}

static int read_element(void) {
  int token, state;
  int id;
  int nnode                     = 0;
  int *node                     = NULL;
  double *matitem               = NULL;
  int nmatitem                  = 0;
  int type                      = -1;
  int flag_type                 = 0; /* flag for TYPE */
  int flag_egrp                 = 0; /* flag for EGRP */
  int flag_matitem              = 0; /* flag for MATITEM */
  int flag_input                = 0; /* flag for INPUT */
  char egrp[HECMW_NAME_LEN + 1] = "";
  enum {
    st_finished,
    st_header_line,
    st_header_line_param,
    st_prepare,
    st_data_include,
    st_data_line_conn,    /* read element ID and connectivity */
    st_data_line_matitem, /* read MATITEM */
    st_data_line_regist,
    st_finalize
  };

  state = st_header_line;
  while (state != st_finished) {
    if (state == st_header_line) {
      if (read_elem_head()) return -1;
      state = st_header_line_param;
    } else if (state == st_header_line_param) {
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_K_TYPE) {
        /* must */
        if (read_elem_param_type(&type)) return -1;
        flag_type = 1;
      } else if (token == HECMW_HECLEX_K_EGRP) {
        /* optional */
        if (read_elem_param_egrp(egrp)) return -1;
        flag_egrp = 1;
      } else if (token == HECMW_HECLEX_K_MATITEM) {
        /* optional */
        if (read_elem_param_nmatitem(&nmatitem)) return -1;
        flag_matitem = 1;
      } else if (token == HECMW_HECLEX_K_INPUT) {
        /* optional */
        if (read_input(HECMW_IO_HEC_E0600)) return -1;
        flag_input = 1;
      } else {
        set_err_token(token, HECMW_IO_HEC_E0600, "Unknown parameter");
        return -1;
      }

      /* check next state */
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_NL) {
        /* check TYPE */
        if (!flag_type) {
          set_err(HECMW_IO_HEC_E0606, "");
          return -1;
        }
        state = st_prepare;
      } else if (token == ',') {
        ; /* continue this state */
      } else {
        set_err_token(token, HECMW_IO_HEC_E0600, "Unknown parameter");
        return -1;
      }
    } else if (state == st_prepare) {
      HECMW_assert(flag_type);
      HECMW_assert(type != -1);

      /* get # of connectivity */
      nnode = HECMW_get_max_node(type);
      HECMW_assert(nnode > 0);
      HECMW_assert(nnode <= HECMW_MAX_NODE_MAX);

      node = HECMW_malloc(sizeof(*node) * nnode);
      if (node == NULL) {
        HECMW_set_error(errno, "");
        return -1;
      }

      /* nmatitem */
      HECMW_assert(nmatitem >= 0);

      if (flag_matitem && nmatitem) {
        HECMW_assert(nmatitem > 0);
        matitem = HECMW_malloc(sizeof(*matitem) * nmatitem);
        if (matitem == NULL) {
          HECMW_set_error(errno, "");
          return -1;
        }
      } else {
        matitem = NULL;
      }

      /* set next state */
      state = flag_input ? st_data_include : st_data_line_conn;
    } else if (state == st_data_include) {
      HECMW_assert(flag_input);
      HECMW_assert(flag_type);

      if (HECMW_heclex_switch_to_include(include_filename)) return -1;
      state = st_data_line_conn;
    } else if (state == st_data_line_conn) {
      HECMW_assert(flag_type);

      if (read_elem_data_conn(&id, nnode, node)) return -1;
      if (HECMW_convert_connectivity(connectivity_type, type, node)) return -1;

      /* check next state */
      token = HECMW_heclex_next_token();
      if (flag_matitem) {
        if (token != ',' && token != HECMW_HECLEX_NL) {
          set_err_token(token, HECMW_IO_HEC_E0600,
                        "',' or NL required after connectivity");
          return -1;
        }
        if (token == ',') {
          token = HECMW_heclex_next_token();
          if (token != HECMW_HECLEX_NL) {
            HECMW_heclex_unput_token();
          }
        }
        state = st_data_line_matitem;
      } else {
        if (token != HECMW_HECLEX_NL) {
          set_err_token(token, HECMW_IO_HEC_E0600, "NL required");
          return -1;
        }
        state = st_data_line_regist;
      }
    } else if (state == st_data_line_matitem) {
      HECMW_assert(flag_matitem);
      HECMW_assert(nmatitem > 0);
      HECMW_assert(matitem);
      HECMW_assert(flag_type);

      if (read_elem_data_mat(nmatitem, matitem)) return -1;
      state = st_data_line_regist;
    } else if (state == st_data_line_regist) {
      HECMW_assert(node);
      HECMW_assert(flag_type);

      /* add element */
      if (HECMW_io_add_elem(id, type, node, nmatitem, matitem) == NULL) {
        return -1;
      }

      /* add element to eroup */
      if (HECMW_io_add_egrp("ALL", 1, &id) < 0) return -1;

      if (flag_egrp) {
        if (HECMW_io_add_egrp(egrp, 1, &id) < 0) {
          return -1;
        }
      }

      /* check next state */
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_INT) {
        state = st_data_line_conn;
      } else {
        state = st_finalize;
      }
      HECMW_heclex_unput_token();
    } else if (state == st_finalize) {
      HECMW_free(node);
      HECMW_free(matitem);

      /* set next state */
      state = st_finished;
    } else {
      HECMW_assert(0);
    }
  }
  HECMW_log(HECMW_LOG_DEBUG, "read_element done");
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_equation_head(int *last_token) {
  int token;

  /* !EQUATION */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_H_EQUATION) {
    set_err_token(token, HECMW_IO_HEC_E0700, "!EQUATION required");
    return -1;
  }

  token = HECMW_heclex_next_token();
  if (token != ',' && token != HECMW_HECLEX_NL) {
    set_err_token(token, HECMW_IO_HEC_E0700,
                  "',' or NL required after !EQUATION");
    return -1;
  }
  *last_token = token;

  return 0;
}

static int read_equation_data_line1(int *neq, double *cnst) {
  int token;
  char *p;

  /* NEQ */
  token = HECMW_heclex_next_token();

  if (token != HECMW_HECLEX_INT && token != HECMW_HECLEX_NAME) {
    set_err_token(token, HECMW_IO_HEC_E0700, "required NEQ");
    return -1;
  }

  if (token == HECMW_HECLEX_NAME) {
    p = HECMW_heclex_get_text();
    if (strcmp(p, "link") == 0 || strcmp(p, "LINK") == 0) {
      *neq  = 2;
      *cnst = 0.0;
    }
    HECMW_heclex_unput_token();
    return 0;
  }

  *neq = HECMW_heclex_get_number();
  if (*neq < 2) {
    set_err(HECMW_IO_HEC_E0701, "");
    return -1;
  }

  /* ',' */
  token = HECMW_heclex_next_token();
  if (token == ',') {
    /* const */
    token = HECMW_heclex_next_token();
    if (token != HECMW_HECLEX_DOUBLE && token != HECMW_HECLEX_INT) {
      set_err_token(token, HECMW_IO_HEC_E0700, "required CONST");
      return -1;
    }
    *cnst = HECMW_heclex_get_number();
  } else {
    *cnst = 0.0;
    HECMW_heclex_unput_token();
  }

  /* NL */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_NL) {
    set_err_token(token, HECMW_IO_HEC_E0700, "NL required after NEQ");
    return -1;
  }

  return 0;
}

static int read_equation_data_line2(int neq, double cnst) {
  int i, token;
  int is_node     = 0;
  int is_ngrp     = 0;
  int is_link     = 0;
  int is_beam     = 0;
  const int NITEM = 100;
  char *p;
  struct hecmw_io_mpcitem *mpcitem;
  bool isAllDof = false;

  mpcitem = HECMW_malloc(sizeof(*mpcitem) * neq);
  if (mpcitem == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  token = HECMW_heclex_next_token();
  if (token == HECMW_HECLEX_NAME) {
    p = HECMW_heclex_get_text();
    if (strcmp(p, "link") == 0 || strcmp(p, "LINK") == 0) {
      is_link = 1;
    }
  }
  HECMW_heclex_unput_token();

  if (is_link == 0) {
    isAllDof = false;
    for (i = 0; i < neq; i++) {
      token = HECMW_heclex_next_token();
      if (i != 0 && token == HECMW_HECLEX_NL) break;

      /* nod */
      if (token == HECMW_HECLEX_INT) {
        if (is_ngrp) {
          set_err(HECMW_IO_HEC_E0702, "");
          return -1;
        }
        mpcitem[i].node = HECMW_heclex_get_number();
        strcpy(mpcitem[i].ngrp, "");
        is_node = 1;
      } else if (token == HECMW_HECLEX_NAME) {
        char *p = HECMW_heclex_get_text();
        if (is_node) {
          set_err(HECMW_IO_HEC_E0702, "");
          return -1;
        }
        if (strlen(p) > HECMW_NAME_LEN) {
          set_err(HECMW_IO_E0001, "");
          return -1;
        }
        strcpy(mpcitem[i].ngrp, p);
        HECMW_toupper(mpcitem[i].ngrp);
        if (HECMW_io_is_reserved_name(mpcitem[i].ngrp)) {
          set_err(HECMW_IO_E0003, "");
          return -1;
        }
        mpcitem[i].node = -1;
        is_ngrp         = 1;
      } else {
        set_err_token(token, HECMW_IO_HEC_E0700, "Node ID or NGRP required");
        return -1;
      }

      /* ',' */
      token = HECMW_heclex_next_token();
      if (token != ',') {
        set_err_token(token, HECMW_IO_HEC_E0700, "',' required after node");
        return -1;
      }

      /* DOF */
      token = HECMW_heclex_next_token();
      if (token != HECMW_HECLEX_INT) {
        set_err(HECMW_IO_HEC_E0703, "");
        return -1;
      }
      mpcitem[i].dof = HECMW_heclex_get_number();
      if (mpcitem[i].dof == 0) {
        isAllDof       = true;
        mpcitem[i].dof = 1;
      }
      if (HECMW_io_check_mpc_dof(mpcitem[i].dof)) {
        set_err(HECMW_IO_HEC_E0703, "");
        return -1;
      }

      /* ',' */
      token = HECMW_heclex_next_token();
      if (token != ',') {
        set_err_token(token, HECMW_IO_HEC_E0700, "',' required after DOF");
        return -1;
      }

      /* A */
      token = HECMW_heclex_next_token();
      if (token != HECMW_HECLEX_DOUBLE && token != HECMW_HECLEX_INT) {
        set_err_token(token, HECMW_IO_HEC_E0700, "A(coefficient) required ");
        return -1;
      }
      mpcitem[i].a = HECMW_heclex_get_number();

      /* ',' or NL */
      token = HECMW_heclex_next_token();
      if (token != ',' && token != HECMW_HECLEX_NL) {
        set_err_token(token, HECMW_IO_HEC_E0700,
                      "',' or NL required after coefficient");
        return -1;
      }
      if (token == ',' && i == NITEM - 1) {
        token = HECMW_heclex_next_token();
        if (token != HECMW_HECLEX_NL) {
          set_err_token(token, HECMW_IO_HEC_E0700, "NL required");
          return -1;
        }
        continue;
      }
      if (token == HECMW_HECLEX_NL) continue;
    }

    /* add */
    if (isAllDof) {
      for (i = 0; i < neq; i++) {
        mpcitem[i].dof = 1;
      }
      if (HECMW_io_add_mpc(neq, mpcitem, cnst) == NULL) return -1;
      for (i = 0; i < neq; i++) {
        mpcitem[i].dof = 2;
      }
      if (HECMW_io_add_mpc(neq, mpcitem, cnst) == NULL) return -1;
      for (i = 0; i < neq; i++) {
        mpcitem[i].dof = 3;
      }
      if (HECMW_io_add_mpc(neq, mpcitem, cnst) == NULL) return -1;

    } else {
      if (HECMW_io_add_mpc(neq, mpcitem, cnst) == NULL) return -1;
    }
    HECMW_free(mpcitem);

    /* link */
  } else if (is_link == 1) {
    token = HECMW_heclex_next_token();

    /* ',' */
    token = HECMW_heclex_next_token();
    if (token != ',') {
      set_err_token(token, HECMW_IO_HEC_E0700, "',' required after DOF");
      return -1;
    }

    token = HECMW_heclex_next_token();
    if (token != HECMW_HECLEX_INT) {
      return -1;
    }
    mpcitem[0].node = HECMW_heclex_get_number();
    strcpy(mpcitem[0].ngrp, "");
    mpcitem[0].a = 1.0;

    /* ',' */
    token = HECMW_heclex_next_token();
    if (token != ',') {
      set_err_token(token, HECMW_IO_HEC_E0700, "',' required after DOF");
      return -1;
    }

    token = HECMW_heclex_next_token();
    if (token != HECMW_HECLEX_INT) {
      return -1;
    }
    mpcitem[1].node = HECMW_heclex_get_number();
    strcpy(mpcitem[1].ngrp, "");
    mpcitem[1].a = -1.0;

    /* add 1 */
    mpcitem[0].dof = 1;
    mpcitem[1].dof = 1;
    if (HECMW_io_add_mpc(neq, mpcitem, cnst) == NULL) return -1;
    /* add 2 */
    mpcitem[0].dof = 2;
    mpcitem[1].dof = 2;
    if (HECMW_io_add_mpc(neq, mpcitem, cnst) == NULL) return -1;
    /* add 3 */
    mpcitem[0].dof = 3;
    mpcitem[1].dof = 3;
    if (HECMW_io_add_mpc(neq, mpcitem, cnst) == NULL) return -1;
    HECMW_free(mpcitem);

    token = HECMW_heclex_next_token();
    if (token != HECMW_HECLEX_NL) {
      return -1;
    }
  }

  return 0;
}

static int read_equation(void) {
  int token, state;
  int neq        = -1;
  double cnst    = 0.0;
  int flag_input = 0; /* flag for INPUT */
  char *p;
  enum {
    ST_FINISHED,
    ST_HEADER_LINE,
    ST_HEADER_LINE_PARAM,
    ST_DATA_INCLUDE,
    ST_DATA_LINE1,
    ST_DATA_LINE2
  };

  state = ST_HEADER_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_HEADER_LINE) {
      if (read_equation_head(&token)) return -1;
      if (token == ',') {
        state = ST_HEADER_LINE_PARAM;
      } else if (token == HECMW_HECLEX_NL) {
        state = ST_DATA_LINE1;
      } else {
        HECMW_assert(0);
      }
    } else if (state == ST_HEADER_LINE_PARAM) {
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_K_INPUT) {
        /* optional */
        if (read_input(HECMW_IO_HEC_E0700)) return -1;
        flag_input = 1;
      } else {
        set_err_token(token, HECMW_IO_HEC_E0700, "Unknown parameter");
        return -1;
      }

      /* check next state */
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_NL) {
        state = flag_input ? ST_DATA_INCLUDE : ST_DATA_LINE1;
      } else {
        set_err_token(token, HECMW_IO_HEC_E0700, "NL required");
        return -1;
      }
    } else if (state == ST_DATA_INCLUDE) {
      HECMW_assert(flag_input);

      if (HECMW_heclex_switch_to_include(include_filename)) return -1;
      state = ST_DATA_LINE1;
    } else if (state == ST_DATA_LINE1) {
      if (read_equation_data_line1(&neq, &cnst)) return -1;
      state = ST_DATA_LINE2;
    } else if (state == ST_DATA_LINE2) {
      HECMW_assert(neq != -1);
      if (read_equation_data_line2(neq, cnst)) return -1;

      /* check next state */
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_INT) {
        state = ST_DATA_LINE1;
      } else if (token == HECMW_HECLEX_NAME) {
        p = HECMW_heclex_get_text();
        if (strcmp(p, "link") == 0 || strcmp(p, "LINK") == 0) {
          state = ST_DATA_LINE1;
        } else {
          state = ST_FINISHED;
        }
      } else {
        state = ST_FINISHED;
      }
      HECMW_heclex_unput_token();
    } else {
      HECMW_assert(0);
    }
  }
  HECMW_log(HECMW_LOG_DEBUG, "read_equation done");
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_header(void) {
  int token, len;
  char *p;
  struct hecmw_io_header *header;

  header = HECMW_malloc(sizeof(struct hecmw_io_header));
  if (header == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  /* !HEADER */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_H_HEADER) {
    set_err_token(token, HECMW_IO_HEC_E0800, "!HEADER required");
    return -1;
  }

  /* get header data */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_HEADER) {
    set_err_token(token, HECMW_IO_HEC_E0800, "TITLE required after !HEADER");
    return -1;
  }
  p = HECMW_heclex_get_text();
  while (*p && *p == ' ') p++;
  if (p == NULL) p                = "";
  len                             = strlen(p);
  if (len > HECMW_HEADER_LEN) len = HECMW_HEADER_LEN;
  strncpy(header->header, p, len);
  header->header[len] = '\0';

  /* Note:
   * NL is ignored by LEX until the end of the header data.
   */

  /* Ignore the rest of the header data */
  while (HECMW_heclex_next_token() == HECMW_HECLEX_HEADER)
    ;
  HECMW_heclex_unput_token();

  /* set */
  HECMW_io_set_header(header);

  HECMW_log(HECMW_LOG_DEBUG, "read_header done");
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_include(void) {
  int token;

  /* !INCLUDE */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_H_INCLUDE) {
    set_err_token(token, HECMW_IO_HEC_E0900, "!INCLUDE required");
    return -1;
  }

  /* ',' */
  token = HECMW_heclex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_HEC_E0900, "',' required after !INCLUDE");
    return -1;
  }

  /* INPUT */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_K_INPUT) {
    set_err_token(token, HECMW_IO_HEC_E0901, "");
    return -1;
  }

  /* =filename */
  if (read_input(HECMW_IO_HEC_E0900)) return -1;

  /* NL */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_NL) {
    set_err_token(token, HECMW_IO_HEC_E0900, "NL required after INPUT value");
    return -1;
  }

  /* include */
  if (HECMW_heclex_switch_to_include(include_filename)) {
    return -1;
  }

  HECMW_log(HECMW_LOG_DEBUG, "read_include done");
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_initial_head(void) {
  int token;

  /* !INITIAL CONDITION */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_H_INITIAL) {
    set_err_token(token, HECMW_IO_HEC_E1000, "!INITIAL CONDITION required");
    return -1;
  }

  token = HECMW_heclex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_HEC_E1001, "");
    return -1;
  }

  return 0;
}

static int read_initial_param_type(int *type) {
  int token;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E1000, "'=' required after TYPE");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_K_TEMPERATURE) {
    set_err_token(token, HECMW_IO_HEC_E1000, "TEMPERATURE required");
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
  token = HECMW_heclex_next_token();
  if (token == HECMW_HECLEX_INT) {
    node = HECMW_heclex_get_number();
    if (node <= 0) {
      set_err(HECMW_IO_HEC_E1002, "");
      return -1;
    }
  } else if (token == HECMW_HECLEX_NAME) {
    ngrp = HECMW_heclex_get_text();
    if (strlen(ngrp) > HECMW_NAME_LEN) {
      set_err(HECMW_IO_E0001, "");
      return -1;
    }
    HECMW_toupper(ngrp);
    ngrp = HECMW_strdup(ngrp);
    if (ngrp == NULL) {
      HECMW_set_error(errno, "");
      return -1;
    }
  } else {
    set_err_token(token, HECMW_IO_HEC_E1000, "Node ID or NGROUP name required");
    return -1;
  }

  /* ',' */
  token = HECMW_heclex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_HEC_E1000, "',' required after node");
    return -1;
  }

  /* VAL */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_DOUBLE && token != HECMW_HECLEX_INT) {
    set_err_token(token, HECMW_IO_HEC_E1000, "VAL required");
    return -1;
  }
  val = HECMW_heclex_get_number();

  /* NL */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_NL) {
    set_err_token(token, HECMW_IO_HEC_E1000, "NL required after VAL");
    return -1;
  }

  /* add */
  HECMW_assert(type != -1);
  if (HECMW_io_add_initial(type, node, ngrp, val) == NULL) {
    return -1;
  };
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
    ST_HEADER_LINE,
    ST_HEADER_LINE_PARAM,
    ST_DATA_INCLUDE,
    ST_DATA_LINE
  };

  state = ST_HEADER_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_HEADER_LINE) {
      if (read_initial_head()) return -1;
      state = ST_HEADER_LINE_PARAM;
    } else if (state == ST_HEADER_LINE_PARAM) {
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_K_TYPE) {
        /* must */
        if (read_initial_param_type(&type)) return -1;
        flag_type = 1;
      } else if (token == HECMW_HECLEX_K_INPUT) {
        /* oprtional */
        if (read_input(HECMW_IO_HEC_E1000)) return -1;
        flag_input = 1;
      } else {
        set_err_token(token, HECMW_IO_HEC_E1000, "Unknown parameter");
        return -1;
      }

      /* check next parameter */
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_NL) {
        if (flag_input) {
          state = ST_DATA_INCLUDE;
        } else {
          state = ST_DATA_LINE;
        }
        /* check */
        if (!flag_type) {
          set_err(HECMW_IO_HEC_E1001, "");
          return -1;
        }
      } else if (token == ',') {
        ; /* continue this state */
      } else {
        set_err_token(token, HECMW_IO_HEC_E1000, "Unknown parameter");
        return -1;
      }
    } else if (state == ST_DATA_INCLUDE) {
      HECMW_assert(flag_input);
      HECMW_assert(flag_type);
      if (HECMW_heclex_switch_to_include(include_filename)) return -1;
      state = ST_DATA_LINE;
    } else if (state == ST_DATA_LINE) {
      HECMW_assert(flag_type);
      if (read_initial_data(type)) return -1;

      /* check next state */
      token = HECMW_heclex_next_token();
      if (token != HECMW_HECLEX_INT && token != HECMW_HECLEX_NAME) {
        state = ST_FINISHED;
      }
      HECMW_heclex_unput_token();
    } else {
      HECMW_assert(0);
    }
  }
  HECMW_log(HECMW_LOG_DEBUG, "read_initial done");
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_matitem_head(int *item, int *last_token) {
  int token;

  /* !ITEM */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_H_ITEM) {
    set_err_token(token, HECMW_IO_HEC_E1100, "!ITEM required");
    return -1;
  }

  /* '=' */
  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E1100, "'=' required after !ITEM");
    return -1;
  }

  /* ITEM */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_INT) {
    set_err_token(token, HECMW_IO_HEC_E1100, "required !ITEM value");
    return -1;
  }
  *item = HECMW_heclex_get_number();
  if (*item <= 0) {
    set_err(HECMW_IO_HEC_E1104, "");
    return -1;
  }

  /* ',' or NL */
  token = HECMW_heclex_next_token();
  if (token != ',' && token != HECMW_HECLEX_NL) {
    set_err_token(token, HECMW_IO_HEC_E1100, "',' or NL after !ITEM value");
    return -1;
  }
  *last_token = token;

  return 0;
}

static int read_matitem_param_subitem(int *subitem) {
  int token;

  /* optional */
  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E1100, "'=' required after SUBITEM");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_INT) {
    set_err_token(token, HECMW_IO_HEC_E1100, "SUBITEM value required");
    return -1;
  }
  *subitem = HECMW_heclex_get_number();
  if (*subitem <= 0) {
    set_err(HECMW_IO_HEC_E1106, "");
    return -1;
  }
  return 0;
}

static int read_matitem_data(int subitem, int depend_temp,
                             struct hecmw_io_matitem *matitem) {
  int i, token;
  double *val;
  double temp = 0.0;
  struct hecmw_io_matsubitem *p, *q, *msitem;

  msitem = HECMW_malloc(sizeof(*msitem));
  if (msitem == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  val = HECMW_malloc(sizeof(*val) * subitem);
  if (val == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  for (i = 0; i < subitem; i++) {
    val[i] = 0.0; /* default value */
  }

  for (i = 0; i < subitem; i++) {
    /* VAL */
    token = HECMW_heclex_next_token();
    if (token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
      val[i] = HECMW_heclex_get_number();
    } else if (token == ',') {
      HECMW_heclex_unput_token();
    } else if (token == HECMW_HECLEX_NL) {
      break;
    } else {
      set_err_token(token, HECMW_IO_HEC_E1100, "VAL or ',' or NL reuqired");
      return -1;
    }

    /* ',' or NL*/
    token = HECMW_heclex_next_token();
    if (token != ',' && token != HECMW_HECLEX_NL) {
      set_err_token(token, HECMW_IO_HEC_E1100, "',' or NL required after VAL");
      return -1;
    }
    if (token == HECMW_HECLEX_NL) break;
  }

  if (token != HECMW_HECLEX_NL) {
    /* TEMP or NL */
    token = HECMW_heclex_next_token();
    if (token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
      temp = HECMW_heclex_get_number();
    } else if (token == HECMW_HECLEX_NL) {
      ;
    } else {
      set_err_token(token, HECMW_IO_HEC_E1100, "Temperature or NL required");
      return -1;
    }
  }

  if (depend_temp) {
    /* previous TEMP */
    q = NULL;
    for (p = matitem->subitem; p; p = (q = p)->next)
      ;
    if (q && temp <= q->temp) {
      set_err(HECMW_IO_HEC_E1107, "");
      return -1;
    }
  }

  if (token != HECMW_HECLEX_NL) {
    /* NL */
    token = HECMW_heclex_next_token();
    if (token != HECMW_HECLEX_NL) {
      set_err_token(token, HECMW_IO_HEC_E1100, "NL required");
      return -1;
    }
  }

  /* set */
  msitem->val  = val;
  msitem->temp = temp;
  msitem->next = NULL;

  q = NULL;
  for (p = matitem->subitem; p; p = (q = p)->next)
    ;
  if (q == NULL) {
    matitem->subitem = msitem;
  } else {
    q->next = msitem;
  }

  return 0;
}

static int read_matitem(struct hecmw_io_matitem *matitem) {
  int token, state;
  int item         = -1;
  int subitem      = 1;
  int flag_item    = 0; /* flag for !ITEM */
  int flag_subitem = 0; /* flag for SUBITEM */
  int depend_temp  = 0;
  enum {
    ST_FINISHED,
    ST_HEADER_LINE,
    ST_HEADER_LINE_PARAM,
    ST_PREPARE,
    ST_DATA_LINE
  };

  HECMW_assert(matitem);

  state = ST_HEADER_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_HEADER_LINE) {
      if (read_matitem_head(&item, &token)) return -1;
      if (token == ',') {
        state = ST_HEADER_LINE_PARAM;
      } else if (token == HECMW_HECLEX_NL) {
        state = ST_PREPARE;
      }
      flag_item = 1;
    } else if (state == ST_HEADER_LINE_PARAM) {
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_K_SUBITEM) {
        if (read_matitem_param_subitem(&subitem)) return -1;
        flag_subitem = 1;
      } else {
        set_err_token(token, HECMW_IO_HEC_E1100, "Unknown parameter");
        return -1;
      }

      /* NL */
      token = HECMW_heclex_next_token();
      if (token != HECMW_HECLEX_NL) {
        set_err_token(token, HECMW_IO_HEC_E1100, "NL required after SUBITEM");
        return -1;
      }

      /* set next state */
      state = ST_PREPARE;
    } else if (state == ST_PREPARE) {
      HECMW_assert(flag_item);
      HECMW_assert(item > 0);
      HECMW_assert(subitem > 0);

      matitem->item    = item;
      matitem->nval    = subitem;
      matitem->subitem = NULL;

      depend_temp = 0;

      /* set next state */
      state = ST_DATA_LINE;
    } else if (state == ST_DATA_LINE) {
      HECMW_assert(flag_item);
      HECMW_assert(item > 0);
      HECMW_assert(subitem > 0);

      if (read_matitem_data(subitem, depend_temp, matitem)) return -1;
      depend_temp = 1;

      /* check next state */
      token = HECMW_heclex_next_token();
      if (token != HECMW_HECLEX_DOUBLE && token != HECMW_HECLEX_INT) {
        state = ST_FINISHED;
      }
      HECMW_heclex_unput_token();
    } else {
      HECMW_assert(0);
    }
  }
  return 0;
}

static int read_material_head(void) {
  int token;

  /* !MATERIAL */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_H_MATERIAL) {
    set_err_token(token, HECMW_IO_HEC_E1100, "!MATERIAL required");
    return -1;
  }

  token = HECMW_heclex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_HEC_E1101, "");
    return -1;
  }

  return 0;
}

static int read_material_param_name(char *name) {
  int token;
  char *p;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E1100, "'=' required after NAME");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_NAME) {
    set_err_token(token, HECMW_IO_HEC_E1100,
                  "NAME must begin with a letter or '_'");
    return -1;
  }
  p = HECMW_heclex_get_text();
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
    set_err(HECMW_IO_HEC_E1102, "%s already exists", name);
    return -1;
  }
  return 0;
}

static int read_material_param_item(int *item) {
  int token;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E1100, "'=' required after ITEM");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_INT) {
    set_err_token(token, HECMW_IO_HEC_E1100, "Invalid ITEM");
    return -1;
  }
  *item = HECMW_heclex_get_number();
  if (*item <= 0) {
    set_err(HECMW_IO_HEC_E1103, "");
    return -1;
  }
  return 0;
}

static int matitem_comp(const void *matitem1, const void *matitem2) {
  const struct hecmw_io_matitem *m1, *m2;

  m1 = matitem1;
  m2 = matitem2;

  if (m1->item == m2->item) return 0;
  if (m1->item < m2->item) {
    return -1;
  } else {
    return 1;
  }
}

static int read_material_data(int item, char *name) {
  int i;
  struct hecmw_io_material *mat;
  struct hecmw_io_matitem *matitem;

  mat = HECMW_malloc(sizeof(*mat));
  if (mat == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  matitem = HECMW_malloc(sizeof(*matitem) * item);
  if (matitem == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }

  for (i = 0; i < item; i++) {
    if (read_matitem(&matitem[i])) return -1;
  }

  /* sort */
  qsort(matitem, item, sizeof(*matitem), matitem_comp);

  /* check !ITEM  value */
  for (i = 0; i < item; i++) {
    if (matitem[i].item != i + 1) {
      HECMW_set_error(HECMW_IO_HEC_E1105, "In MATERIAL %s", name);
      return -1;
    }
  }

  /* set */
  strcpy(mat->name, name);
  mat->nitem = item;
  mat->item  = matitem;
  mat->next  = NULL;

  /* add */
  if (HECMW_io_add_mat(name, mat) == NULL) return -1;

  return 0;
}

static int read_material(void) {
  int token, state;
  int item                      = 1;
  int flag_name                 = 0; /* flag for NAME */
  int flag_item                 = 0; /* flag for ITEM */
  int flag_input                = 0; /* flag for INPUT */
  char name[HECMW_NAME_LEN + 1] = "";
  enum {
    ST_FINISHED,
    ST_HEADER_LINE,
    ST_HEADER_LINE_PARAM,
    ST_DATA_INCLUDE,
    ST_DATA_LINE
  };

  state = ST_HEADER_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_HEADER_LINE) {
      if (read_material_head()) return -1;
      state = ST_HEADER_LINE_PARAM;
    } else if (state == ST_HEADER_LINE_PARAM) {
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_K_NAME) {
        /* must */
        if (read_material_param_name(name)) return -1;
        flag_name = 1;
      } else if (token == HECMW_HECLEX_K_ITEM) {
        /* optioanal */
        if (read_material_param_item(&item)) return -1;
        flag_item = 1;
      } else if (token == HECMW_HECLEX_K_INPUT) {
        /* oprtional */
        if (read_input(HECMW_IO_HEC_E1100)) return -1;
        flag_input = 1;
      } else {
        set_err_token(token, HECMW_IO_HEC_E1100, "Unknown parameter");
        return -1;
      }

      /* check next parameter */
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_NL) {
        /* check */
        if (!flag_name) {
          set_err(HECMW_IO_HEC_E1101, "");
          return -1;
        }
        state = flag_input ? ST_DATA_INCLUDE : ST_DATA_LINE;
      } else if (token == ',') {
        ; /* continue this state */
      } else {
        set_err_token(token, HECMW_IO_HEC_E1100, "Unknown parameter");
        return -1;
      }
    } else if (state == ST_DATA_INCLUDE) {
      HECMW_assert(flag_input);
      HECMW_assert(flag_name);
      if (HECMW_heclex_switch_to_include(include_filename)) return -1;
      state = ST_DATA_LINE;
    } else if (state == ST_DATA_LINE) {
      HECMW_assert(flag_name);
      if (read_material_data(item, name)) return -1;
      state = ST_FINISHED;
    } else {
      HECMW_assert(0);
    }
  }
  HECMW_log(HECMW_LOG_DEBUG, "read_material done");
  return 0;
}

/*----------------------------------------------------------------------------*/
#if 0
static int
read_ncopy(void)
{
	fprintf(stderr, "!NCOPY has not implemented yet\n");
	HECMW_abort(HECMW_comm_get_comm());
	return 0;
}


static int
read_nfill(void)
{
	fprintf(stderr, "!NFILL has not implemented yet\n");
	HECMW_abort(HECMW_comm_get_comm());
	return 0;
}


static int
read_ngen(void)
{
	fprintf(stderr, "!NGEN has not implemented yet\n");
	HECMW_abort(HECMW_comm_get_comm());
	return 0;
}
#endif

/*----------------------------------------------------------------------------*/

static int read_ngrp_head(void) {
  int token;

  /* !NGROUP */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_H_NGROUP) {
    set_err_token(token, HECMW_IO_HEC_E1500, "!NGROUP required");
    return -1;
  }

  token = HECMW_heclex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_HEC_E1500, "',' required after !NGROUP");
    return -1;
  }

  return 0;
}

static int read_ngrp_param_ngrp(char *ngrp) {
  int token;
  char *p;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E1500, "'=' required after NGRP");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_NAME) {
    set_err_token(token, HECMW_IO_HEC_E1500,
                  "NGRP must begin with a letter or '_'");
    return -1;
  }
  p = HECMW_heclex_get_text();
  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }
  strcpy(ngrp, p);
  HECMW_toupper(ngrp);
  if (HECMW_io_is_reserved_name(ngrp)) {
    set_err(HECMW_IO_E0003, "");
    return -1;
  }
  if (strcmp(ngrp, "EQUATION_BLOCK") == 0) {
    HECMW_set_error(HECMW_IO_E0003, "Reserved name: %s", ngrp);
    return -1;
  }
  if (strcmp(ngrp, "ALL") == 0) {
    HECMW_set_error(HECMW_IO_E0003, "Reserved name: %s", ngrp);
    return -1;
  }
  return 0;
}

static int read_ngrp_data(char *ngrp) {
  int i, n, *node, token;
  struct hecmw_io_id *head, *prev, *p, *q;

  n    = 0;
  prev = NULL;
  head = NULL;
  while (1) {
    struct hecmw_io_id *id;

    token = HECMW_heclex_next_token();
    if (n != 0 && token == HECMW_HECLEX_NL) break;

    id = HECMW_malloc(sizeof(*id));
    if (id == NULL) {
      HECMW_set_error(errno, "");
      return -1;
    }

    /* nodX */
    if (token != HECMW_HECLEX_INT) {
      set_err_token(token, HECMW_IO_HEC_E1500, "Node ID required");
      return -1;
    }
    id->id   = HECMW_heclex_get_number();
    id->next = NULL;
    if (head == NULL) {
      head = id;
    } else {
      prev->next = id;
    }
    prev = id;
    n++;

    /* ',' or NL */
    token = HECMW_heclex_next_token();
    if (token != ',' && token != HECMW_HECLEX_NL) {
      set_err_token(token, HECMW_IO_HEC_E1500,
                    "',' or NL required after node ID");
      return -1;
    }
    if (token == HECMW_HECLEX_NL) break;
  }
  HECMW_assert(head);
  HECMW_assert(n > 0);

  /* add node to group */
  node = HECMW_malloc(sizeof(*node) * n);
  if (node == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }
  i = 0;
  p = head;
  while (p) {
    node[i++] = p->id;
    q         = p;
    p         = p->next;
    HECMW_free(q);
  }
  if (HECMW_io_add_ngrp(ngrp, n, node) < 0) return -1;
  HECMW_free(node);

  return 0;
}

static int read_ngrp_generate(char *ngrp) {
  int i, n, id, *node, token;
  int nod1, nod2, nod3;

  while (1) {
    /* nod1 */
    token = HECMW_heclex_next_token();
    if (token != HECMW_HECLEX_INT) {
      set_err_token(token, HECMW_IO_HEC_E1500, "nod1 required");
      return -1;
    }
    nod1 = HECMW_heclex_get_number();
    if (nod1 <= 0) {
      set_err(HECMW_IO_HEC_E1502, "");
      return -1;
    }

    /* ',' */
    token = HECMW_heclex_next_token();
    if (token != ',') {
      set_err_token(token, HECMW_IO_HEC_E1500, "',' required after nod1");
      return -1;
    }

    /* nod2 */
    token = HECMW_heclex_next_token();
    if (token != HECMW_HECLEX_INT) {
      set_err_token(token, HECMW_IO_HEC_E1500, "nod2 required");
      return -1;
    }
    nod2 = HECMW_heclex_get_number();
    if (nod2 <= 0) {
      set_err(HECMW_IO_HEC_E1502, "");
      return -1;
    }

    /* ',' or NL */
    token = HECMW_heclex_next_token();
    if (token == ',') {
      /* nod3 */
      token = HECMW_heclex_next_token();
      if (token != HECMW_HECLEX_INT) {
        set_err_token(token, HECMW_IO_HEC_E1500, "Increment required");
        return -1;
      }
      nod3 = HECMW_heclex_get_number();
      if (nod3 <= 0) {
        set_err(HECMW_IO_HEC_E1502, "");
        return -1;
      }

      /* NL */
      token = HECMW_heclex_next_token();
      if (token != HECMW_HECLEX_NL) {
        set_err_token(token, HECMW_IO_HEC_E1500, "NL required after increment");
        return -1;
      }
    } else if (token == HECMW_HECLEX_NL) {
      nod3 = 1;
    } else {
      set_err_token(token, HECMW_IO_HEC_E1500, "',' or NL required after nod2");
      return -1;
    }
    HECMW_assert(token == HECMW_HECLEX_NL);

    /* make node */
    if (nod1 > nod2) {
      set_err(HECMW_IO_HEC_E1503,
              "Cannot generate between %d and %d with an increment of %d", nod1,
              nod2, nod3);
      return -1;
    }
    if ((nod2 - nod1) % nod3) {
      set_err(HECMW_IO_HEC_E1503,
              "Cannot generate between %d and %d with an increment of %d", nod1,
              nod2, nod3);
      return -1;
    }

    n    = (nod2 - nod1) / nod3 + 1;
    node = HECMW_malloc(sizeof(int) * n);
    if (node == NULL) {
      HECMW_set_error(errno, "");
      return -1;
    }

    i = 0;
    for (id = nod1; id <= nod2; id += nod3) {
      node[i++] = id;
    }
    HECMW_assert(i == n);
    if (HECMW_io_add_ngrp(ngrp, n, node) < 0) return -1;
    HECMW_free(node);

    /* check next state */
    token = HECMW_heclex_next_token();
    if (token != HECMW_HECLEX_INT) {
      HECMW_heclex_unput_token();
      break;
    }
    HECMW_heclex_unput_token();
  }
  return 0;
}

static int read_ngroup(void) {
  int token, state;
  int flag_ngrp                 = 0; /* flag for NGRP */
  int flag_generate             = 0; /* flag for GENERATE */
  int flag_input                = 0; /* flag for INPUT */
  char ngrp[HECMW_NAME_LEN + 1] = "";
  enum {
    st_finished,
    st_header_line,
    st_header_line_param,
    st_data_include,
    st_data_line,
    st_data_line_generate
  };

  state = st_header_line;
  while (state != st_finished) {
    if (state == st_header_line) {
      if (read_ngrp_head()) return -1;
      state = st_header_line_param;
    } else if (state == st_header_line_param) {
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_K_NGRP) {
        /* must */
        if (read_ngrp_param_ngrp(ngrp)) return -1;
        flag_ngrp = 1;
      } else if (token == HECMW_HECLEX_K_GENERATE) {
        /* oprtional */
        flag_generate = 1;
      } else if (token == HECMW_HECLEX_K_INPUT) {
        /* oprtional */
        if (read_input(HECMW_IO_HEC_E1500)) return -1;
        flag_input = 1;
      } else {
        set_err_token(token, HECMW_IO_HEC_E1500, "Unknown parameter");
        return -1;
      }

      /* check next parameter */
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_NL) {
        if (flag_input) {
          state = st_data_include;
        } else if (flag_generate) {
          state = st_data_line_generate;
        } else {
          state = st_data_line;
        }
        /* check */
        if (!flag_ngrp) {
          set_err(HECMW_IO_HEC_E1501, "");
          return -1;
        }
      } else if (token == ',') {
        ; /* continue this state */
      } else {
        set_err_token(token, HECMW_IO_HEC_E1500, "Unknown parameter");
        return -1;
      }
    } else if (state == st_data_include) {
      HECMW_assert(flag_input);
      HECMW_assert(flag_ngrp);
      if (HECMW_heclex_switch_to_include(include_filename)) return -1;
      state = flag_generate ? st_data_line_generate : st_data_line;
    } else if (state == st_data_line) {
      HECMW_assert(flag_ngrp);
      if (read_ngrp_data(ngrp)) return -1;

      /* check next state */
      token = HECMW_heclex_next_token();
      if (token != HECMW_HECLEX_INT) {
        state = st_finished;
      }
      HECMW_heclex_unput_token();
    } else if (state == st_data_line_generate) {
      HECMW_assert(flag_generate);
      HECMW_assert(flag_ngrp);
      if (read_ngrp_generate(ngrp)) return -1;
      state = st_finished;
    } else {
      HECMW_assert(0);
    }
  }
  HECMW_log(HECMW_LOG_DEBUG, "read_ngroup done");
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_node_head(int *last_token) {
  int token;

  /* !NODE */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_H_NODE) {
    set_err_token(token, HECMW_IO_HEC_E1600, "!NODE required");
    return -1;
  }

  token = HECMW_heclex_next_token();
  if (token != ',' && token != HECMW_HECLEX_NL) {
    set_err_token(token, HECMW_IO_HEC_E1600, "',' or NL required after !NODE");
    return -1;
  }
  *last_token = token;

  return 0;
}

static int read_node_param_system(int *system) {
  int token;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E1600, "'=' required after SYSTEM");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != 'C' && token != 'R') {
    set_err_token(token, HECMW_IO_HEC_E1600, "Invalid SYSTEM");
    return -1;
  }
  *system = token;
  return 0;
}

static int read_node_param_ngrp(char *ngrp) {
  int token;
  char *p;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E1600, "'=' required after NGRP");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_NAME) {
    set_err_token(token, HECMW_IO_HEC_E1600,
                  "NGRP must begin with a letter or '_'");
    return -1;
  }
  p = HECMW_heclex_get_text();
  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }
  strcpy(ngrp, p);
  HECMW_toupper(ngrp);
  if (HECMW_io_is_reserved_name(ngrp)) {
    set_err(HECMW_IO_E0003, "");
    return -1;
  }
  if (strcmp(ngrp, "ALL") == 0) {
    HECMW_set_error(HECMW_IO_E0003, "Reserved name: %s", ngrp);
    return -1;
  }
  return 0;
}

static int read_node_data(int *id_arg, double *x_arg, double *y_arg,
                          double *z_arg) {
  int id, token;
  double x, y, z;

  /* node ID */
  token = HECMW_heclex_next_token();
  if (token == HECMW_HECLEX_INT) {
    id = HECMW_heclex_get_number();
    if (id <= 0) {
      set_err(HECMW_IO_HEC_E1601, "");
      return -1;
    }
  } else {
    set_err(HECMW_IO_HEC_E1601, "");
    return -1;
  }

  /* ',' */
  token = HECMW_heclex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_HEC_E1600, "',' required after nood ID");
    return -1;
  }

  x = y = z = 0.0;
  while (1) {
    /* X */
    token = HECMW_heclex_next_token();
    if (token == HECMW_HECLEX_NL) break;
    if (token == ',') {
      HECMW_heclex_unput_token();
    } else if (token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
      x = HECMW_heclex_get_number();
    } else {
      set_err_token(token, HECMW_IO_HEC_E1600, "X required");
      return -1;
    }

    /* ',' */
    token = HECMW_heclex_next_token();
    if (token == HECMW_HECLEX_NL) break;
    if (token != ',') {
      set_err_token(token, HECMW_IO_HEC_E1600, "',' required after X");
      return -1;
    }

    /* Y */
    token = HECMW_heclex_next_token();
    if (token == HECMW_HECLEX_NL) break;
    if (token == ',') {
      HECMW_heclex_unput_token();
    } else if (token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
      y = HECMW_heclex_get_number();
    } else {
      set_err_token(token, HECMW_IO_HEC_E1600, "Y required");
      return -1;
    }

    /* ',' */
    token = HECMW_heclex_next_token();
    if (token == HECMW_HECLEX_NL) break;
    if (token != ',') {
      set_err_token(token, HECMW_IO_HEC_E1600, "',' required after Y");
      return -1;
    }

    /* Z */
    token = HECMW_heclex_next_token();
    if (token == HECMW_HECLEX_NL) break;
    if (token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
      z = HECMW_heclex_get_number();
    } else {
      set_err_token(token, HECMW_IO_HEC_E1600, "Z required");
      return -1;
    }

    /* ',' or NL */
    token = HECMW_heclex_next_token();
    if (token == HECMW_HECLEX_NL) break;
    if (token == ',') {
      token = HECMW_heclex_next_token();
      if (token != HECMW_HECLEX_NL) {
        set_err_token(token, HECMW_IO_HEC_E1600, "NL required after Z");
        return -1;
      }
    }

    break;
  }

  *id_arg = id;
  *x_arg  = x;
  *y_arg  = y;
  *z_arg  = z;

  return 0;
}

static int read_node_convert_coord(int system, double *x, double *y,
                                   double *z) {
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

  /* reflect !SYSTEM */
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
  int flag_ngrp                 = 0; /* flag for NGRP */
  int flag_input                = 0; /* flag for INPUT */
  char ngrp[HECMW_NAME_LEN + 1] = "";
  enum {
    ST_FINISHED,
    ST_HEADER_LINE,
    ST_HEADER_LINE_PARAM,
    ST_DATA_INCLUDE,
    ST_DATA_LINE
  };

  state = ST_HEADER_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_HEADER_LINE) {
      if (read_node_head(&token)) return -1;
      if (token == HECMW_HECLEX_NL) {
        state = ST_DATA_LINE;
      } else if (token == ',') {
        state = ST_HEADER_LINE_PARAM;
      } else {
        HECMW_assert(0);
      }
    } else if (state == ST_HEADER_LINE_PARAM) {
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_K_SYSTEM) {
        /* optional */
        if (read_node_param_system(&system)) return -1;
        flag_system = 1;
      } else if (token == HECMW_HECLEX_K_NGRP) {
        /* optional */
        if (read_node_param_ngrp(ngrp)) return -1;
        flag_ngrp = 1;
      } else if (token == HECMW_HECLEX_K_INPUT) {
        /* optional */
        if (read_input(HECMW_IO_HEC_E1600)) return -1;
        flag_input = 1;
      } else {
        set_err_token(token, HECMW_IO_HEC_E1600, "Unknown parameter");
        return -1;
      }

      /* check next parameter */
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_NL) {
        state = flag_input ? ST_DATA_INCLUDE : ST_DATA_LINE;
      } else if (token == ',') {
        ; /* continue this state */
      } else {
        set_err_token(token, HECMW_IO_HEC_E1600, "Unknown parameter");
        return -1;
      }
    } else if (state == ST_DATA_INCLUDE) {
      HECMW_assert(flag_input);
      if (HECMW_heclex_switch_to_include(include_filename)) return -1;
      state = ST_DATA_LINE;
    } else if (state == ST_DATA_LINE) {
      int id;
      double x, y, z;

      if (read_node_data(&id, &x, &y, &z)) return -1;
      if (read_node_convert_coord(system, &x, &y, &z)) return -1;

      /* add node */
      if (HECMW_io_add_node(id, x, y, z) == NULL) return -1;

      /* add node to group */
      if (HECMW_io_add_ngrp("ALL", 1, &id) < 0) return -1;

      if (flag_ngrp) {
        if (HECMW_io_add_ngrp(ngrp, 1, &id) < 0) return -1;
      }

      /* check next state */
      token = HECMW_heclex_next_token();
      if (token != HECMW_HECLEX_INT) {
        state = ST_FINISHED;
      } else {
        state = ST_DATA_LINE;
      }
      HECMW_heclex_unput_token();
    } else {
      HECMW_assert(0);
    }
  }
  HECMW_log(HECMW_LOG_DEBUG, "read_node done");
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_section_head(void) {
  int token;

  /* !SECTION */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_H_SECTION) {
    set_err_token(token, HECMW_IO_HEC_E1700, "!SECTION required");
    return -1;
  }

  token = HECMW_heclex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_HEC_E1700, "',' required after !SECTION");
    return -1;
  }

  return 0;
}

static int read_section_param_type(int *type) {
  int token;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E1700, "'=' required after TYPE");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token == HECMW_HECLEX_K_SOLID) {
    *type = HECMW_SECT_TYPE_SOLID;
  } else if (token == HECMW_HECLEX_K_SHELL) {
    *type = HECMW_SECT_TYPE_SHELL;
  } else if (token == HECMW_HECLEX_K_BEAM) {
    *type = HECMW_SECT_TYPE_BEAM;
  } else if (token == HECMW_HECLEX_K_INTERFACE) {
    *type = HECMW_SECT_TYPE_INTERFACE;
  } else {
    set_err_token(token, HECMW_IO_HEC_E1700, "Invalid  TYPE");
    return -1;
  }
  return 0;
}

static int read_section_param_egrp(char *egrp) {
  int token;
  char *p;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E1700, "'=' reuqired after EGRP");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_NAME) {
    set_err_token(token, HECMW_IO_HEC_E1700,
                  "EGRP must begin with a letter or '_'");
    return -1;
  }
  p = HECMW_heclex_get_text();
  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }
  strcpy(egrp, p);
  HECMW_toupper(egrp);
  if (HECMW_io_is_reserved_name(egrp)) {
    set_err(HECMW_IO_E0003, "");
    return -1;
  }
  return 0;
}

static int read_section_param_material(char *material) {
  int token;
  char *p;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E1700, "'=' reuqired after MATERIAL");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_NAME) {
    set_err_token(token, HECMW_IO_HEC_E1700,
                  "MATERIAL must begin with a letter or '_'");
    return -1;
  }
  p = HECMW_heclex_get_text();
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

#if 0
static int
read_section_param_composite(int *composite)
{
	int token;

	*composite = 1;	/* default value */
	token = HECMW_heclex_next_token();
	if(token == '=') {
		token = HECMW_heclex_next_token();
		if(token != HECMW_HECLEX_INT) {
			set_err_token(token, HECMW_IO_HEC_E1700, "COMPOSITE value reuqired");
			return -1;
		}
		*composite = HECMW_heclex_get_number();
	} else {
		HECMW_heclex_unput_token();
	}
	if(*composite <= 0) {
		set_err(HECMW_IO_HEC_E1704, "");
		return -1;
	}
	return 0;
}
#endif

static int read_section_param_secopt(int *secopt_arg) {
  int token, secopt;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E1700, "'=' required after SECOPT");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_INT) {
    set_err_token(token, HECMW_IO_HEC_E1700, "SECOPT value reuqired");
    return -1;
  }
  secopt = HECMW_heclex_get_number();
  if (secopt != HECMW_SECT_OPT_PSTRESS && secopt != HECMW_SECT_OPT_PSTRAIN &&
      secopt != HECMW_SECT_OPT_ASYMMETRY &&
      secopt != HECMW_SECT_OPT_PSTRESS_RI &&
      secopt != HECMW_SECT_OPT_PSTRAIN_RI &&
      secopt != HECMW_SECT_OPT_ASYMMETRY_RI) {
    set_err_token(token, HECMW_IO_HEC_E1700, "Invalid SECOPT");
    return -1;
  }

  *secopt_arg = secopt;

  return 0;
}

static int read_section_solid(union hecmw_io_section_item *sect_item) {
  int token;
  double thickness;

  /* THICKNESS */
  token = HECMW_heclex_next_token();
  if (token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
    thickness = HECMW_heclex_get_number();

    /* NL */
    token = HECMW_heclex_next_token();
    if (token != HECMW_HECLEX_NL) {
      set_err_token(token, HECMW_IO_HEC_E1700, "NL required after THICKNESS");
      return -1;
    }
  } else {
    thickness = 1.0;
    HECMW_heclex_unput_token();
  }

  if (thickness <= 0.0) {
    set_err(HECMW_IO_HEC_E1705, "");
    return -1;
  }

  /* set */
  sect_item->solid.thickness = thickness;

  return 0;
}

static int read_section_shell(union hecmw_io_section_item *sect_item) {
  double thickness;
  int token, integpoints;

  /* THICKNESS */
  token = HECMW_heclex_next_token();
  if (token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
    thickness = HECMW_heclex_get_number();
  } else {
    set_err_token(token, HECMW_IO_HEC_E1700, "THICKNESS reuiqred");
    return -1;
  }
  if (thickness <= 0.0) {
    set_err(HECMW_IO_HEC_E1705, "");
    return -1;
  }

  /* ',' */
  token = HECMW_heclex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_HEC_E1700, "',' required after THICKNESS");
    return -1;
  }

  /* INTEGPOINTS */
  token = HECMW_heclex_next_token();
  if (token == HECMW_HECLEX_INT) {
    integpoints = HECMW_heclex_get_number();
  } else {
    set_err_token(token, HECMW_IO_HEC_E1700, "INTEGPOINTS required");
    return -1;
  }
  if (integpoints <= 0) {
    set_err(HECMW_IO_HEC_E1706, "");
    return -1;
  }

  /* NL */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_NL) {
    set_err_token(token, HECMW_IO_HEC_E1700, "NL required after INTEGPOINTS");
    return -1;
  }

  /* set */
  sect_item->shell.thickness   = thickness;
  sect_item->shell.integpoints = integpoints;

  return 0;
}

static int read_section_beam(union hecmw_io_section_item *sect_item) {
  double nx, ny, nz, area, Iyy, Izz, Jx;
  int token;

  /* vx */
  token = HECMW_heclex_next_token();
  if (token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
    nx = HECMW_heclex_get_number();
  } else {
    set_err_token(token, HECMW_IO_HEC_E1700, "vx reuiqred");
    return -1;
  }

  /* ',' */
  token = HECMW_heclex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_HEC_E1700, "',' required after vx");
    return -1;
  }

  /* vy */
  token = HECMW_heclex_next_token();
  if (token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
    ny = HECMW_heclex_get_number();
  } else {
    set_err_token(token, HECMW_IO_HEC_E1700, "vy reuiqred");
    return -1;
  }

  /* ',' */
  token = HECMW_heclex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_HEC_E1700, "',' required after vy");
    return -1;
  }

  /* vz */
  token = HECMW_heclex_next_token();
  if (token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
    nz = HECMW_heclex_get_number();
  } else {
    set_err_token(token, HECMW_IO_HEC_E1700, "vz reuiqred");
    return -1;
  }

  /* ',' */
  token = HECMW_heclex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_HEC_E1700, "',' required after vz");
    return -1;
  }

  /* area */
  token = HECMW_heclex_next_token();
  if (token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
    area = HECMW_heclex_get_number();
  } else {
    set_err_token(token, HECMW_IO_HEC_E1700, "area required");
    return -1;
  }
  if (area <= 0) {
    set_err(HECMW_IO_HEC_E1707, "");
    return -1;
  }

  /* ',' */
  token = HECMW_heclex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_HEC_E1700, "',' required after vz");
    return -1;
  }

  /* Iyy */
  token = HECMW_heclex_next_token();
  if (token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
    Iyy = HECMW_heclex_get_number();
  } else {
    set_err_token(token, HECMW_IO_HEC_E1700, "Iyy reuiqred");
    return -1;
  }
  if (Iyy <= 0) {
    set_err(HECMW_IO_HEC_E1708, "");
    return -1;
  }

  /* ',' */
  token = HECMW_heclex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_HEC_E1700, "',' required after Iyy");
    return -1;
  }

  /* Izz */
  token = HECMW_heclex_next_token();
  if (token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
    Izz = HECMW_heclex_get_number();
  } else {
    set_err_token(token, HECMW_IO_HEC_E1700, "Izz reuiqred");
    return -1;
  }
  if (Izz <= 0) {
    set_err(HECMW_IO_HEC_E1709, "");
    return -1;
  }

  /* ',' */
  token = HECMW_heclex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_HEC_E1700, "',' required after Izz");
    return -1;
  }

  /* Jx */
  token = HECMW_heclex_next_token();
  if (token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
    Jx = HECMW_heclex_get_number();
  } else {
    set_err_token(token, HECMW_IO_HEC_E1700, "Jx reuiqred");
    return -1;
  }
  if (Jx <= 0) {
    set_err(HECMW_IO_HEC_E1710, "");
    return -1;
  }

  /* NL */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_NL) {
    set_err_token(token, HECMW_IO_HEC_E1700, "NL required after Jx");
    return -1;
  }

  /* set */
  sect_item->beam.vxyz[0] = nx;
  sect_item->beam.vxyz[1] = ny;
  sect_item->beam.vxyz[2] = nz;
  sect_item->beam.area    = area;
  sect_item->beam.Iyy     = Iyy;
  sect_item->beam.Izz     = Izz;
  sect_item->beam.Jx      = Jx;

  return 0;
}

static int read_section_interface(union hecmw_io_section_item *sect_item) {
  int token;
  double thickness;
  double gapcon  = 0.0;
  double gaprad1 = 0.0;
  double gaprad2 = 0.0;

  while (1) {
    /* THICKNESS */
    token = HECMW_heclex_next_token();
    if (token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
      thickness = HECMW_heclex_get_number();
    } else {
      set_err_token(token, HECMW_IO_HEC_E1700, "THICKNESS required");
      return -1;
    }
    if (thickness <= 0.0) {
      set_err(HECMW_IO_HEC_E1705, "");
      return -1;
    }

    /* ',' or NL */
    token = HECMW_heclex_next_token();
    if (token != ',' && token != HECMW_HECLEX_NL) {
      set_err_token(token, HECMW_IO_HEC_E1700,
                    "',' or NL reuqired after THICKNESS");
      return -1;
    }
    if (token == HECMW_HECLEX_NL) break;

    /* GAPCON */
    token = HECMW_heclex_next_token();
    if (token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
      gapcon = HECMW_heclex_get_number();
    } else if (token == ',') {
      HECMW_heclex_unput_token();
    } else if (token == HECMW_HECLEX_NL) {
      break;
    } else {
      set_err_token(token, HECMW_IO_HEC_E1700, "GAPCON reuiqred");
      return -1;
    }

    /* ',' or NL */
    token = HECMW_heclex_next_token();
    if (token != ',' && token != HECMW_HECLEX_NL) {
      set_err_token(token, HECMW_IO_HEC_E1700,
                    "',' or NL reuiqred after GAPCON");
      return -1;
    }
    if (token == HECMW_HECLEX_NL) break;

    /* GAPRAD1 */
    token = HECMW_heclex_next_token();
    if (token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
      gaprad1 = HECMW_heclex_get_number();
    } else if (token == ',') {
      HECMW_heclex_unput_token();
    } else if (token == HECMW_HECLEX_NL) {
      break;
    } else {
      set_err_token(token, HECMW_IO_HEC_E1700, "GAPRAD1 reuiqred");
      return -1;
    }

    /* ',' or NL */
    token = HECMW_heclex_next_token();
    if (token != ',' && token != HECMW_HECLEX_NL) {
      set_err_token(token, HECMW_IO_HEC_E1700,
                    "',' or NL reuqired after GAPRAD1");
      return -1;
    }
    if (token == HECMW_HECLEX_NL) break;

    /* GAPRAD2 */
    token = HECMW_heclex_next_token();
    if (token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
      gaprad2 = HECMW_heclex_get_number();
    } else if (token == HECMW_HECLEX_NL) {
      HECMW_heclex_unput_token();
    } else {
      set_err_token(token, HECMW_IO_HEC_E1700, "GAPRAD2 reuiqred");
      return -1;
    }

    /* NL */
    token = HECMW_heclex_next_token();
    if (token != HECMW_HECLEX_NL) {
      set_err_token(token, HECMW_IO_HEC_E1700, "NL required after GAPRAD2");
      return -1;
    }
    break;
  }

  /* set */
  sect_item->interface.thickness = thickness;
  sect_item->interface.gapcon    = gapcon;
  sect_item->interface.gaprad1   = gaprad1;
  sect_item->interface.gaprad2   = gaprad2;

  return 0;
}

static int read_section(void) {
  int token, state;
  int type      = -1;
  int secopt    = 0;
  int composite = -1;
  union hecmw_io_section_item sect_item;
  int flag_type                     = 0; /* flag for TYPE */
  int flag_egrp                     = 0; /* flag for EGRP */
  int flag_material                 = 0; /* flag for MATERIAL */
  int flag_composite                = 0; /* flag for COMPOSITE */
  int flag_secopt                   = 0; /* flag for SECOPT */
  int flag_input                    = 0; /* flag for INPUT */
  char egrp[HECMW_NAME_LEN + 1]     = "";
  char material[HECMW_NAME_LEN + 1] = "ALL";
  enum {
    ST_FINISHED,
    ST_HEADER_LINE,
    ST_HEADER_LINE_PARAM,
    ST_DATA_INCLUDE,
    ST_DATA_LINE_SOLID,
    ST_DATA_LINE_SHELL,
    ST_DATA_LINE_BEAM,
    ST_DATA_LINE_INTERFACE,
    ST_DATA_LINE_REGIST
  };

  state = ST_HEADER_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_HEADER_LINE) {
      if (read_section_head()) return -1;
      state = ST_HEADER_LINE_PARAM;
    } else if (state == ST_HEADER_LINE_PARAM) {
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_K_TYPE) {
        /* must */
        if (read_section_param_type(&type)) return -1;
        flag_type = 1;
      } else if (token == HECMW_HECLEX_K_EGRP) {
        /* must */
        if (read_section_param_egrp(egrp)) return -1;
        flag_egrp = 1;
      } else if (token == HECMW_HECLEX_K_MATERIAL) {
        /* optional */
        if (flag_composite) {
          set_err(HECMW_IO_HEC_E1703, "");
          return -1;
        }
        if (read_section_param_material(material)) return -1;
        flag_material = 1;
#if 0
			} else if(token == HECMW_HECLEX_K_COMPOSITE) {
				/* optional */
				if(flag_material) {
					set_err(HECMW_IO_HEC_E1703, "");
					return -1;
				}
				if(read_section_param_composite(&composite)) return -1;
				flag_composite = 1;
#endif
      } else if (token == HECMW_HECLEX_K_SECOPT) {
        /* optional */
        if (read_section_param_secopt(&secopt)) return -1;
        flag_secopt = 1;
      } else if (token == HECMW_HECLEX_K_INPUT) {
        /* optional */
        if (read_input(HECMW_IO_HEC_E1700)) return -1;
        flag_input = 1;
      } else {
        set_err_token(token, HECMW_IO_HEC_E1700, "Unknown parameter");
        return -1;
      }

      /* check next parameter */
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_NL) {
        /* check */
        if (!flag_type) {
          set_err(HECMW_IO_HEC_E1701, "");
          return -1;
        }
        if (!flag_egrp) {
          set_err(HECMW_IO_HEC_E1702, "");
          return -1;
        }
        /* set next state */
        if (flag_input) {
          state = ST_DATA_INCLUDE;
        } else if (type == HECMW_SECT_TYPE_SOLID) {
          state = ST_DATA_LINE_SOLID;
        } else if (type == HECMW_SECT_TYPE_SHELL) {
          state = ST_DATA_LINE_SHELL;
        } else if (type == HECMW_SECT_TYPE_BEAM) {
          state = ST_DATA_LINE_BEAM;
        } else if (type == HECMW_SECT_TYPE_INTERFACE) {
          state = ST_DATA_LINE_INTERFACE;
        } else {
          HECMW_assert(0);
        }
      } else if (token == ',') {
        ; /* continue this state */
      } else {
        set_err_token(token, HECMW_IO_HEC_E1700, "Unknown parameter");
        return -1;
      }
    } else if (state == ST_DATA_INCLUDE) {
      HECMW_assert(flag_input);

      if (HECMW_heclex_switch_to_include(include_filename)) return -1;

      /* set next state */
      if (type == HECMW_SECT_TYPE_SOLID) {
        state = ST_DATA_LINE_SOLID;
      } else if (type == HECMW_SECT_TYPE_SHELL) {
        state = ST_DATA_LINE_SHELL;
      } else if (type == HECMW_SECT_TYPE_BEAM) {
        state = ST_DATA_LINE_BEAM;
      } else if (type == HECMW_SECT_TYPE_INTERFACE) {
        state = ST_DATA_LINE_INTERFACE;
      } else {
        HECMW_assert(0);
      }
    } else if (state == ST_DATA_LINE_SOLID) {
      HECMW_assert(flag_egrp);
      HECMW_assert(flag_type);
      HECMW_assert(type == HECMW_SECT_TYPE_SOLID);

      if (read_section_solid(&sect_item)) return -1;
      state = ST_DATA_LINE_REGIST;
    } else if (state == ST_DATA_LINE_SHELL) {
      HECMW_assert(flag_egrp);
      HECMW_assert(flag_type);
      HECMW_assert(type == HECMW_SECT_TYPE_SHELL);

      if (read_section_shell(&sect_item)) return -1;
      state = ST_DATA_LINE_REGIST;
    } else if (state == ST_DATA_LINE_BEAM) {
      HECMW_assert(flag_egrp);
      HECMW_assert(flag_type);
      HECMW_assert(type == HECMW_SECT_TYPE_BEAM);

      if (read_section_beam(&sect_item)) return -1;
      state = ST_DATA_LINE_REGIST;
    } else if (state == ST_DATA_LINE_INTERFACE) {
      HECMW_assert(flag_egrp);
      HECMW_assert(flag_type);
      HECMW_assert(type == HECMW_SECT_TYPE_INTERFACE);

      if (read_section_interface(&sect_item)) return -1;
      state = ST_DATA_LINE_REGIST;
    } else if (state == ST_DATA_LINE_REGIST) {
      struct hecmw_io_section sect;
      HECMW_assert(flag_type);
      HECMW_assert(flag_egrp);

      /* set */
      strcpy(sect.egrp, egrp);
      strcpy(sect.material, material);
      sect.composite = composite;
      sect.secopt    = secopt;
      sect.type      = type;
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
  HECMW_log(HECMW_LOG_DEBUG, "read_section done");
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_sgrp_head(void) {
  int token;

  /* !SGROUP */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_H_SGROUP) {
    set_err_token(token, HECMW_IO_HEC_E1800, "!SGROUP required");
    return -1;
  }

  token = HECMW_heclex_next_token();
  if (token != ',') {
    set_err_token(token, HECMW_IO_HEC_E1800, "',' required after !SGROUP");
    return -1;
  }

  return 0;
}

static int read_sgrp_param_sgrp(char *sgrp) {
  int token;
  char *p;

  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E1800, "'=' required after SGRP");
    return -1;
  }
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_NAME) {
    set_err_token(token, HECMW_IO_HEC_E1800,
                  "SGRP must begin with a letter or '_'");
    return -1;
  }
  p = HECMW_heclex_get_text();
  if (strlen(p) > HECMW_NAME_LEN) {
    set_err(HECMW_IO_E0001, "");
    return -1;
  }
  strcpy(sgrp, p);
  HECMW_toupper(sgrp);
  if (HECMW_io_is_reserved_name(sgrp)) {
    set_err(HECMW_IO_E0003, "");
    return -1;
  }
  return 0;
}

static int read_sgrp_data(char *sgrp) {
  int i, n, *elem, *surf, token;
  struct hecmw_io_id *elem_head, *surf_head, *elem_prev, *surf_prev;
  struct hecmw_io_id *eid, *sid, *pe, *qe, *ps, *qs;

  n         = 0;
  elem_head = surf_head = NULL;
  elem_prev = surf_prev = NULL;
  while (1) {
    token = HECMW_heclex_next_token();
    if (n != 0 && token == HECMW_HECLEX_NL) break;

    eid = HECMW_malloc(sizeof(*eid));
    if (eid == NULL) {
      HECMW_set_error(errno, "");
      return -1;
    }
    eid->next = NULL;

    /* elemX */
    if (token != HECMW_HECLEX_INT) {
      set_err_token(token, HECMW_IO_HEC_E1800, "Element ID required");
      return -1;
    }
    eid->id = HECMW_heclex_get_number();

    if (elem_head == NULL) {
      elem_head = elem_prev = eid;
    } else {
      elem_prev->next = eid;
      elem_prev       = eid;
    }

    /* ',' */
    token = HECMW_heclex_next_token();
    if (token != ',') {
      set_err_token(token, HECMW_IO_HEC_E1800, "',' reuqired after element ID");
      return -1;
    }

    sid = HECMW_malloc(sizeof(*sid));
    if (sid == NULL) {
      HECMW_set_error(errno, "");
      return -1;
    }
    sid->next = NULL;

    /* lsufX */
    token = HECMW_heclex_next_token();
    if (token != HECMW_HECLEX_INT) {
      set_err_token(token, HECMW_IO_HEC_E1800, "Surface ID required");
      return -1;
    }
    sid->id = HECMW_heclex_get_number();

    if (surf_head == NULL) {
      surf_head = surf_prev = sid;
    } else {
      surf_prev->next = sid;
      surf_prev       = sid;
    }

    n++;

    /* ',' or NL */
    token = HECMW_heclex_next_token();
    if (token != ',' && token != HECMW_HECLEX_NL) {
      set_err_token(token, HECMW_IO_HEC_E1800,
                    "',' or NL required after surface ID");
      return -1;
    }
    if (token == HECMW_HECLEX_NL) break;
  }

  if (n > 0) {
    /* add elem and surf to group */
    elem = HECMW_malloc(sizeof(*elem) * n);
    if (elem == NULL) {
      HECMW_set_error(errno, "");
      return -1;
    }
    surf = HECMW_malloc(sizeof(*surf) * n);
    if (surf == NULL) {
      HECMW_set_error(errno, "");
      return -1;
    }
    i  = 0;
    qe = qs = NULL;
    pe      = elem_head;
    ps      = surf_head;
    for (i = 0; i < n; i++) {
      HECMW_assert(pe);
      HECMW_assert(ps);
      elem[i] = pe->id;
      surf[i] = ps->id;
      qe      = pe;
      qs      = ps;
      pe      = pe->next;
      ps      = ps->next;
      HECMW_free(qe);
      HECMW_free(qs);
    }
    if (HECMW_io_add_sgrp(sgrp, n, elem, surf) < 0) return -1;

    HECMW_free(elem);
    HECMW_free(surf);
  }

  return 0;
}

static int read_sgroup(void) {
  int token, state;
  int flag_sgrp                 = 0; /* flag for SGRP */
  int flag_input                = 0; /* flag for INPUT */
  char sgrp[HECMW_NAME_LEN + 1] = "";
  enum {
    ST_FINISHED,
    ST_HEADER_LINE,
    ST_HEADER_LINE_PARAM,
    ST_DATA_INCLUDE,
    ST_DATA_LINE
  };

  state = ST_HEADER_LINE;
  while (state != ST_FINISHED) {
    if (state == ST_HEADER_LINE) {
      if (read_sgrp_head()) return -1;
      state = ST_HEADER_LINE_PARAM;
    } else if (state == ST_HEADER_LINE_PARAM) {
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_K_SGRP) {
        /* must */
        if (read_sgrp_param_sgrp(sgrp)) return -1;
        flag_sgrp = 1;
      } else if (token == HECMW_HECLEX_K_INPUT) {
        /* oprtional */
        if (read_input(HECMW_IO_HEC_E1800)) return -1;
        flag_input = 1;
      } else {
        set_err_token(token, HECMW_IO_HEC_E1800, "Unknown parameter");
        return -1;
      }

      /* check next parameter */
      token = HECMW_heclex_next_token();
      if (token == HECMW_HECLEX_NL) {
        /* check */
        if (!flag_sgrp) {
          set_err(HECMW_IO_HEC_E1801, "");
          return -1;
        }
        state = flag_input ? ST_DATA_INCLUDE : ST_DATA_LINE;
      } else if (token == ',') {
        ; /* continue this state */
      } else {
        set_err_token(token, HECMW_IO_HEC_E1800, "Unknown parameter");
        return -1;
      }
    } else if (state == ST_DATA_INCLUDE) {
      HECMW_assert(flag_input);
      HECMW_assert(flag_sgrp);
      if (HECMW_heclex_switch_to_include(include_filename)) return -1;
      state = ST_DATA_LINE;
    } else if (state == ST_DATA_LINE) {
      HECMW_assert(flag_sgrp);
      if (read_sgrp_data(sgrp)) return -1;

      /* check next state */
      token = HECMW_heclex_next_token();
      if (token != HECMW_HECLEX_INT) {
        state = ST_FINISHED;
      } else {
        state = ST_DATA_LINE;
      }
      HECMW_heclex_unput_token();
    } else {
      HECMW_assert(0);
    }
  }
  HECMW_log(HECMW_LOG_DEBUG, "read_sgroup done");
  return 0;
}

/*----------------------------------------------------------------------------*/
#if 0
static int
read_system_head(void)
{
	int token;

	/* !SYSTEM */
	token = HECMW_heclex_next_token();
	if(token != HECMW_HECLEX_H_SYSTEM) {
		set_err_token(token, HECMW_IO_HEC_E1900, "!SYSTEM required");
		return -1;
	}

	/* NL */
	token = HECMW_heclex_next_token();
	if(token != HECMW_HECLEX_NL) {
		set_err_token(token, HECMW_IO_HEC_E1900, "NL required after !SYSTEM");
		return -1;
	}

	return 0;
}


static int
read_system_data_line1a(struct hecmw_system_param *system, int *last_token)
{
	int token;

	/* Xa */
	token = HECMW_heclex_next_token();
	if(token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
		system->xa = HECMW_heclex_get_number();
	} else if(token == ',') {
		system->xa = 0.0;
		HECMW_heclex_unput_token();
	} else {
		set_err_token(token, HECMW_IO_HEC_E1900, "Xa required");
		return -1;
	}

	/* ',' */
	token = HECMW_heclex_next_token();
	if(token != ',') {
		set_err_token(token, HECMW_IO_HEC_E1900, "',' required after Xa");
		return -1;
	}

	/* Ya */
	token = HECMW_heclex_next_token();
	if(token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
		system->ya = HECMW_heclex_get_number();
	} else if(token == ',') {
		system->ya = 0.0;
		HECMW_heclex_unput_token();
	} else {
		set_err_token(token, HECMW_IO_HEC_E1900, "Ya required");
		return -1;
	}

	/* ',' */
	token = HECMW_heclex_next_token();
	if(token != ',') {
		set_err_token(token, HECMW_IO_HEC_E1900, "',' required after Ya");
		return -1;
	}

	/* Za */
	token = HECMW_heclex_next_token();
	if(token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
		system->za = HECMW_heclex_get_number();
	} else if(token == ',' || token == HECMW_HECLEX_NL) {
		system->za = 0.0;
		HECMW_heclex_unput_token();
	} else {
		set_err_token(token, HECMW_IO_HEC_E1900, "Za required");
		return -1;
	}

	/* ',' or NL */
	token = HECMW_heclex_next_token();
	if(token != ',' && token != HECMW_HECLEX_NL) {
		set_err_token(token, HECMW_IO_HEC_E1900, "',' or NL required after Za");
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
	token = HECMW_heclex_next_token();
	if(token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
		system->xb = HECMW_heclex_get_number();
	} else if(token == ',') {
		system->xb = 0.0;
		HECMW_heclex_unput_token();
	} else {
		set_err_token(token, HECMW_IO_HEC_E1900, "Xb required");
		return -1;
	}

	/* ',' */
	token = HECMW_heclex_next_token();
	if(token != ',') {
		set_err_token(token, HECMW_IO_HEC_E1900, "',' required after Xb");
		return -1;
	}

	/* Yb */
	token = HECMW_heclex_next_token();
	if(token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
		system->yb = HECMW_heclex_get_number();
	} else if(token == ',') {
		system->yb = 0.0;
		HECMW_heclex_unput_token();
	} else {
		set_err_token(token, HECMW_IO_HEC_E1900, "Yb required");
		return -1;
	}

	/* ',' */
	token = HECMW_heclex_next_token();
	if(token != ',') {
		set_err_token(token, HECMW_IO_HEC_E1900, "',' required after Yb");
		return -1;
	}

	/* Zb */
	token = HECMW_heclex_next_token();
	if(token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
		system->zb = HECMW_heclex_get_number();
	} else if(token == HECMW_HECLEX_NL) {
		system->zb = 0.0;
		HECMW_heclex_unput_token();
	} else {
		set_err_token(token, HECMW_IO_HEC_E1900, "Zb required");
		return -1;
	}

	/*NL */
	token = HECMW_heclex_next_token();
	if(token != HECMW_HECLEX_NL) {
		set_err_token(token, HECMW_IO_HEC_E1900, "NL required after Zb");
		return -1;
	}

	return 0;
}


static int
read_system_data_line2(struct hecmw_system_param *system)
{
	int token;

	/* Xc */
	token = HECMW_heclex_next_token();
	if(token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
		system->xc = HECMW_heclex_get_number();
	} else if(token == ',') {
		system->xc = 0.0;
		HECMW_heclex_unput_token();
	} else {
		set_err_token(token, HECMW_IO_HEC_E1900, "Xc required");
		return -1;
	}

	/* ',' */
	token = HECMW_heclex_next_token();
	if(token != ',') {
		set_err_token(token, HECMW_IO_HEC_E1900, "',' required after Xc");
		return -1;
	}

	/* Yc */
	token = HECMW_heclex_next_token();
	if(token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
		system->yc = HECMW_heclex_get_number();
	} else if(token == ',') {
		system->yc = 0.0;
		HECMW_heclex_unput_token();
	} else {
		set_err_token(token, HECMW_IO_HEC_E1900, "Yc required");
		return -1;
	}

	/* ',' */
	token = HECMW_heclex_next_token();
	if(token != ',') {
		set_err_token(token, HECMW_IO_HEC_E1900, "',' required after Yc");
		return -1;
	}

	/* Zc */
	token = HECMW_heclex_next_token();
	if(token == HECMW_HECLEX_DOUBLE || token == HECMW_HECLEX_INT) {
		system->zc = HECMW_heclex_get_number();
	} else if(token == HECMW_HECLEX_NL) {
		system->zc = 0.0;
		HECMW_heclex_unput_token();
	} else {
		set_err_token(token, HECMW_IO_HEC_E1900, "Zc required");
		return -1;
	}

	/* NL */
	token = HECMW_heclex_next_token();
	if(token != HECMW_HECLEX_NL) {
		set_err_token(token, HECMW_IO_HEC_E1900, "NL required after Zc");
		return -1;
	}

	return 0;
}


static int
read_system(void)
{
	int token,state;
	struct hecmw_system_param *system;
	enum {
		ST_FINISHED,
		ST_HEADER_LINE,
		ST_DATA_LINE1,
		ST_DATA_LINE2
	};

	system = HECMW_malloc(sizeof(*system));
	if(system == NULL) {
		HECMW_set_error(errno, "");
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

	state = ST_HEADER_LINE;
	while(state != ST_FINISHED) {
		if(state == ST_HEADER_LINE) {
			if(read_system_head()) return -1;
			/* check next state */
			token = HECMW_heclex_next_token();
			if(token != HECMW_HECLEX_DOUBLE && token != HECMW_HECLEX_INT && token != ',') {
				/* clear !SYSTEM */
				HECMW_free(system);
				system = NULL;
				state = ST_FINISHED;
			} else {
				state = ST_DATA_LINE1;
			}
			HECMW_heclex_unput_token();
		} else if(state == ST_DATA_LINE1) {
			if(read_system_data_line1a(system, &token)) return -1;
			if(token == HECMW_HECLEX_NL) {
				state = ST_FINISHED;
				continue;
			}
			HECMW_assert(token == ',');

			if(read_system_data_line1b(system)) return -1;
			token = HECMW_heclex_next_token();
			if(token != HECMW_HECLEX_DOUBLE && token != HECMW_HECLEX_INT && token != ',') {
				state = ST_FINISHED;
			} else {
				state = ST_DATA_LINE2;
			}
			HECMW_heclex_unput_token();
		} else if(state == ST_DATA_LINE2) {
			if(read_system_data_line2(system)) return -1;
			state = ST_FINISHED;
		} else {
			HECMW_assert(0);
		}
	}

	/* set */
	HECMW_io_set_system(system);

	HECMW_log(HECMW_LOG_DEBUG, "read_system done");
	return 0;
}
#endif

/*----------------------------------------------------------------------------*/

static int read_zero(void) {
  int token;
  double zero;
  struct hecmw_io_zero *new_zero;

  /* !ZERO */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_H_ZERO) {
    set_err_token(token, HECMW_IO_HEC_E2000, "!ZERO required");
    return -1;
  }

  /* NL */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_NL) {
    set_err_token(token, HECMW_IO_HEC_E2000, "NL reqyured after !ZERO");
    return -1;
  }

  /* ZERO */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_INT && token != HECMW_HECLEX_DOUBLE) {
    set_err_token(token, HECMW_IO_HEC_E2000, "ZERO required");
    return -1;
  }
  zero = HECMW_heclex_get_number();

  /* NL */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_NL) {
    set_err_token(token, HECMW_IO_HEC_E2000, "NL reqyured after ZERO");
    return -1;
  }

  new_zero = HECMW_malloc(sizeof(*new_zero));
  if (new_zero == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }
  new_zero->zero = zero;

  /* set */
  HECMW_io_set_zero(new_zero);

  HECMW_log(HECMW_LOG_DEBUG, "read_zero done");
  return 0;
}

/*----------------------------------------------------------------------------*/

static int read_connectivity(void) {
  int token, type;

  /* !CONNECTIVITY */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_H_CONNECTIVITY) {
    set_err_token(token, HECMW_IO_HEC_E0200, "!CONNECTIVITY required");
    return -1;
  }

  /* , or NL */
  token = HECMW_heclex_next_token();
  if (token != ',' && token != HECMW_HECLEX_NL) {
    set_err_token(token, HECMW_IO_HEC_E0200,
                  "',' or NL reqyured after !CONNECTIVITY");
    return -1;
  }
  if (token == HECMW_HECLEX_NL) {
    connectivity_type = HECMW_CONNTYPE_HECMW; /* set default value */
    return 0;
  }

  /* TYPE */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_K_TYPE) {
    set_err_token(token, HECMW_IO_HEC_E0200, "TYPE required");
    return -1;
  }

  /* = */
  token = HECMW_heclex_next_token();
  if (token != '=') {
    set_err_token(token, HECMW_IO_HEC_E0200, "'=' reqyured after TYPE");
    return -1;
  }

  /* TYPE value */
  token = HECMW_heclex_next_token();
  switch (token) {
    case HECMW_HECLEX_K_HECMW:
      type = HECMW_CONNTYPE_HECMW;
      break;
    case HECMW_HECLEX_K_ABAQUS:
      type = HECMW_CONNTYPE_ABAQUS;
      break;
    case HECMW_HECLEX_K_NASTRAN:
      type = HECMW_CONNTYPE_NASTRAN;
      break;
    default:
      set_err_token(token, HECMW_IO_HEC_E0200, "Unsupported connectivity TYPE");
      return -1;
  }

  /* NL */
  token = HECMW_heclex_next_token();
  if (token != HECMW_HECLEX_NL) {
    set_err_token(token, HECMW_IO_HEC_E0200, "NL reqyured after TYPE value");
    return -1;
  }

  connectivity_type = type;

  HECMW_log(HECMW_LOG_DEBUG, "read_connectivity done");
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
    {HECMW_HECLEX_H_AMPLITUDE, read_amplitude},
    {HECMW_HECLEX_H_CONNECTIVITY, read_connectivity},
    {HECMW_HECLEX_H_CONTACT_PAIR, read_contact_pair},
    {HECMW_HECLEX_H_EMBED_PAIR, read_contact_pair},
    /*	{ HECMW_HECLEX_H_ECOPY,     read_ecopy     }, */
    /*	{ HECMW_HECLEX_H_EGEN,      read_egen      }, */
    {HECMW_HECLEX_H_EGROUP, read_egroup},
    {HECMW_HECLEX_H_ELEMENT, read_element},
    {HECMW_HECLEX_H_EQUATION, read_equation},
    {HECMW_HECLEX_H_HEADER, read_header},
    {HECMW_HECLEX_H_INCLUDE, read_include},
    {HECMW_HECLEX_H_INITIAL, read_initial},
    {HECMW_HECLEX_H_MATERIAL, read_material},
    /*	{ HECMW_HECLEX_H_NCOPY,     read_ncopy     }, */
    /*	{ HECMW_HECLEX_H_NFILL,     read_nfill     }, */
    /*	{ HECMW_HECLEX_H_NGEN,      read_ngen      }, */
    {HECMW_HECLEX_H_NGROUP, read_ngroup},
    {HECMW_HECLEX_H_NODE, read_node},
    {HECMW_HECLEX_H_SECTION, read_section},
    {HECMW_HECLEX_H_SGROUP, read_sgroup},
    /*	{ HECMW_HECLEX_H_SYSTEM,    read_system    }, */
    {HECMW_HECLEX_H_ZERO, read_zero},
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
  int token;
  ReadFunc func;

  while ((token = HECMW_heclex_next_token())) {
    if (token == HECMW_HECLEX_NL) continue;
    if (token == HECMW_HECLEX_H_END) {
      /* stop reading */
      return 0;
    }
    func = get_read_func(token);
    if (func == NULL) {
      char *p = HECMW_heclex_get_text();
      if (p[0] == '!') {
        set_err(HECMW_IO_HEC_E0099, "");
      } else {
        set_err(HECMW_IO_HEC_E0098, "");
      }
      return -1;
    }
    HECMW_heclex_unput_token(); /* unput !XXXX */
    if ((*func)()) return -1;
  }
  return 0;
}

/* read only. Not make hecmwST_local_mesh */
int HECMW_read_entire_mesh(const char *filename) {
  FILE *fp;

  HECMW_log(HECMW_LOG_DEBUG, "Start to read HECMW-ENTIRE mesh");

  if (filename == NULL) {
    HECMW_set_error(
        HECMW_IO_E0001,
        "Not specified filename for HECMW-ENTIRE mesh input routine");
    return -1;
  }
  HECMW_log(HECMW_LOG_DEBUG, "HECMW-ENTIRE mesh file is '%s'", filename);

  if (strlen(filename) > HECMW_FILENAME_LEN) {
    HECMW_set_error(HECMW_IO_E0002, "");
    return -1;
  }

  strcpy(grid_filename, filename);
  HECMW_io_set_gridfile(grid_filename);

  if ((fp = fopen(filename, "r")) == NULL) {
    HECMW_set_error(HECMW_IO_HEC_E0001, "File: %s, %s", filename,
                    strerror(errno));
    return -1;
  }

  if (HECMW_heclex_set_input(fp)) return -1;

  HECMW_log(HECMW_LOG_DEBUG, "Parsing...");
  if (parse()) {
    return -1;
  }

  if (fclose(fp)) {
    HECMW_set_error(HECMW_IO_HEC_E0002, "File: %s, %s", filename,
                    strerror(errno));
    return -1;
  }

  strcpy(grid_filename, "Unknown");

  return 0;
}

struct hecmwST_local_mesh *HECMW_get_entire_mesh(const char *filename) {
  struct hecmwST_local_mesh *local_mesh;

  if (HECMW_io_init()) return NULL;
  if (HECMW_io_pre_process()) return NULL;
  if (HECMW_read_entire_mesh(filename)) return NULL;
  if (HECMW_io_post_process()) return NULL;
  local_mesh = HECMW_io_make_local_mesh();
  if (local_mesh == NULL) return NULL;
  if (HECMW_io_finalize()) return NULL;

  strcpy(grid_filename, "Unknown");

  return local_mesh;
}
