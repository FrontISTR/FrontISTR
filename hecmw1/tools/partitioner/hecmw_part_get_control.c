/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

#include "hecmw_util.h"
#include "hecmw_io.h"

#include "hecmw_partlex.h"
#include "hecmw_part_define.h"
#include "hecmw_part_struct.h"
#include "hecmw_part_get_control.h"

static char ctrl_file_name[HECMW_FILENAME_LEN] = "\0";
static int  args_subdomain = 0;

/*================================================================================================*/

extern int HECMW_part_set_subdomains(int n_domain) {
  args_subdomain = n_domain;
  return 0;
}

extern int HECMW_part_set_ctrl_file_name(char *fname) {
  if (fname == NULL) {
    HECMW_set_error(HECMW_PART_E_INV_ARG, "'fname' is NULL");
    goto error;
  }
  if (strlen(fname) > HECMW_FILENAME_LEN) {
    HECMW_set_error(HECMW_PART_E_TOO_LONG_FNAME,
                    "control file for partitioner");
    goto error;
  }

  strcpy(ctrl_file_name, fname);

  return 0;

error:
  return -1;
}

/*============================================================================*/
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  partitioning type < TYPE={ NODE-BASED | ELEMENT-BASED } >                 */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int part_cont_type(void) {
  int token;

  /* '=' */
  token = HECMW_partlex_next_token();
  if (token != '=') {
    HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_TYPE_NOEQ));
    HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
              "part_cont_type", HECMW_partlex_get_text());
    return -1;
  }

  /* { NODE-BASED | ELEMENT-BASED } */
  token = HECMW_partlex_next_token();
  switch (token) {
    case HECMW_PARTLEX_V_NODE_BASED: /* TYPE=NODE-BASED */
      return HECMW_PART_TYPE_NODE_BASED;

    case HECMW_PARTLEX_V_ELEMENT_BASED: /* TYPE=ELEMENT-BASED */
      return HECMW_PART_TYPE_ELEMENT_BASED;

    default:
      HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_TYPE_INVAL));
      HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
                "part_cont_type", HECMW_partlex_get_text());
      return -1;
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  partitioning method < METHOD={ RCB | KMETIS | PMETIS | USER } >                  */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int part_cont_method(void) {
  int token;

  /* '=' */
  token = HECMW_partlex_next_token();
  if (token != '=') {
    HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_METHOD_NOEQ));
    HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
              "part_cont_method", HECMW_partlex_get_text());
    return -1;
  }

  /* { RCB | KMETIS | PMETIS | USER } */
  token = HECMW_partlex_next_token();
  switch (token) {
    case HECMW_PARTLEX_V_RCB: /* METHOD=RCB */
      return HECMW_PART_METHOD_RCB;

    case HECMW_PARTLEX_V_KMETIS: /* METHOD=KMETIS */
      return HECMW_PART_METHOD_KMETIS;

    case HECMW_PARTLEX_V_PMETIS: /* METHOD=PMETIS */
      return HECMW_PART_METHOD_PMETIS;

    case HECMW_PARTLEX_V_USER: /* METHOD=USER */
      return HECMW_PART_METHOD_USER;

    default:
      HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_METHOD_INVAL));
      HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
                "part_cont_method", HECMW_partlex_get_text());
      return -1;
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  number of sub-domains < DOMAIN=n >                                        */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int part_cont_domain(void) {
  int token, domain;

  /* '=' */
  token = HECMW_partlex_next_token();
  if (token != '=') {
    HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_DOMAIN_NOEQ));
    HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
              "part_cont_domain", HECMW_partlex_get_text());
    return -1;
  }

  /* n */
  token = HECMW_partlex_next_token();
  switch (token) {
    case HECMW_PARTLEX_INT: /* DOMAIN=n */
      domain = (int)HECMW_partlex_get_number();
      if (domain < 1) {
        HECMW_log(HECMW_LOG_ERROR,
                  HECMW_strmsg(HECMW_PART_E_CTRL_DOMAIN_INVAL));
        HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
                  "part_cont_domain", HECMW_partlex_get_text());
        return -1;
      }
      return domain;

    default:
      HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_DOMAIN_INVAL));
      HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
                "part_cont_domain", HECMW_partlex_get_text());
      return -1;
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  depth of overlapping zone < DEPTH=n >                                     */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int part_cont_depth(void) {
  int token, depth;

  /* '=' */
  token = HECMW_partlex_next_token();
  if (token != '=') {
    HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_DEPTH_NOEQ));
    HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
              "part_cont_depth", HECMW_partlex_get_text());
    return -1;
  }

  /* n */
  token = HECMW_partlex_next_token();
  switch (token) {
    case HECMW_PARTLEX_INT: /* DEPTH=n */
      depth = (int)HECMW_partlex_get_number();
      if (depth < 1) {
        HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_DEPTH_INVAL));
        HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
                  "part_cont_depth", HECMW_partlex_get_text());
        return -1;
      }
      return depth;

    default:
      HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_DEPTH_INVAL));
      HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
                "part_cont_depth", HECMW_partlex_get_text());
      return -1;
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  ucd file name < UCD=filename >                                            */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int part_cont_ucd(char *name) {
  char *p;
  int token, is_print_ucd = 0;

  /* '=' */
  token = HECMW_partlex_next_token();
  if (token != '=') {
    HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_UCD_NOEQ));
    HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
              "part_cont_ucd", HECMW_partlex_get_text());
    return -1;
  }

  /* filename */
  token = HECMW_partlex_next_token();
  switch (token) {
    case HECMW_PARTLEX_NAME:     /* UCD=filename */
    case HECMW_PARTLEX_FILENAME: /* UCD=filename */
      p = HECMW_partlex_get_text();
      if (strlen(p) > HECMW_FILENAME_LEN) {
        HECMW_log(HECMW_LOG_ERROR,
                  HECMW_strmsg(HECMW_PART_E_CTRL_UCD_TOO_LONG));
        HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
                  "part_cont_ucd", HECMW_partlex_get_text());
        return -1;
      }
      strcpy(name, p);
      is_print_ucd = 1;

      return is_print_ucd;

    default:
      HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_UCD_INVAL));
      HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
                "part_cont_ucd", HECMW_partlex_get_text());
      return -1;
  }

  HECMW_assert(0);

  return -1;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  partitioning contact < CONTACT={ DEFAULT | AGGREGATE | DISTRIBUTE | SIMPLE }
 * > */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int part_cont_contact(void) {
  int token;

  /* '=' */
  token = HECMW_partlex_next_token();
  if (token != '=') {
    HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_METHOD_NOEQ));
    HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
              "part_cont_method", HECMW_partlex_get_text());
    return -1;
  }

  /* { DEFAULT | AGGREGATE | DISTRIBUTE | SIMPLE } */
  token = HECMW_partlex_next_token();
  switch (token) {
    case HECMW_PARTLEX_V_DEFAULT: /* CONTACT=DEFAULT */
      return HECMW_PART_CONTACT_DEFAULT;

    case HECMW_PARTLEX_V_AGGREGATE: /* CONTACT=AGGREGATE */
      return HECMW_PART_CONTACT_AGGREGATE;

    case HECMW_PARTLEX_V_DISTRIBUTE: /* CONTACT=DISTRIBUTE */
      return HECMW_PART_CONTACT_DISTRIBUTE;

    case HECMW_PARTLEX_V_SIMPLE: /* CONTACT=SIMPLE */
      return HECMW_PART_CONTACT_SIMPLE;

    default:
      HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_CONTACT_INVAL));
      HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
                "part_cont_method", HECMW_partlex_get_text());
      return -1;
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  part file name < PART=filename >                                            */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int part_cont_part(char *name) {
  char *p;
  int token, is_print_part = 0;

  /* '=' */
  token = HECMW_partlex_next_token();
  if (token != '=') {
    HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_PART_NOEQ));
    HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
              "part_cont_part", HECMW_partlex_get_text());
    return -1;
  }

  /* filename */
  token = HECMW_partlex_next_token();
  switch (token) {
    case HECMW_PARTLEX_NAME:     /* PART=filename */
    case HECMW_PARTLEX_FILENAME: /* PART=filename */
      p = HECMW_partlex_get_text();
      if (strlen(p) > HECMW_FILENAME_LEN) {
        HECMW_log(HECMW_LOG_ERROR,
                  HECMW_strmsg(HECMW_PART_E_CTRL_PART_TOO_LONG));
        HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
                  "part_cont_part", HECMW_partlex_get_text());
        return -1;
      }
      strcpy(name, p);
      is_print_part = 1;

      return is_print_part;

    default:
      HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_PART_INVAL));
      HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
                "part_cont_part", HECMW_partlex_get_text());
      return -1;
  }

  HECMW_assert(0);

  return -1;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  partitioning directions for RCB partitioning < x, y, z >                  */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int part_cont_rcb_divs(int n_domain) {
  int base, rest, n;
  int i;

  for (i = 0; i < n_domain; i++) {
    n    = (int)pow(2, i);
    base = n_domain / n;
    rest = n_domain % n;

    if ((base == 1) && (rest == 0)) return i;
    if ((base == 0) && (rest == 1)) return 0;
    if ((base == 1) && (rest > 0)) return -1;
  }

  HECMW_assert(0);

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int part_cont_rcb_opt(struct hecmw_part_cont_data *cont_data) {
  int token;
  int i;

  if (cont_data->n_domain == 1) {
    cont_data->n_rcb_div = 0;
    cont_data->rcb_axis  = NULL;
    return 0;
  }

  /* number of decompositions */
  cont_data->n_rcb_div = part_cont_rcb_divs(cont_data->n_domain);
  if (cont_data->n_rcb_div < -1) {
    HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_DOMAIN_POW));
    HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%d)", __FILE__, __LINE__,
              "part_cont_rcb_opt", cont_data->n_domain);
    return -1;
  }

  /* partitioning directions */
  cont_data->rcb_axis = (int *)HECMW_malloc(sizeof(int) * cont_data->n_rcb_div);
  if (cont_data->rcb_axis == NULL) {
    HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(errno));
    HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
              "part_cont_rcb_opt", "cont_data->rcb_axis");
    return -1;
  }

  for (i = 0; i < cont_data->n_rcb_div; i++) {
    token = HECMW_partlex_next_token();

    switch (token) {
      case 'x':
        cont_data->rcb_axis[i] = HECMW_PART_RCB_X_AXIS;
        break;

      case 'y':
        cont_data->rcb_axis[i] = HECMW_PART_RCB_Y_AXIS;
        break;

      case 'z':
        cont_data->rcb_axis[i] = HECMW_PART_RCB_Z_AXIS;
        break;

      default:
        HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_RCB_INVAL));
        HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
                  "part_cont_rcb_opt", HECMW_partlex_get_text());
        HECMW_free(cont_data->rcb_axis);
        return -1;
    }

    token = HECMW_partlex_next_token();
    if (token != ',' && token != HECMW_PARTLEX_NL && token) {
      HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
                "part_cont_rcb_opt", HECMW_partlex_get_text());
      HECMW_free(cont_data->rcb_axis);
      return -1;
    }

    if (i + 1 == cont_data->n_rcb_div) {
      if (token == ',') {
        HECMW_log(HECMW_LOG_WARN, HECMW_strmsg(HECMW_PART_W_CTRL_RCB_MANY_DIR));
        while ((token = HECMW_partlex_next_token()) != HECMW_PARTLEX_NL) {
        }
      }
      break;
    } else {
      if (token == HECMW_PARTLEX_NL || !token) {
        HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_RCB_FEW_DIR));
        HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s", __FILE__, __LINE__,
                  "part_cont_rcb_opt");
        HECMW_free(cont_data->rcb_axis);
        return -1;
      }
    }
  }

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  control data for partitioner < !PARTITION >                               */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int part_cont_partition(struct hecmw_part_cont_data *cont_data) {
  int token;

  cont_data->type         = -1;
  cont_data->method       = -1;
  cont_data->n_domain     = -1;
  cont_data->depth        = -1;
  cont_data->is_print_ucd = -1;
  cont_data->contact      = -1;
  cont_data->is_print_part= -1;
  strcpy(cont_data->ucd_file_name, "\0");
  strcpy(cont_data->part_file_name, "\0");

  while ((token = HECMW_partlex_next_token()) != HECMW_PARTLEX_NL ||
         (token = HECMW_partlex_next_token())) {
    switch (token) {
      case HECMW_PARTLEX_K_TYPE: /* TYPE */
        cont_data->type = part_cont_type();
        if (cont_data->type < 0) return -1;
        break;

      case HECMW_PARTLEX_K_METHOD: /* METHOD */
        cont_data->method = part_cont_method();
        if (cont_data->method < 0) return -1;
        break;

      case HECMW_PARTLEX_K_DOMAIN: /* DOMAIN */
        cont_data->n_domain = part_cont_domain();
        if (cont_data->n_domain < 0) return -1;
        break;

      case HECMW_PARTLEX_K_DEPTH: /* DEPTH */
        cont_data->depth = part_cont_depth();
        if (cont_data->depth < 0) return -1;
        break;

      case HECMW_PARTLEX_K_UCD: /* UCD */
        cont_data->is_print_ucd = part_cont_ucd(cont_data->ucd_file_name);
        if (cont_data->is_print_ucd < 0) return -1;
        break;

      case HECMW_PARTLEX_K_CONTACT: /* CONTACT */
        cont_data->contact = part_cont_contact();
        if (cont_data->contact < 0) return -1;
        break;

      case HECMW_PARTLEX_K_PART: /* PART */
        cont_data->is_print_part = part_cont_part(cont_data->part_file_name);
        if (cont_data->is_print_part < 0) return -1;
        break;

      default:
        HECMW_log(HECMW_LOG_ERROR, "%s %s (%s)",
                  HECMW_strmsg(HECMW_PART_E_INVALID_TOKEN),
                  "in control data for partitioner", HECMW_partlex_get_text());
        HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
                  "part_cont_partition", HECMW_partlex_get_text());
        return -1;
    }

    token = HECMW_partlex_next_token();
    if (token != ',' && token != HECMW_PARTLEX_NL && token != EOF) {
      HECMW_log(HECMW_LOG_ERROR, "%s %s (%s)",
                HECMW_strmsg(HECMW_PART_E_INVALID_TOKEN),
                "in control file for partitioner", HECMW_partlex_get_text());
      HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
                "part_cont_partition", HECMW_partlex_get_text());
      return -1;
    } else {
      if (token == HECMW_PARTLEX_NL || token == EOF) break;
    }
  }

  /* check data */
  if (cont_data->type < 0) {
    HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_NO_TYPE));
    HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s", __FILE__, __LINE__,
              "part_cont_partition");
    return -1;
  }

#ifndef HECMW_PART_WITH_METIS
  if (cont_data->method == HECMW_PART_METHOD_PMETIS) {
    HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_NODEF_PMETIS));
    HECMW_log(HECMW_LOG_DEBUG, "%s:%d", __FILE__, __LINE__);
    return -1;
  }

  if (cont_data->method == HECMW_PART_METHOD_KMETIS) {
    HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_NODEF_KMETIS));
    HECMW_log(HECMW_LOG_DEBUG, "%s:%d", __FILE__, __LINE__);
    return -1;
  }
#endif

  HECMW_assert(cont_data->type == HECMW_PART_TYPE_NODE_BASED ||
               cont_data->type == HECMW_PART_TYPE_ELEMENT_BASED);

  if (cont_data->method < 0) {
    HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_NO_METHOD));
    HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s", __FILE__, __LINE__,
              "part_cont_partition");
    return -1;
  }
  HECMW_assert(cont_data->method == HECMW_PART_METHOD_RCB ||
               cont_data->method == HECMW_PART_METHOD_KMETIS ||
               cont_data->method == HECMW_PART_METHOD_PMETIS);

  if (cont_data->n_domain <= 0) {
    HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_NO_DOMAIN));
    HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s", __FILE__, __LINE__,
              "part_cont_partition");
    return -1;
  }

  if (cont_data->depth < 0) {
    cont_data->depth = 1;
  }

  if (cont_data->is_print_ucd < 0) {
    cont_data->is_print_ucd = 0;
  }

  if (cont_data->is_print_part < 0) {
    cont_data->is_print_part = 0;
  }

  if (cont_data->method == HECMW_PART_METHOD_USER) {
    cont_data->is_print_part = 0;
  }

  /* partitioning directions ( option for RCB ) */
  if (cont_data->method == HECMW_PART_METHOD_RCB) {
    if (token == EOF) {
      HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_RCB_NODIR));
      HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s", __FILE__, __LINE__,
                "part_cont_partition");
      return -1;
    }
    part_cont_rcb_opt(cont_data);
  }

  if (token == EOF) {
    HECMW_log(HECMW_LOG_ERROR, "%s (%s)",
              HECMW_strmsg(HECMW_PART_E_INVALID_EOF),
              "control file for partitioner");
    HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s", __FILE__, __LINE__,
              "part_cont_partition");
    return -1;
  } else {
    return 0;
  }
}

/*----------------------------------------------------------------------------*/
/*  get control data                                                          */
/*----------------------------------------------------------------------------*/
static int part_get_control(struct hecmw_part_cont_data *cont_data) {
  int token, stat;
  FILE *fp;

  if (args_subdomain){
    cont_data->type         = HECMW_PART_TYPE_NODE_BASED;
    cont_data->method       = HECMW_PART_METHOD_KMETIS; /*HECMW_PART_METHOD_RCB; HECMW_PART_METHOD_PMETIS;*/
    cont_data->n_domain     = args_subdomain;
    cont_data->depth        = 1;
    cont_data->is_print_ucd = 0;
    cont_data->contact      = HECMW_PART_CONTACT_DEFAULT;
    cont_data->is_print_part= 0;
    return 0;
  }

  /* open file */
  if (strlen(ctrl_file_name) == 0) {
    HECMW_log(HECMW_LOG_ERROR, "%s (%s)",
              HECMW_strmsg(HECMW_PART_E_NULL_POINTER),
              "control file name for partitioner is not set");
    HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s", __FILE__, __LINE__,
              "part_get_control");
    return -1;
  }
  if ((fp = fopen(ctrl_file_name, "r")) == NULL) {
    HECMW_log(HECMW_LOG_ERROR, "%s (%s)", HECMW_strmsg(errno),
              "control file for partitioner");
    HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s", __FILE__, __LINE__,
              "part_get_control");
    return -1;
  }

  /* initialize lex */
  stat = HECMW_partlex_set_input(fp);
  if (stat) {
    return -1;
  }

  /* get control data */
  while ((token = HECMW_partlex_next_token())) {
    switch (token) {
      case HECMW_PARTLEX_H_PARTITION: /* !PARTITION */
        token = HECMW_partlex_next_token();
        switch (token) {
          case ',': /* normal */
            if (part_cont_partition(cont_data)) return -1;
            break;

          case HECMW_PARTLEX_NL: /* no option */
            HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_CTRL_NO_TYPE));
            HECMW_log(HECMW_LOG_ERROR,
                      HECMW_strmsg(HECMW_PART_E_CTRL_NO_METHOD));
            HECMW_log(HECMW_LOG_ERROR,
                      HECMW_strmsg(HECMW_PART_E_CTRL_NO_DOMAIN));
            HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s", __FILE__, __LINE__,
                      "part_cont_partition");
            return -1;

          default: /* invalid delimiter */
            HECMW_log(HECMW_LOG_ERROR, "%s (%s)",
                      HECMW_strmsg(HECMW_PART_E_INVALID_TOKEN),
                      HECMW_partlex_get_text());
            HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
                      "part_cont_partition", HECMW_partlex_get_text());
            return -1;
        }
        break;

      case HECMW_PARTLEX_NL: /* new line */
        break;

      case 'x':
      case 'y':
      case 'z':
        switch (cont_data->method) {
          case HECMW_PART_METHOD_RCB:
            if (cont_data->n_domain > 1) {
              HECMW_log(HECMW_LOG_ERROR,
                        HECMW_strmsg(HECMW_PART_E_INVALID_TOKEN));
              HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
                        "part_cont_partition", HECMW_partlex_get_text());
              return -1;
            } else {
              break;
            }

          case HECMW_PART_METHOD_KMETIS:
          case HECMW_PART_METHOD_PMETIS:
          case HECMW_PART_METHOD_USER:
            HECMW_log(HECMW_LOG_WARN,
                      HECMW_strmsg(HECMW_PART_W_CTRL_DIR_WORCB));
            break;

          default:
            HECMW_log(HECMW_LOG_ERROR,
                      HECMW_strmsg(HECMW_PART_E_INVALID_TOKEN));
            HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
                      "part_cont_partition", HECMW_partlex_get_text());
            return -1;
        }
        while ((token = HECMW_partlex_next_token()) &&
               token != HECMW_PARTLEX_NL) {
        }
        goto finalize;

      default:
        HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(HECMW_PART_E_INVALID_TOKEN));
        HECMW_log(HECMW_LOG_DEBUG, "%s:%d:%s (%s)", __FILE__, __LINE__,
                  "part_cont_partition", HECMW_partlex_get_text());
        return -1;
    }
  }

/* close file */
finalize:
  if (fclose(fp)) {
    HECMW_log(HECMW_LOG_ERROR, "%s (%s)", HECMW_strmsg(HECMW_PART_E_FILE_CLOSE),
              "control file for partitioner");
    return -1;
  }

  return 0;
}

/*============================================================================*/

extern struct hecmw_part_cont_data *HECMW_part_get_control(void) {
  struct hecmw_part_cont_data *cont_data;

  /* allocate structure for control data */
  cont_data = (struct hecmw_part_cont_data *)HECMW_malloc(
      sizeof(struct hecmw_part_cont_data));
  if (cont_data == NULL) {
    HECMW_log(HECMW_LOG_ERROR, HECMW_strmsg(errno));
    return NULL;
  }
  cont_data->rcb_axis = NULL;

  /* get control data via file */
  if (part_get_control(cont_data)) {
    HECMW_part_free_control(cont_data);
    return NULL;
  }

  return cont_data;
}

/*================================================================================================*/

extern void HECMW_part_free_control(struct hecmw_part_cont_data *cont_data) {
  if (cont_data->rcb_axis) {
    HECMW_free(cont_data->rcb_axis);
  }
  HECMW_free(cont_data);
}
