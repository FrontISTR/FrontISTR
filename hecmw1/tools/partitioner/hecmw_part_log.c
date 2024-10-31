/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "hecmw_msgno.h"
#include "hecmw_error.h"
#include "hecmw_malloc.h"
#include "hecmw_part_define.h"
#include "hecmw_part_log.h"

static int is_init = 0;

static int n_node_g = 0;

static int n_elem_g = 0;

static int *n_node = NULL;

static int *nn_internal = NULL;

static int *n_elem = NULL;

static int *ne_internal = NULL;

static int *n_neighbor_pe = NULL;

static int n_domain = 0;

static int depth = 0;

static long long int n_edge = 0;

static int n_edgecut = 0;

static char part_type[HECMW_NAME_LEN] = "";

static char part_method[HECMW_NAME_LEN] = "";

static char part_contact[HECMW_NAME_LEN] = "";

/*================================================================================================*/

static void clean_log(void) {
  HECMW_free(n_node);
  HECMW_free(n_elem);
  HECMW_free(nn_internal);
  HECMW_free(ne_internal);
  HECMW_free(n_neighbor_pe);

  n_node      = NULL;
  n_elem      = NULL;
  nn_internal = NULL;
  ne_internal = NULL;

  is_init   = 0;
  n_node_g  = 0;
  n_elem_g  = 0;
  n_domain  = 0;
  depth     = 0;
  n_edge    = 0;
  n_edgecut = 0;
  memset(part_type, '\0', HECMW_NAME_LEN);
  memset(part_method, '\0', HECMW_NAME_LEN);
}

extern int HECMW_part_init_log(int _n_domain) {
  if (is_init) {
    HECMW_log(HECMW_LOG_WARN, HECMW_strmsg(HECMW_PART_W_LOG_INIT_ALREADY));
    return 0;
  } else {
    is_init = 1;
  }

  if (_n_domain < 1) {
    HECMW_set_error(HECMW_PART_E_INVALID_NDOMAIN, "%d", _n_domain);
    goto error;
  }

  n_domain = _n_domain;

  n_node = (int *)HECMW_calloc(n_domain, sizeof(int));
  if (n_node == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  n_elem = (int *)HECMW_calloc(n_domain, sizeof(int));
  if (n_elem == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  nn_internal = (int *)HECMW_calloc(n_domain, sizeof(int));
  if (nn_internal == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  ne_internal = (int *)HECMW_calloc(n_domain, sizeof(int));
  if (ne_internal == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  n_neighbor_pe = (int *)HECMW_calloc(n_domain, sizeof(int));
  if (n_neighbor_pe == NULL) {
    HECMW_set_error(errno, "");
    goto error;
  }

  return 0;

error:
  return -1;
}

/*================================================================================================*/

extern int HECMW_part_set_log_part_type(int _part_type) {
  if (is_init == 0) {
    HECMW_set_error(HECMW_PART_E_LOG_INIT_NOT_YET, "");
    goto error;
  }

  switch (_part_type) {
    case HECMW_PART_TYPE_NODE_BASED:
      strcpy(part_type, "NODE-BASED");
      break;

    case HECMW_PART_TYPE_ELEMENT_BASED:
      strcpy(part_type, "ELEMENT-BASED");
      break;

    default:
      HECMW_set_error(HECMW_PART_E_INVALID_PTYPE, "");
      goto error;
  }

  return 0;

error:
  return -1;
}

extern int HECMW_part_set_log_part_method(int _part_method) {
  if (is_init == 0) {
    HECMW_set_error(HECMW_PART_E_LOG_INIT_NOT_YET, "");
    return -1;
  }

  switch (_part_method) {
    case HECMW_PART_METHOD_RCB:
      strcpy(part_method, "RCB");
      break;

    case HECMW_PART_METHOD_KMETIS:
      strcpy(part_method, "kMETIS");
      break;

    case HECMW_PART_METHOD_PMETIS:
      strcpy(part_method, "pMETIS");
      break;

    case HECMW_PART_METHOD_USER:
      strcpy(part_method, "USER");
      break;

    default:
      HECMW_set_error(HECMW_PART_E_INVALID_PMETHOD, "");
      return -1;
  }

  return 0;
}

extern int HECMW_part_set_log_part_depth(int _depth) {
  if (is_init == 0) {
    HECMW_set_error(HECMW_PART_E_LOG_INIT_NOT_YET, "");
    return -1;
  }

  if (_depth < 1) {
    HECMW_set_error(HECMW_PART_E_INVALID_PDEPTH, "");
    return -1;
  }

  depth = _depth;

  return 0;
}

extern int HECMW_part_set_log_part_contact(int _part_contact) {
  if (is_init == 0) {
    HECMW_set_error(HECMW_PART_E_LOG_INIT_NOT_YET, "");
    return -1;
  }

  switch (_part_contact) {
    case HECMW_PART_CONTACT_AGGREGATE:
      strcpy(part_contact, "AGGREGATE");
      break;

    case HECMW_PART_CONTACT_DISTRIBUTE:
      strcpy(part_contact, "DISTRIBUTE");
      break;

    case HECMW_PART_CONTACT_SIMPLE:
      strcpy(part_contact, "SIMPLE");
      break;

    case HECMW_PART_CONTACT_DEFAULT:
      strcpy(part_contact, "DEFAULT");
      break;

    default:
      strcpy(part_contact, "not set");
      break;
  }

  return 0;
}

extern int HECMW_part_set_log_n_edgecut(long long int _n_edge, int _n_edgecut) {
  if (is_init == 0) {
    HECMW_set_error(HECMW_PART_E_LOG_INIT_NOT_YET, "");
    return -1;
  }

  if (_n_edgecut < 0) {
    HECMW_set_error(HECMW_PART_E_NEDGECUT_LOWER, "%d", _n_edgecut);
    return -1;
  }

  if (_n_edge < 1) {
    HECMW_set_error(HECMW_PART_E_NEDGECUTA_LOWER, "%lld", _n_edge);
    return -1;
  }

  n_edge    = _n_edge;
  n_edgecut = _n_edgecut;

  return 0;
}

extern int HECMW_part_set_log_n_node_g(int _n_node_g) {
  if (is_init == 0) {
    HECMW_set_error(HECMW_PART_E_LOG_INIT_NOT_YET, "");
    return -1;
  }

  if (_n_node_g < 1) {
    HECMW_set_error(HECMW_PART_E_NNODE_MIN, "%d", _n_node_g);
    return -1;
  }

  n_node_g = _n_node_g;

  return 0;
}

extern int HECMW_part_set_log_n_elem_g(int _n_elem_g) {
  if (is_init == 0) {
    HECMW_set_error(HECMW_PART_E_LOG_INIT_NOT_YET, "");
    return -1;
  }

  if (_n_elem_g < 1) {
    HECMW_set_error(HECMW_PART_E_NELEM_MIN, "%d", _n_elem_g);
    return -1;
  }

  n_elem_g = _n_elem_g;

  return 0;
}

extern int HECMW_part_set_log_n_node(int domain, int _n_node) {
  if (is_init == 0) {
    HECMW_set_error(HECMW_PART_E_LOG_INIT_NOT_YET, "");
    return -1;
  }

  if (domain < 0) {
    HECMW_set_error(HECMW_PART_E_DOMAIN_MIN, "domain");
    return -1;
  }
  if (domain >= n_domain) {
    HECMW_set_error(HECMW_PART_E_DOMAIN_MAX, "domain");
    return -1;
  }

  if (_n_node < 1) {
    HECMW_set_error(HECMW_PART_E_NNODE_MIN, "_n_node");
    return -1;
  }

  n_node[domain] = _n_node;

  return 0;
}

extern int HECMW_part_set_log_n_elem(int domain, int _n_elem) {
  if (is_init == 0) {
    HECMW_set_error(HECMW_PART_E_LOG_INIT_NOT_YET, "");
    return -1;
  }

  if (domain < 0) {
    HECMW_set_error(HECMW_PART_E_DOMAIN_MIN, "%d", domain);
    return -1;
  }
  if (domain >= n_domain) {
    HECMW_set_error(HECMW_PART_E_DOMAIN_MAX, "%d", domain);
    return -1;
  }

  if (_n_elem < 1) {
    HECMW_set_error(HECMW_PART_E_NELEM_MIN, "%d", _n_elem);
    return -1;
  }

  n_elem[domain] = _n_elem;

  return 0;
}

extern int HECMW_part_set_log_nn_internal(int domain, int _nn_internal) {
  if (is_init == 0) {
    HECMW_set_error(HECMW_PART_E_LOG_INIT_NOT_YET, "");
    return -1;
  }

  if (domain < 0) {
    HECMW_set_error(HECMW_PART_E_DOMAIN_MIN, "%d", domain);
    return -1;
  }
  if (domain >= n_domain) {
    HECMW_set_error(HECMW_PART_E_DOMAIN_MAX, "%d", domain);
    return -1;
  }

  if (_nn_internal < 0) {
    HECMW_set_error(HECMW_PART_E_NNINT_MIN, "%d", domain);
    return -1;
  }
  if (_nn_internal > n_node[domain]) {
    HECMW_set_error(HECMW_PART_E_NNINT_MAX, "%d", domain);
    return -1;
  }

  nn_internal[domain] = _nn_internal;

  return 0;
}

extern int HECMW_part_set_log_ne_internal(int domain, int _ne_internal) {
  if (is_init == 0) {
    HECMW_set_error(HECMW_PART_E_LOG_INIT_NOT_YET, "");
    return -1;
  }

  if (domain < 0) {
    HECMW_set_error(HECMW_PART_E_DOMAIN_MIN, "%d", domain);
    return -1;
  }
  if (domain >= n_domain) {
    HECMW_set_error(HECMW_PART_E_DOMAIN_MAX, "%d", domain);
    return -1;
  }

  if (_ne_internal < 0) {
    HECMW_set_error(HECMW_PART_E_NEINT_MIN, "%d", _ne_internal);
    return -1;
  }

  if (_ne_internal > n_elem[domain]) {
    HECMW_set_error(HECMW_PART_E_NEINT_MAX, "%d", _ne_internal);
    return -1;
  }

  ne_internal[domain] = _ne_internal;

  return 0;
}

extern int HECMW_part_set_log_n_neighbor_pe(int domain, int _n_neighbor_pe) {
  if (is_init == 0) {
    HECMW_set_error(HECMW_PART_E_LOG_INIT_NOT_YET, "");
    return -1;
  }

  if (domain < 0) {
    HECMW_set_error(HECMW_PART_E_DOMAIN_MIN, "domain");
    return -1;
  }
  if (domain >= n_domain) {
    HECMW_set_error(HECMW_PART_E_DOMAIN_MAX, "domain");
    return -1;
  }

  if (_n_neighbor_pe < 0) {
    HECMW_set_error(HECMW_PART_E_NNEIGHBORPE_LOWER, "_n_neighbor_pe");
    return -1;
  }

  n_neighbor_pe[domain] = _n_neighbor_pe;

  return 0;
}

/*================================================================================================*/

extern int HECMW_part_print_log(void) {
  int i;
  FILE *fp;

  if (is_init == 0) {
    HECMW_set_error(HECMW_PART_E_LOG_INIT_NOT_YET, "");
    return -1;
  }

  if ((fp = fopen(HECMW_PART_LOG_NAME, "w")) == NULL) {
    HECMW_set_error(errno, "log file for partitioner");
    return -1;
  }

  fprintf(fp, "# log file for partitioner ( %s )\n", HECMW_get_date());
  fprintf(fp, "\n");

  /* conditions */
  fprintf(fp, "# conditions\n");
  fprintf(fp, "number of sub-domains : %d\n", n_domain);
  fprintf(fp, "partitioning type     : %s\n", part_type);
  fprintf(fp, "partitioning method   : %s\n", part_method);
  fprintf(fp, "depth of overlapping  : %d\n", depth);
  fprintf(fp, "contact partitioning  : %s\n", part_contact);
  if (n_domain == 1) {
    fprintf(fp, "number of edgecut     : ----- / -----\n");
  } else {
    fprintf(fp, "number of edgecut     : %d / %lld\n", n_edgecut, n_edge);
  }

  fprintf(fp, "\n");

  /* information of entire mesh */
  fprintf(fp, "# information of entire mesh\n");
  fprintf(fp, "number of nodes       : %d\n", n_node_g);
  fprintf(fp, "number of elements    : %d\n", n_elem_g);

  fprintf(fp, "\n");

  /* information of distributed mesh */
  fprintf(fp, "# information of distributed mesh\n");
  fprintf(fp, "domain,       nodes,  int. nodes,       elems,  int. elems, neighbor PE\n");

  for (i = 0; i < n_domain; i++) {
    fprintf(fp, "%6d %12d %12d %12d %12d %12d\n", i, n_node[i], nn_internal[i],
            n_elem[i], ne_internal[i], n_neighbor_pe[i]);
  }

  /* close file */
  if (fclose(fp)) {
    HECMW_set_error(HECMW_PART_E_FILE_CLOSE, "log file for partitioner");
    return -1;
  }

  return 0;
}

extern void HECMW_part_finalize_log(void) { clean_log(); }
