/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
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
#include "hecmw_dist_free.h"
#include "hecmw_common.h"
#include "hecmw_reorder.h"
#include "hecmw_map_int.h"
#include "hecmw_set_int.h"
#include "hecmw_hash.h"

#define HECMW_FLAG_VERSION 5

static int global_node_ID_max = -1;
static int global_elem_ID_max = -1;

/* temporary data structures */
static struct hecmw_io_header *_head;
static struct hecmw_io_initial *_init;
static struct hecmw_io_amplitude *_amp;
static struct hecmw_map_int *_node;
static struct hecmw_map_int *_elem;
static struct hecmw_io_egrp *_egrp;
static struct hecmw_io_ngrp *_ngrp;
static struct hecmw_io_sgrp *_sgrp;
static struct hecmw_io_mpc *_mpc;
static struct hecmw_io_material *_mat;
static struct hecmw_io_section *_sect;
static struct hecmw_system_param *_system;
static struct hecmw_io_zero *_zero;
static struct hecmw_io_contact *_contact;

static char grid_filename[HECMW_FILENAME_LEN + 1] = "Unknown";

/*----------------------------------------------------------------------------*/

static void do_logging(int loglv, int msgno, const char *fmt, va_list ap) {
  if (loglv == HECMW_LOG_ERROR) {
    HECMW_set_verror(msgno, fmt, ap);
  } else {
    HECMW_print_vmsg(loglv, msgno, fmt, ap);
  }
}

static void set_err(int msgno, const char *fmt, ...) {
  va_list ap;

  va_start(ap, fmt);
  do_logging(HECMW_LOG_ERROR, msgno, fmt, ap);
  va_end(ap);
}

static void set_warn(int msgno, const char *fmt, ...) {
  va_list ap;

  va_start(ap, fmt);
  do_logging(HECMW_LOG_WARN, msgno, fmt, ap);
  va_end(ap);
}

/*----------------------------------------------------------------------------*/

static int get_gid2lid_node(int gid) {
  size_t clocal;
  int ret;
  HECMW_assert(_node);
  ret = HECMW_map_int_key2local(_node, gid, &clocal);
  HECMW_assert(ret == HECMW_SUCCESS);
  return clocal + 1;
}

static int get_gid2lid_elem(int gid) {
  size_t clocal;
  int ret;
  HECMW_assert(_elem);
  ret = HECMW_map_int_key2local(_elem, gid, &clocal);
  HECMW_assert(ret == HECMW_SUCCESS);
  return clocal + 1;
}

/*----------------------------------------------------------------------------*/

static int make_surf_key(int elem_id, int surf_id) {
  /* return elem_id*10 + surf_id - 1; */
  if (surf_id < 4) {
    return elem_id * 3 + surf_id - 1;
  } else {
    return -(elem_id * 3 + surf_id - 4);
  }
}

static void decode_surf_key(int key, int *elem_id, int *surf_id) {
  /* *elem_id = key/10; */
  /* *surf_id = key%10 + 1; */
  if (key > 0) {
    *elem_id = key / 3;
    *surf_id = key % 3 + 1;
  } else {
    *elem_id = (-key) / 3;
    *surf_id = (-key) % 3 + 4;
  }
}

static int clear(void) {
  if (HECMW_io_free_all()) return -1;

  strcpy(grid_filename, "Unknown");

  _head    = NULL;
  _init    = NULL;
  _amp     = NULL;
  _node    = NULL;
  _elem    = NULL;
  _egrp    = NULL;
  _ngrp    = NULL;
  _sgrp    = NULL;
  _mpc     = NULL;
  _mat     = NULL;
  _sect    = NULL;
  _system  = NULL;
  _zero    = NULL;
  _contact = NULL;

  return 0;
}

/*------------------------------------------------------------------------------
  print functions
*/

static void print_header(FILE *fp) {
  HECMW_assert(fp);

  fprintf(fp, "HEADER:\n");
  fprintf(fp, "%s\n", _head ? _head->header : "none");
  fprintf(fp, "END of HEADER\n");
}

static void print_amp(FILE *fp) {
  struct hecmw_io_amplitude *p;

  HECMW_assert(fp);

  fprintf(fp, "AMPLITUDE:\n");
  for (p = _amp; p; p = p->next) {
    struct hecmw_io_amplitude_item *item;
    fprintf(fp, "NAME: %s, DEFINITION: %d, TIME: %d, VALUE: %d\n", p->name,
            p->type_def, p->type_time, p->type_val);
    for (item = p->item; item; item = item->next) {
      fprintf(fp, "VAL: %E, T: %E\n", item->val, item->table);
    }
  }
  fprintf(fp, "END of AMPLITUDE\n");
}

static void print_init(FILE *fp) {
  struct hecmw_io_initial *p;

  HECMW_assert(fp);

  fprintf(fp, "INITIAL CONDITION:\n");
  for (p = _init; p; p = p->next) {
    fprintf(fp, "TYPE: %d, NODE: %d, NGRP: %s, VAL: %E\n", p->type, p->node,
            (*p->ngrp != '\0') ? p->ngrp : "none", p->val);
  }
  fprintf(fp, "END of INITIAL CONDITION\n");
}

static void print_node(FILE *fp) {
  int seq, id;
  struct hecmw_io_node *p;

  HECMW_assert(fp);
  HECMW_assert(_node);

  fprintf(fp, "NODE:\n");

  seq = 1;
  HECMW_map_int_iter_init(_node);
  while (HECMW_map_int_iter_next(_node, &id, (void **)&p))
    fprintf(fp, "Node %d: ID=%d: %E  %E  %E\n", seq++, id, p->x, p->y, p->z);

  fprintf(fp, "END of NODE\n");
}

static void print_elem(FILE *fp) {
  int seq, n, id, j;
  struct hecmw_io_element *p;

  HECMW_assert(fp);
  HECMW_assert(_elem);

  fprintf(fp, "ELEMENT:\n");

  seq = 1;
  HECMW_map_int_iter_init(_elem);
  while (HECMW_map_int_iter_next(_elem, &id, (void **)&p)) {
    fprintf(fp, "Element %d: ID=%d: TYPE=%d: ", seq++, id, p->type);
    n = HECMW_get_max_node(p->type);
    for (j = 0; j < n; j++) {
      fprintf(fp, "%d ", p->node[j]);
    }
    fprintf(fp, ": MATITEM: ");
    if (p->nmatitem == 0) {
      fprintf(fp, "none");
    } else {
      for (j = 0; j < p->nmatitem; j++) {
        fprintf(fp, "%E ", p->matitem[j]);
      }
    }
    fprintf(fp, "\n");
  }

  fprintf(fp, "END of ELEMENT\n");
}

static void print_ngrp(FILE *fp) {
  int i;
  const int NITEM = 10;
  struct hecmw_io_ngrp *p;

  HECMW_assert(fp);

  fprintf(fp, "NGROUP:\n");
  for (p = _ngrp; p; p = p->next) {
    int id;
    fprintf(fp, "NAME=%s:\n", p->name);
    HECMW_set_int_iter_init(p->node);
    for (i = 0; HECMW_set_int_iter_next(p->node, &id); i++) {
      fprintf(fp, "%d %c", id, (i + 1) % NITEM ? ' ' : '\n');
    }
    if (i % NITEM) {
      fprintf(fp, "\n");
    }
  }
  fprintf(fp, "END of NGROUP\n");
}

static void print_egrp(FILE *fp) {
  int i;
  const int NITEM = 10;
  struct hecmw_io_egrp *p;

  HECMW_assert(fp);

  fprintf(fp, "EGROUP:\n");
  for (p = _egrp; p; p = p->next) {
    int id;
    fprintf(fp, "NAME=%s:\n", p->name);
    HECMW_set_int_iter_init(p->elem);
    for (i = 0; HECMW_set_int_iter_next(p->elem, &id); i++) {
      fprintf(fp, "%d %c", id, (i + 1) % NITEM ? ' ' : '\n');
    }
    if (i % NITEM) {
      fprintf(fp, "\n");
    }
  }
  fprintf(fp, "END of EGROUP\n");
}

static void print_sgrp(FILE *fp) {
  int i;
  const int NITEM = 10;
  struct hecmw_io_sgrp *p;

  HECMW_assert(fp);

  fprintf(fp, "SGROUP:\n");
  for (p = _sgrp; p; p = p->next) {
    int id, eid, sid;
    fprintf(fp, "NAME=%s:\n", p->name);
    HECMW_set_int_iter_init(p->item);
    for (i = 0; HECMW_set_int_iter_next(p->item, &id); i++) {
      decode_surf_key(id, &eid, &sid);
      fprintf(fp, "%d %d %c", eid, sid, (i + 1) % NITEM ? ' ' : '\n');
    }
    if (i % NITEM) {
      fprintf(fp, "\n");
    }
  }
  fprintf(fp, "END of SGROUP\n");
}

static void print_sect(FILE *fp) {
  struct hecmw_io_section *p;

  HECMW_assert(fp);

  fprintf(fp, "SECTION:\n");
  for (p = _sect; p; p = p->next) {
    fprintf(fp, "EGRP: %s, MATERIAL: %s, COMPOSITE: %d, SECOPT: %d\n", p->egrp,
            p->material, p->composite, p->secopt);
    if (p->type == HECMW_SECT_TYPE_SOLID) {
      fprintf(fp, "TYPE: SOLID, THICKNESS: %E\n", p->sect.solid.thickness);
    } else if (p->type == HECMW_SECT_TYPE_SHELL) {
      fprintf(fp, "TYPE: SHELL, THICKNESS: %E, INTEGPOINTS: %d\n",
              p->sect.shell.thickness, p->sect.shell.integpoints);
    } else if (p->type == HECMW_SECT_TYPE_BEAM) {
      fprintf(fp, "TYPE: BEAM, Reference vector: %E %E %E, Iyy: %E\n",
              p->sect.beam.vxyz[0], p->sect.beam.vxyz[1], p->sect.beam.vxyz[2],
              p->sect.beam.Iyy);
    } else if (p->type == HECMW_SECT_TYPE_INTERFACE) {
      fprintf(fp,
              "TYPE: INTERFACE, THICKNESS: %E, "
              "GAPCON: %E, GAPRAD1: %E, GAPRAD2: %E\n",
              p->sect.interface.thickness, p->sect.interface.gapcon,
              p->sect.interface.gaprad1, p->sect.interface.gaprad2);
    }
  }
  fprintf(fp, "END of SECTION\n");
}

static void print_mat(FILE *fp) {
  int i, j;
  struct hecmw_io_material *p;

  HECMW_assert(fp);

  fprintf(fp, "MATERIAL:\n");
  for (p = _mat; p; p = p->next) {
    fprintf(fp, "NAME: %s\n", p->name);
    for (i = 0; i < p->nitem; i++) {
      struct hecmw_io_matitem *item = &p->item[i];
      struct hecmw_io_matsubitem *p;
      fprintf(fp, "ITEM=%d, SUBITEM=%d:\n", item->item, item->nval);
      for (p = item->subitem; p; p = p->next) {
        fprintf(fp, "VAL: ");
        for (j = 0; j < item->nval; j++) {
          fprintf(fp, "%E ", p->val[j]);
        }
        fprintf(fp, "TEMP: %E\n", p->temp);
      }
    }
  }
  fprintf(fp, "END of MATERIAL\n");
}

static void print_mpc(FILE *fp) {
  int i;
  struct hecmw_io_mpc *p;

  HECMW_assert(fp);

  fprintf(fp, "EQUATION:\n");
  for (p = _mpc; p; p = p->next) {
    fprintf(fp, "NEQ: %d\n", p->neq);
    for (i = 0; i < p->neq; i++) {
      struct hecmw_io_mpcitem *item = &p->item[i];
      fprintf(fp, "ngrp: %s, nod: %d, DOF: %d, A: %E\n",
              (item->node == -1) ? item->ngrp : "(none)", item->node, item->dof,
              item->a);
    }
  }
  fprintf(fp, "END of EQUATION\n");
}

static void print_system(FILE *fp) {
  struct hecmw_system_param param;

  HECMW_assert(fp);

  if (_system) {
    param = *_system;
  } else {
    memset(&param, 0, sizeof(param));
  }

  fprintf(fp, "SYSTEM:\n");
  fprintf(fp, "%E %E %E\n", param.xa, param.ya, param.za);
  fprintf(fp, "%E %E %E\n", param.xb, param.yb, param.zb);
  fprintf(fp, "%E %E %E\n", param.xc, param.yc, param.zc);
  fprintf(fp, "END of SYSTEM\n");
}

static void print_zero(FILE *fp) {
  HECMW_assert(fp);

  fprintf(fp, "ZERO:\n");
  fprintf(fp, "%E\n", _zero ? _zero->zero : 0.0);
  fprintf(fp, "END of ZERO\n");
}

static void print_contact(FILE *fp) {
  int i;
  struct hecmw_io_contact *p;

  HECMW_assert(fp);

  fprintf(fp, "CONTACT PAIR:\n");
  for (p = _contact; p; p = p->next) {
    fprintf(fp, "NAME=%s, ", p->name);
    if (p->type == HECMW_CONTACT_TYPE_NODE_SURF) {
      fprintf(fp, "TYPE=NODE-SURF, ");
    } else if (p->type == HECMW_CONTACT_TYPE_SURF_SURF) {
      fprintf(fp, "TYPE=SURF-SURF, ");
    } else if (p->type == HECMW_CONTACT_TYPE_NODE_ELEM) {
      fprintf(fp, "TYPE=NODE-ELEM, ");
    }
    fprintf(fp, "SLAVE_GRP=%s, MASTER_GRP=%s\n", p->slave_grp, p->master_grp);
  }
  fprintf(fp, "END of CONTACT PAIR\n");
}

void HECMW_io_print_all(FILE *fp) {
  HECMW_assert(fp);

  print_header(fp);
  fprintf(fp, "\n");
  print_zero(fp);
  fprintf(fp, "\n");
  print_init(fp);
  fprintf(fp, "\n");
  print_amp(fp);
  fprintf(fp, "\n");
  print_system(fp);
  fprintf(fp, "\n");
  print_node(fp);
  fprintf(fp, "\n");
  print_elem(fp);
  fprintf(fp, "\n");
  print_ngrp(fp);
  fprintf(fp, "\n");
  print_egrp(fp);
  fprintf(fp, "\n");
  print_sgrp(fp);
  fprintf(fp, "\n");
  print_sect(fp);
  fprintf(fp, "\n");
  print_mat(fp);
  fprintf(fp, "\n");
  print_mpc(fp);
  fprintf(fp, "\n");
  print_contact(fp);
  fprintf(fp, "\n");
}

/*------------------------------------------------------------------------------
  free
*/

static int free_header(struct hecmw_io_header *header) {
  if (header == NULL) return 0;

  HECMW_free(header);
  return 0;
}

static int free_zero(struct hecmw_io_zero *zero) {
  if (zero == NULL) return 0;

  HECMW_free(zero);
  return 0;
}

static int free_node(struct hecmw_map_int *node) {
  if (node == NULL) return 0;

  HECMW_map_int_finalize(node);
  HECMW_free(node);
  return 0;
}

static int free_elem(struct hecmw_map_int *elem) {
  if (elem == NULL) return 0;

  HECMW_map_int_finalize(elem);
  HECMW_free(elem);
  return 0;
}

static int free_ngrp(struct hecmw_io_ngrp *ngrp) {
  struct hecmw_io_ngrp *p, *q;

  for (p = ngrp; p; p = q) {
    q = p->next;
    HECMW_set_int_finalize(p->node);
    HECMW_free(p->node);
    HECMW_free(p);
  }
  return 0;
}

static int free_egrp(struct hecmw_io_egrp *egrp) {
  struct hecmw_io_egrp *p, *q;

  for (p = egrp; p; p = q) {
    q = p->next;
    HECMW_set_int_finalize(p->elem);
    HECMW_free(p->elem);
    HECMW_free(p);
  }
  return 0;
}

static int free_sgrp(struct hecmw_io_sgrp *sgrp) {
  struct hecmw_io_sgrp *p, *q;

  for (p = sgrp; p; p = q) {
    q = p->next;
    HECMW_set_int_finalize(p->item);
    HECMW_free(p->item);
    HECMW_free(p);
  }
  return 0;
}

static int free_mpc(struct hecmw_io_mpc *mpc) {
  struct hecmw_io_mpc *p, *q;

  for (p = mpc; p; p = q) {
    q = p->next;
    HECMW_free(p->item);
    HECMW_free(p);
  }
  return 0;
}

static int free_amp(struct hecmw_io_amplitude *amp) {
  struct hecmw_io_amplitude *p, *q;
  struct hecmw_io_amplitude_item *pp, *qq;

  for (p = amp; p; p = q) {
    q = p->next;
    for (pp = p->item; pp; pp = qq) {
      qq = pp->next;
      HECMW_free(pp);
    }
    HECMW_free(p);
  }
  return 0;
}

static int free_init(struct hecmw_io_initial *init) {
  struct hecmw_io_initial *p, *q;

  for (p = init; p; p = q) {
    q = p->next;
    HECMW_free(p);
  }
  return 0;
}

static int free_material(struct hecmw_io_material *mat) {
  int i;
  struct hecmw_io_material *p, *q;
  struct hecmw_io_matsubitem *pp, *qq;

  for (p = mat; p; p = q) {
    q = p->next;
    for (i = 0; i < p->nitem; i++) {
      for (pp = p->item[i].subitem; pp; pp = qq) {
        qq = pp->next;
        HECMW_free(pp->val);
        HECMW_free(pp);
      }
    }
    HECMW_free(p->item);
    HECMW_free(p);
  }
  return 0;
}

static int free_sect(struct hecmw_io_section *sect) {
  struct hecmw_io_section *p, *q;

  for (p = sect; p; p = q) {
    q = p->next;
    HECMW_free(p);
  }
  return 0;
}

static int free_system(struct hecmw_system_param *system) {
  if (system == NULL) return 0;

  HECMW_free(system);
  return 0;
}

static int free_contact(struct hecmw_io_contact *contact) {
  struct hecmw_io_contact *p, *q;

  for (p = contact; p; p = q) {
    q = p->next;
    HECMW_free(p);
  }
  return 0;
}

int HECMW_io_free_all(void) {
  if (free_header(_head)) return -1;
  if (free_zero(_zero)) return -1;
  if (free_node(_node)) return -1;
  if (free_elem(_elem)) return -1;
  if (free_ngrp(_ngrp)) return -1;
  if (free_egrp(_egrp)) return -1;
  if (free_sgrp(_sgrp)) return -1;
  if (free_mpc(_mpc)) return -1;
  if (free_amp(_amp)) return -1;
  if (free_init(_init)) return -1;
  if (free_material(_mat)) return -1;
  if (free_sect(_sect)) return -1;
  if (free_system(_system)) return -1;
  if (free_contact(_contact)) return -1;
  return 0;
}

/*----------------------------------------------------------------------------*/

int HECMW_io_set_gridfile(char *gridfile) {
  if (gridfile == NULL) gridfile = "";

  strcpy(grid_filename, gridfile);
  return 0;
}

struct hecmw_io_amplitude *HECMW_io_add_amp(const char *name, int definition,
                                            int time, int value, double val,
                                            double t) {
  static struct hecmw_io_amplitude *prev_amp = NULL;
  struct hecmw_io_amplitude *p;
  struct hecmw_io_amplitude_item *item;

  if (name == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_io_add_amp(): name");
    return NULL;
  }
  if (strlen(name) > HECMW_NAME_LEN) {
    set_err(HECMW_ALL_E0101, "HECMW_io_add_amp(): name too long");
    return NULL;
  }

  if (prev_amp != NULL && strcmp(prev_amp->name, name) == 0) {
    p = prev_amp;
  } else {
    p = HECMW_malloc(sizeof(*p));
    if (p == NULL) {
      set_err(errno, "");
      return NULL;
    }
    strcpy(p->name, name);
    p->next = NULL;
    p->item = NULL;
    p->last = NULL;

    if (prev_amp == NULL) {
      _amp = p;
    } else {
      prev_amp->next = p;
    }
    prev_amp = p;
  }
  p->type_def  = definition;
  p->type_time = time;
  p->type_val  = value;

  item = HECMW_malloc(sizeof(*item));
  if (item == NULL) {
    set_err(errno, "");
    return NULL;
  }
  item->next  = NULL;
  item->val   = val;
  item->table = t;

  if (p->last == NULL) {
    p->item = item;
    p->last = item;
  } else {
    p->last->next = item;
    p->last       = item;
  }

  return p;
}

struct hecmw_io_initial *HECMW_io_get_initial(int node) {
  struct hecmw_io_initial *p;

  for (p = _init; p; p = p->next) {
    if (p->node == node) return p;
  }
  return NULL;
}

struct hecmw_io_initial *HECMW_io_add_initial(int type, int node,
                                              const char *ngrp, double val) {
  static struct hecmw_io_initial *prev_init = NULL;
  struct hecmw_io_initial *p;

  if (ngrp == NULL && node <= 0) {
    set_err(HECMW_ALL_E0101, "HECMW_io_add_initial(): ngrp,node");
    return NULL;
  }

  p = HECMW_malloc(sizeof(*p));
  if (p == NULL) {
    set_err(errno, "");
    return NULL;
  }

  if (ngrp) {
    strcpy(p->ngrp, ngrp);
  }
  p->type = type;
  p->node = ngrp ? -1 : node;
  p->val  = val;
  p->next = NULL;

  if (prev_init == NULL) {
    _init = p;
  } else {
    prev_init->next = p;
  }
  prev_init = p;

  return p;
}

struct hecmw_io_element *HECMW_io_get_elem(int id) {
  HECMW_assert(_elem);

  return (struct hecmw_io_element *)HECMW_map_int_get(_elem, id);
}

int HECMW_io_get_n_elem(void) {
  HECMW_assert(_elem);

  return HECMW_map_int_nval(_elem);
}

int HECMW_io_get_elem_max_id(void) {
  int id, max = 0;
  struct hecmw_io_element *val;

  HECMW_assert(_elem);

  HECMW_map_int_iter_init(_elem);
  while (HECMW_map_int_iter_next(_elem, &id, (void **)&val)) {
    if (id > max) {
      max = id;
    }
  }
  return max;
}

static void free_io_elem(void *io_elem) {
  struct hecmw_io_element *p;
  p = (struct hecmw_io_element *)io_elem;
  HECMW_free(p->node);
  HECMW_free(p->matitem);
  HECMW_free(p);
}

struct hecmw_io_element *HECMW_io_add_elem(int id, int type, int *node,
                                           int nmatitem, double *matitem) {
  int nnode;
  int *new_node;
  double *new_matitem;
  struct hecmw_io_element *new_elem;

  if (node == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_io_add_elem(): node");
    return NULL;
  }
  if (nmatitem < 0) {
    set_err(HECMW_ALL_E0101, "HECMW_io_add_elem(): nmatitem");
    return NULL;
  }

  /* get # of connectivity */
  nnode = HECMW_get_max_node(type);
  HECMW_assert(nnode > 0);
  HECMW_assert(nnode <= HECMW_MAX_NODE_MAX);

  new_node = HECMW_malloc(sizeof(*new_node) * nnode);
  if (new_node == NULL) {
    set_err(errno, "");
    return NULL;
  }
  memcpy(new_node, node, sizeof(*new_node) * nnode);

  /* material */
  new_matitem = NULL;
  if (nmatitem > 0) {
    new_matitem = HECMW_malloc(sizeof(*new_matitem) * nmatitem);
    if (new_matitem == NULL) {
      set_err(errno, "");
      return NULL;
    }
    memcpy(new_matitem, matitem, sizeof(*new_matitem) * nmatitem);
  }

  new_elem = HECMW_malloc(sizeof(*new_elem));
  if (new_elem == NULL) {
    set_err(errno, "");
    return NULL;
  }
  new_elem->type       = type;
  new_elem->node       = new_node;
  new_elem->nmatitem   = nmatitem;
  new_elem->matitem    = new_matitem;
  new_elem->mpc_matid  = -1;
  new_elem->mpc_sectid = -1;

  if (_elem == NULL) {
    _elem = (struct hecmw_map_int *)HECMW_malloc(sizeof(struct hecmw_map_int));
    if (_elem == NULL) {
      set_err(errno, "");
      return NULL;
    }
    if (HECMW_map_int_init(_elem, free_io_elem)) {
      set_err(errno, "");
      return NULL;
    }
  }

  if (HECMW_map_int_add(_elem, id, new_elem)) {
    set_err(errno, "");
    return NULL;
  }

  if (id > global_elem_ID_max) {
    global_elem_ID_max = id;
  }

  return new_elem;
}

struct hecmw_io_egrp *HECMW_io_get_egrp(const char *name) {
  extern hecmw_hash_p *hash_eg;

  return (struct hecmw_io_egrp *)hecmw_hash_p_get(hash_eg, name);
}

struct hecmw_io_id_array *HECMW_io_get_elem_in_egrp(const char *name) {
  int nval, i, eid;
  struct hecmw_io_egrp *egrp;
  struct hecmw_io_id_array *id = NULL;

  egrp = HECMW_io_get_egrp(name);
  if (egrp == NULL) goto error;

  nval = HECMW_set_int_nval(_egrp->elem);
  HECMW_assert(nval > 0);

  id = HECMW_malloc(sizeof(*id));
  if (id == NULL) {
    set_err(errno, "");
    goto error;
  }

  id->id = HECMW_malloc(sizeof(*id->id) * nval);
  if (id->id == NULL) {
    set_err(errno, "");
    goto error;
  }

  id->n = nval;

  HECMW_set_int_iter_init(_egrp->elem);
  for (i = 0; HECMW_set_int_iter_next(_egrp->elem, &eid); i++) {
    id->id[i] = eid;
  }

  HECMW_assert(i == nval);

  return id;

error:
  if (id) {
    if (id->id) HECMW_free(id->id);
    HECMW_free(id);
  }
  return NULL;
}

int HECMW_io_add_egrp(const char *name, int nelem, int *elem) {
  int i;
  static struct hecmw_io_egrp *cp = NULL;
  struct hecmw_io_egrp *p;
  extern hecmw_hash_p *hash_eg;

  if (name == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_io_add_egrp(): name");
    return -1;
  }
  if (elem == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_io_add_egrp(): elem");
    return -1;
  }
  if (nelem <= 0) {
    set_err(HECMW_ALL_E0101, "HECMW_io_add_egrp(): nelem");
    return -1;
  }

  /**  for(p=_egrp; p; p=p->next) {
          if(strcmp(p->name, name) == 0) break;
          q = p;
  }**/

  p = (struct hecmw_io_egrp *)hecmw_hash_p_get(hash_eg, name);
  if (p == NULL) {
    p = HECMW_malloc(sizeof(struct hecmw_io_egrp));
    if (p == NULL) {
      set_err(errno, "");
      return -1;
    }
    strcpy(p->name, name);
    p->elem =
        (struct hecmw_set_int *)HECMW_malloc(sizeof(struct hecmw_set_int));
    if (p->elem == NULL) {
      set_err(errno, "");
      return -1;
    }
    if (HECMW_set_int_init(p->elem)) {
      set_err(errno, "");
      return -1;
    }
    p->next = NULL;

    if (cp != NULL) {
      cp->next = p;
    } else {
      _egrp = p;
    }
    cp = p;
  }

  for (i = 0; i < nelem; i++) {
    if (HECMW_set_int_add(p->elem, elem[i])) {
      set_err(errno, "");
      return -1;
    }
  }

  if (HECMW_set_int_is_empty(p->elem)) {
    /* new group && ignored all */
    HECMW_assert(nelem == 0);
    HECMW_set_int_finalize(p->elem);
    HECMW_free(p->elem);
    HECMW_free(p);
    return 0;
  }

  if (hecmw_hash_p_put(hash_eg, name, (void *)p) == 0) {
    printf("HECMW HASH TABLE PUT ERROR\n");
    return -1;
  }

  /* if(cp != NULL) { */
  /* 	cp->next = p; */
  /* } else if (strcmp(cp->name, name) != 0) { */
  /* 	_egrp = p; */
  /* } */
  /* cp = p; */

  return nelem;
}

struct hecmw_io_node *HECMW_io_get_node(int id) {
  HECMW_assert(_node);

  return (struct hecmw_io_node *)HECMW_map_int_get(_node, id);
}

int HECMW_io_get_n_node(void) {
  HECMW_assert(_node);

  return HECMW_map_int_nval(_node);
}

static void free_io_node(void *io_node) {
  struct hecmw_io_node *p;
  p = (struct hecmw_io_node *)io_node;
  /* nothing to do on members */
  HECMW_free(p);
}

struct hecmw_io_node *HECMW_io_add_node(int id, double x, double y, double z) {
  struct hecmw_io_node *new_node;

  new_node = HECMW_malloc(sizeof(*new_node));
  if (new_node == NULL) {
    set_err(errno, "");
    return NULL;
  }
  new_node->x = x;
  new_node->y = y;
  new_node->z = z;

  if (_node == NULL) {
    _node = (struct hecmw_map_int *)HECMW_malloc(sizeof(struct hecmw_map_int));
    if (_node == NULL) {
      set_err(errno, "");
      return NULL;
    }
    if (HECMW_map_int_init(_node, free_io_node)) {
      set_err(errno, "");
      return NULL;
    }
  }

  if (HECMW_map_int_add(_node, id, new_node)) {
    set_err(errno, "");
    return NULL;
  }

  if (id > global_node_ID_max) {
    global_node_ID_max = id;
  }

  return new_node;
}

int HECMW_io_get_nnode_in_ngrp(const char *name) {
  struct hecmw_io_ngrp *p;

  if (name == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_io_get_nnode_in_ngrp(): name");
    return -1;
  }

  for (p = _ngrp; p; p = p->next) {
    if (strcmp(p->name, name) == 0) break;
  }
  if (p == NULL) return 0;

  return HECMW_set_int_nval(p->node);
}

/*
int
HECMW_io_remove_node(int id)
{
        return HECMW_map_int_del(_node, id);
}
*/

struct hecmw_io_ngrp *HECMW_io_get_ngrp(const char *name) {
  /* struct hecmw_io_ngrp *p; */
  extern hecmw_hash_p *hash_ng;

  if (name == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_io_get_ngrp(): name");
    return NULL;
  }

  /* for(p=_ngrp; p; p=p->next) { */
  /* 	if(strcmp(p->name, name) == 0) break; */
  /* } */
  /* return p; */

  return (struct hecmw_io_ngrp *)hecmw_hash_p_get(hash_ng, name);
}

struct hecmw_io_id_array *HECMW_io_get_node_in_ngrp(const char *name) {
  int n, i, nid;
  struct hecmw_io_ngrp *ngrp;
  struct hecmw_io_id_array *id;

  if (name == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_io_get_node_in_ngrp(): name");
    return NULL;
  }

  ngrp = HECMW_io_get_ngrp(name);
  if (ngrp == NULL) return NULL;

  id = HECMW_malloc(sizeof(*id));
  if (id == NULL) {
    set_err(errno, "");
    return NULL;
  }

  n = HECMW_set_int_nval(ngrp->node);
  HECMW_assert(n > 0);

  id->id = HECMW_malloc(sizeof(*id->id) * n);
  if (id->id == NULL) {
    set_err(errno, "");
    return NULL;
  }

  id->n = n;
  HECMW_set_int_iter_init(ngrp->node);
  for (i = 0; HECMW_set_int_iter_next(ngrp->node, &nid); i++) {
    id->id[i] = nid;
  }

  HECMW_assert(i == n);

  return id;
}

int HECMW_io_add_ngrp(const char *name, int nnode, int *node) {
  int i;
  static struct hecmw_io_ngrp *prev_ngrp = NULL;
  struct hecmw_io_ngrp *p;
  extern hecmw_hash_p *hash_ng;

  if (name == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_io_add_ngrp(): name");
    return -1;
  }
  if (node == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_io_add_ngrp(): node");
    return -1;
  }
  if (nnode <= 0) {
    set_err(HECMW_ALL_E0101, "HECMW_io_add_ngrp(): nnode");
    return -1;
  }

  p = (struct hecmw_io_ngrp *)hecmw_hash_p_get(hash_ng, name);

  /* if(prev_ngrp != NULL && strcmp(prev_ngrp->name, name) == 0) { */
  /* 	p = prev_ngrp; */
  /* } else { */
  if (p == NULL) {
    p = HECMW_malloc(sizeof(*p));
    if (p == NULL) {
      set_err(errno, "");
      return -1;
    }
    strcpy(p->name, name);
    p->node =
        (struct hecmw_set_int *)HECMW_malloc(sizeof(struct hecmw_set_int));
    if (p->node == NULL) {
      set_err(errno, "");
      return -1;
    }
    if (HECMW_set_int_init(p->node)) {
      set_err(errno, "");
      return -1;
    }
    p->next = NULL;

    if (prev_ngrp == NULL) {
      _ngrp = p;
    } else {
      prev_ngrp->next = p;
    }
    prev_ngrp = p;
  }

  for (i = 0; i < nnode; i++) {
    if (HECMW_set_int_add(p->node, node[i])) {
      set_err(errno, "");
      return -1;
    }
  }

  if (HECMW_set_int_is_empty(p->node)) {
    /* new group && ignored all */
    HECMW_set_int_finalize(p->node);
    HECMW_free(p->node);
    HECMW_free(p);
    return 0;
  }

  if (hecmw_hash_p_put(hash_ng, name, (void *)p) == 0) {
    printf("HECMW HASH TABLE PUT ERROR\n");
    return -1;
  }

  /* if(prev_ngrp == NULL) { */
  /* 	_ngrp = p; */
  /* } else if(strcmp(prev_ngrp->name, name) != 0) { */
  /* 	prev_ngrp->next = p; */
  /* } */
  /* prev_ngrp = p; */

  return nnode;
}

static int HECMW_io_remove_node_in_ngrp(int node) {
  struct hecmw_io_ngrp *p, *q, *next;

  q = NULL;
  for (p = _ngrp; p; p = next) {
    HECMW_set_int_del(p->node, node);
    if (HECMW_set_int_is_empty(p->node)) {
      /* no node in this group */
      if (q == NULL) {
        _ngrp = p->next;
      } else {
        q->next = p->next;
      }
      next = p->next;
      HECMW_set_int_finalize(p->node);
      HECMW_free(p->node);
      HECMW_free(p);
    } else {
      q    = p;
      next = p->next;
    }
  }
  return 0;
}

static struct hecmw_io_sgrp *HECMW_io_get_sgrp(const char *name) {
  /* struct hecmw_io_sgrp *p;   */
  extern hecmw_hash_p *hash_sg;

  if (name == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_io_get_sgrp(): name");
    return NULL;
  }

  /* for(p=_sgrp; p; p=p->next) { */
  /* 	if(strcmp(p->name, name) == 0) break; */
  /* } */
  /* return p; */

  return (struct hecmw_io_sgrp *)hecmw_hash_p_get(hash_sg, name);
}

int HECMW_io_add_sgrp(const char *name, int n_item, int *elem, int *surf) {
  int i;
  static struct hecmw_io_sgrp *prev_sgrp = NULL;
  struct hecmw_io_sgrp *p;
  extern hecmw_hash_p *hash_sg;

  if (name == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_add_sgrp(): name");
    return -1;
  }
  if (elem == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_add_sgrp(): elem");
    return -1;
  }
  if (surf == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_add_sgrp(): surf");
    return -1;
  }
  if (n_item <= 0) {
    set_err(HECMW_ALL_E0101, "HECMW_add_sgrp(): n_item");
    return -1;
  }

  p = (struct hecmw_io_sgrp *)hecmw_hash_p_get(hash_sg, name);
  if (p == NULL) {
    /* if(prev_sgrp != NULL && strcmp(prev_sgrp->name, name) == 0) { */
    /* 	p = prev_sgrp; */
    /* } else { */
    p = HECMW_malloc(sizeof(*p));
    if (p == NULL) {
      set_err(errno, "");
      return -1;
    }
    strcpy(p->name, name);
    p->item =
        (struct hecmw_set_int *)HECMW_malloc(sizeof(struct hecmw_set_int));
    if (p->item == NULL) {
      set_err(errno, "");
      return -1;
    }
    if (HECMW_set_int_init(p->item)) {
      set_err(errno, "");
      return -1;
    }
    p->next = NULL;

    if (prev_sgrp == NULL) {
      _sgrp = p;
    } else {
      prev_sgrp->next = p;
    }
    prev_sgrp = p;
  }

  for (i = 0; i < n_item; i++) {
    if (HECMW_set_int_add(p->item, make_surf_key(elem[i], surf[i]))) {
      set_err(errno, "");
      return -1;
    }
  }

  if (HECMW_set_int_is_empty(p->item)) {
    /* new group && ignored all */
    HECMW_set_int_finalize(p->item);
    HECMW_free(p->item);
    HECMW_free(p);
    return 0;
  }

  if (hecmw_hash_p_put(hash_sg, name, (void *)p) == 0) {
    printf("HECMW HASH TABLE PUT ERROR\n");
    return -1;
  }

  /* if(prev_sgrp == NULL) { */
  /* 	_sgrp = p; */
  /* } else if(strcmp(prev_sgrp->name, name) != 0) { */
  /* 	prev_sgrp->next = p; */
  /* } */
  /* prev_sgrp = p; */

  return n_item;
}

struct hecmw_io_mpc *HECMW_io_add_mpc(int neq,
                                      const struct hecmw_io_mpcitem *mpcitem,
                                      double cnst) {
  int i;
  static struct hecmw_io_mpc *prev_mpc = NULL;
  struct hecmw_io_mpc *p;
  struct hecmw_io_mpcitem *item;

  if (neq <= 0) {
    set_err(HECMW_ALL_E0101, "HECMW_add_mpc(): neq");
    return NULL;
  }
  if (mpcitem == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_add_mpc(): mpcitem");
    return NULL;
  }

  p = HECMW_malloc(sizeof(*p));
  if (p == NULL) {
    set_err(errno, "");
    return NULL;
  }

  item = HECMW_malloc(sizeof(*item) * neq);
  if (item == NULL) {
    set_err(errno, "");
    return NULL;
  }

  for (i = 0; i < neq; i++) {
    const struct hecmw_io_mpcitem *src = &mpcitem[i];
    struct hecmw_io_mpcitem *dst       = &item[i];
    HECMW_assert((src->node == -1) ? (strlen(src->ngrp) > 0) : 1);
    HECMW_assert((src->node != -1) ? (src->node > 0) : 1);
    HECMW_assert(!HECMW_io_check_mpc_dof(src->dof));
    strcpy(dst->ngrp, src->ngrp);
    dst->node = src->node;
    dst->dof  = src->dof;
    dst->a    = src->a;
  }

  p->neq  = neq;
  p->cnst = cnst;
  p->item = item;
  p->next = NULL;

  if (prev_mpc == NULL) {
    _mpc = p;
  } else {
    prev_mpc->next = p;
  }
  prev_mpc = p;

  return p;
}

struct hecmw_io_section *HECMW_io_add_sect(struct hecmw_io_section *sect) {
  static struct hecmw_io_section *prev_sect = NULL;
  struct hecmw_io_section *p;

  if (sect == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_io_add_sect(): sect");
    return NULL;
  }

  p = HECMW_malloc(sizeof(*p));
  if (p == NULL) {
    set_err(errno, "");
    return NULL;
  }

  *p      = *sect;
  p->next = NULL;

  if (prev_sect == NULL) {
    _sect = p;
  } else {
    prev_sect->next = p;
  }
  prev_sect = p;

  return p;
}

struct hecmw_io_material *HECMW_io_get_mat(const char *name) {
  struct hecmw_io_material *p;
  extern hecmw_hash_p *hash_mat;

  if (name == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_io_get_mat(): name");
    return NULL;
  }

  p = (struct hecmw_io_material *)hecmw_hash_p_get(hash_mat, name);

  /* for(p=_mat; p; p=p->next) {
          if(strcmp(p->name, name) == 0) break;
  }*/

  return p;
}

struct hecmw_io_material *HECMW_io_add_mat(const char *name,
                                           struct hecmw_io_material *mat) {
  static struct hecmw_io_material *prev_mat = NULL;
  struct hecmw_io_material *p;
  extern hecmw_hash_p *hash_mat;

  if (mat == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_io_add_mat(): mat");
    return NULL;
  }

  p = (struct hecmw_io_material *)hecmw_hash_p_get(hash_mat, name);
  if (p == NULL) {
    if (hecmw_hash_p_put(hash_mat, name, (void *)mat) == 0) {
      printf("HECMW HASH TABLE PUT ERROR\n");
      return NULL;
    } else {
      if (prev_mat == NULL) {
        _mat = mat;
      } else {
        prev_mat->next = mat;
      }
      prev_mat = mat;
    }
  }

  return mat;
}

void HECMW_io_set_header(struct hecmw_io_header *header) {
  if (header == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_io_set_header(): header");
    return;
  }

  if (_head) {
    HECMW_free(_head);
    set_warn(HECMW_IO_W1010, "");
  }
  _head = header;
}

struct hecmw_system_param *HECMW_io_get_system(void) {
  return _system;
}

void HECMW_io_set_system(struct hecmw_system_param *system) {
  HECMW_free(_system);
  _system = system; /* allow NULL */
}

void HECMW_io_set_zero(struct hecmw_io_zero *zero) {
  if (_zero) {
    HECMW_free(_zero);
    set_warn(HECMW_IO_W1011, "");
  }
  _zero = zero;
}

struct hecmw_io_contact *HECMW_io_add_contact(const char *name, int type,
                                              const char *slave_grp,
                                              const char *master_grp) {
  static struct hecmw_io_contact *prev_contact = NULL;
  struct hecmw_io_contact *p;

  if (slave_grp == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_io_add_contact(): slave_grp");
    return NULL;
  }
  if (master_grp == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_io_add_contact(): master_grp");
    return NULL;
  }

  p = (struct hecmw_io_contact *)HECMW_malloc(sizeof(*p));
  if (p == NULL) {
    set_err(HECMW_ALL_E0101, "HECMW_io_add_contact(): contact");
    return NULL;
  }

  strcpy(p->name, name);
  p->type = type;
  strcpy(p->slave_grp, slave_grp);
  strcpy(p->slave_orisgrp, slave_grp);
  strcpy(p->master_grp, master_grp);
  p->next = NULL;

  if (prev_contact == NULL) {
    _contact = p;
  } else {
    prev_contact->next = p;
  }
  prev_contact = p;

  return p;
}

/*------------------------------------------------------------------------------
  convert to hecmwST_local_mesh
*/

static int setup_flags(struct hecmwST_local_mesh *mesh) {
  HECMW_assert(mesh);

  mesh->hecmw_flag_adapt       = 0;
  mesh->hecmw_flag_initcon     = 0;
  mesh->hecmw_flag_parttype    = HECMW_FLAG_PARTTYPE_UNKNOWN;
  mesh->hecmw_flag_partdepth   = 1;
  mesh->hecmw_flag_version     = HECMW_FLAG_VERSION;
  mesh->hecmw_flag_partcontact = HECMW_FLAG_PARTCONTACT_UNKNOWN;

  return 0;
}

static int setup_gridfile(struct hecmwST_local_mesh *mesh) {
  HECMW_assert(mesh);

  strcpy(mesh->gridfile, grid_filename);

  return 0;
}

static int setup_files(struct hecmwST_local_mesh *mesh) {
  HECMW_assert(mesh);

  mesh->hecmw_n_file = 0;
  mesh->files        = NULL;

  return 0;
}

static int setup_header(struct hecmwST_local_mesh *mesh) {
  char *p;

  HECMW_assert(mesh);

  p = _head ? _head->header : "";
  strcpy(mesh->header, p);

  return 0;
}

static int setup_zero(struct hecmwST_local_mesh *mesh) {
  HECMW_assert(mesh);

  mesh->zero_temp = 0.0; /* default value */
  if (_zero) {
    mesh->zero_temp = _zero->zero;
  }

  return 0;
}

static int setup_init(struct hecmwST_local_mesh *mesh) {
  int i, n;
  size_t size;
  struct hecmw_io_initial *p;

  HECMW_assert(mesh);
  HECMW_assert(mesh->n_node > 0);

  /* initialize */
  mesh->node_init_val_index = NULL;
  mesh->node_init_val_item  = NULL;

  n = 0;
  for (p = _init; p; p = p->next) {
    n++;
  }
  HECMW_log(HECMW_LOG_DEBUG, "setup_init: n = %d", n);

  if (n == 0) {
    mesh->hecmw_flag_initcon = 0;
    return 0;
  }
  mesh->hecmw_flag_initcon = 1;

  /* allocate index initialized with 0 */
  mesh->node_init_val_index =
      HECMW_calloc(mesh->n_node + 1, sizeof(*mesh->node_init_val_index));
  if (mesh->node_init_val_index == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* set index to 1 where initial condition is specified */
  for (p = _init; p; p = p->next) {
    int lid                        = get_gid2lid_node(p->node);
    mesh->node_init_val_index[lid] = 1;
  }

  /* integrate index */
  for (i = 0; i < mesh->n_node; i++) {
    mesh->node_init_val_index[i + 1] += mesh->node_init_val_index[i];
  }

  /* allocate item */
  size                     = sizeof(*mesh->node_init_val_item) * n;
  mesh->node_init_val_item = HECMW_malloc(size);
  if (mesh->node_init_val_item == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* put values into correct places */
  for (p = _init; p; p = p->next) {
    int lid                         = get_gid2lid_node(p->node);
    int index                       = mesh->node_init_val_index[lid] - 1;
    mesh->node_init_val_item[index] = p->val;
  }

  return 0;
}

static int setup_node(struct hecmwST_local_mesh *mesh) {
  int i, id;
  size_t size;
  struct hecmw_io_node *p;

  HECMW_assert(mesh);

  /* initiallize */
  mesh->n_node             = 0;
  mesh->nn_internal        = 0;
  mesh->node_internal_list = NULL;
  mesh->node_ID            = NULL;
  mesh->global_node_ID     = NULL;
  mesh->node               = NULL;
  mesh->n_dof              = 0;
  mesh->n_dof_grp          = 0;
  mesh->n_dof_tot          = 0;
  mesh->node_dof_index     = NULL;
  mesh->node_dof_item      = NULL;
  mesh->node_val_index     = NULL;
  mesh->node_val_item      = NULL;

  /* n_node */
  mesh->n_node = HECMW_io_get_n_node();
  if (mesh->n_node == 0) {
    return 0;
  }

  /* n_node_gross */
  mesh->n_node_gross = mesh->n_node;

  /* nn_middle */
  mesh->nn_middle = mesh->n_node;

  /* nn_internal */
  mesh->nn_internal = mesh->n_node;

  /* node_internal_list */
  size = sizeof(*mesh->node_internal_list) * mesh->nn_internal;
  mesh->node_internal_list = HECMW_malloc(size);
  if (mesh->node_internal_list == NULL) {
    set_err(errno, "");
    return -1;
  }
  /* node_ID */
  size          = sizeof(*mesh->node_ID) * mesh->n_node * 2;
  mesh->node_ID = HECMW_malloc(size);
  if (mesh->node_ID == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* global_node_ID */
  size                 = sizeof(*mesh->global_node_ID) * mesh->n_node;
  mesh->global_node_ID = HECMW_malloc(size);
  if (mesh->global_node_ID == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* node */
  size       = sizeof(*mesh->node) * mesh->n_node * 3;
  mesh->node = HECMW_malloc(size);
  if (mesh->node == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* set node_internal_list, node_ID, global_node_ID, node */
  HECMW_map_int_iter_init(_node);
  for (i = 0; HECMW_map_int_iter_next(_node, &id, (void **)&p); i++) {
    mesh->node_internal_list[i] = i + 1;
    mesh->node_ID[2 * i]        = i + 1;
    mesh->node_ID[2 * i + 1]    = 0;
    mesh->global_node_ID[i]     = id;
    mesh->node[3 * i]           = p->x;
    mesh->node[3 * i + 1]       = p->y;
    mesh->node[3 * i + 2]       = p->z;
  }

  HECMW_assert(i == mesh->n_node);

  return 0;
}

static int setup_elem(struct hecmwST_local_mesh *mesh) {
  int i, j, n, id;
  size_t size, ncon;
  struct hecmw_io_element *p;

  HECMW_assert(mesh);

  /* initialize */
  mesh->n_elem             = 0;
  mesh->ne_internal        = 0;
  mesh->elem_internal_list = NULL;
  mesh->elem_ID            = NULL;
  mesh->global_elem_ID     = NULL;
  mesh->elem_node_index    = NULL;
  mesh->elem_node_item     = NULL;
  mesh->elem_type          = NULL;
  mesh->n_elem_type        = 0;
  mesh->elem_type_index    = NULL;
  mesh->elem_type_item     = NULL;
  mesh->section_ID         = NULL;
  mesh->n_elem_mat_ID      = 0;
  mesh->elem_mat_ID_index  = NULL;
  mesh->elem_mat_ID_item   = NULL;
  mesh->elem_mat_int_index = NULL;
  mesh->elem_mat_int_val   = NULL;
  mesh->elem_val_index     = NULL;
  mesh->elem_val_item      = NULL;

  /* n_elem */
  mesh->n_elem = HECMW_io_get_n_elem();
  HECMW_assert(mesh->n_elem > 0);

  /* n_elem_gross */
  mesh->n_elem_gross = mesh->n_elem;

  /* ne_internal */
  mesh->ne_internal = mesh->n_elem;

  /* elem_internal_list */
  size = sizeof(*mesh->elem_internal_list) * mesh->ne_internal;
  mesh->elem_internal_list = HECMW_malloc(size);
  if (mesh->elem_internal_list == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* elem_ID */
  size          = sizeof(*mesh->elem_ID) * mesh->n_elem * 2;
  mesh->elem_ID = HECMW_malloc(size);
  if (mesh->elem_ID == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* global_elem_ID */
  size                 = sizeof(*mesh->global_elem_ID) * mesh->n_elem;
  mesh->global_elem_ID = HECMW_malloc(size);
  if (mesh->global_elem_ID == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* elem_type */
  size            = sizeof(*mesh->elem_type) * mesh->n_elem;
  mesh->elem_type = HECMW_malloc(size);
  if (mesh->elem_type == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* elem_node_index */
  size                  = sizeof(*mesh->elem_node_index) * (mesh->n_elem + 1);
  mesh->elem_node_index = HECMW_malloc(size);
  if (mesh->elem_node_index == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* count total # of connectivity and set elem_node_index */
  ncon                     = 0;
  mesh->elem_node_index[0] = 0;

  HECMW_map_int_iter_init(_elem);
  for (i = 0; HECMW_map_int_iter_next(_elem, &id, (void **)&p); i++) {
    n = HECMW_get_max_node(p->type);

    HECMW_assert(n > 0);
    ncon += n;
    HECMW_assert(i + 1 <= mesh->n_elem);
    mesh->elem_node_index[i + 1] = mesh->elem_node_index[i] + n;
    HECMW_assert(mesh->elem_node_index[i + 1] <= ncon);
  }

  HECMW_assert(i == mesh->n_elem);
  HECMW_assert(mesh->elem_node_index[mesh->n_elem] == ncon);

  /* elem_node_item */
  size                 = sizeof(*mesh->elem_node_item) * ncon;
  mesh->elem_node_item = HECMW_malloc(size);
  if (mesh->elem_node_item == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* set elem_ID, global_elem_ID, elem_internal_list, elem_type */
  HECMW_map_int_iter_init(_elem);
  for (i = 0; HECMW_map_int_iter_next(_elem, &id, (void **)&p); i++) {
    int start;

    /* connectivity */
    n     = mesh->elem_node_index[i + 1] - mesh->elem_node_index[i];
    start = mesh->elem_node_index[i];
    for (j = 0; j < n; j++) {
      HECMW_assert(start + j <= mesh->elem_node_index[mesh->n_elem]);
      mesh->elem_node_item[start + j] = get_gid2lid_node(p->node[j]);
    }

    mesh->elem_ID[2 * i]        = i + 1;
    mesh->elem_ID[2 * i + 1]    = 0;
    mesh->global_elem_ID[i]     = id;
    mesh->elem_internal_list[i] = i + 1;
    mesh->elem_type[i]          = p->type;
  }

  HECMW_assert(i == mesh->n_elem);

  return 0;
}

static int setup_ngrp(struct hecmwST_local_mesh *mesh) {
  int i, j, nngrp, nnode;
  size_t size;
  struct hecmwST_node_grp *ngrp;
  struct hecmw_io_ngrp *p;

  HECMW_assert(mesh);

  ngrp = HECMW_malloc(sizeof(*ngrp));
  if (ngrp == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* initialize */
  ngrp->n_grp        = 0;
  ngrp->n_bc         = 0;
  ngrp->grp_name     = NULL;
  ngrp->grp_index    = NULL;
  ngrp->grp_item     = NULL;
  ngrp->bc_grp_ID    = NULL;
  ngrp->bc_grp_type  = NULL;
  ngrp->bc_grp_index = NULL;
  ngrp->bc_grp_dof   = NULL;
  ngrp->bc_grp_val   = NULL;

  nngrp = nnode = 0;
  for (p = _ngrp; p; p = p->next) {
    nnode += HECMW_set_int_nval(p->node);
    nngrp++;
  }
  ngrp->n_grp = nngrp;

  if (ngrp->n_grp <= 0) {
    mesh->node_group = ngrp;
    return 0;
  }

  /* grp_name */
  size           = sizeof(*ngrp->grp_name) * ngrp->n_grp;
  ngrp->grp_name = HECMW_malloc(size);
  if (ngrp->grp_name == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* grp_index */
  size            = sizeof(*ngrp->grp_index) * (ngrp->n_grp + 1);
  ngrp->grp_index = HECMW_malloc(size);
  if (ngrp->grp_index == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* grp_item */
  size           = sizeof(*ngrp->grp_item) * nnode;
  ngrp->grp_item = HECMW_malloc(size);
  if (ngrp->grp_item == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* set */
  ngrp->grp_index[0] = 0;
  for (i = 0, p = _ngrp; p; p = p->next, i++) {
    int start = ngrp->grp_index[i], nid;

    HECMW_set_int_iter_init(p->node);
    for (j = 0; HECMW_set_int_iter_next(p->node, &nid); j++) {
      HECMW_assert(start + j < nnode);

      ngrp->grp_item[start + j] = get_gid2lid_node(nid);

      HECMW_assert(ngrp->grp_item[start + j] > 0);
    }

    ngrp->grp_index[i + 1] = ngrp->grp_index[i] + j;
    ngrp->grp_name[i]      = HECMW_strdup(p->name);
    if (ngrp->grp_name[i] == NULL) {
      set_err(errno, "");
      return -1;
    }
  }

  HECMW_assert(ngrp->grp_index[ngrp->n_grp] == nnode);

  mesh->node_group = ngrp;

  return 0;
}

static int setup_egrp(struct hecmwST_local_mesh *mesh) {
  int i, j, negrp, nelem;
  size_t size;
  struct hecmwST_elem_grp *egrp;
  struct hecmw_io_egrp *p;

  HECMW_assert(mesh);

  egrp = HECMW_malloc(sizeof(*egrp));
  if (egrp == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* initialize */
  egrp->n_grp        = 0;
  egrp->n_bc         = 0;
  egrp->grp_name     = NULL;
  egrp->grp_index    = NULL;
  egrp->grp_item     = NULL;
  egrp->bc_grp_ID    = NULL;
  egrp->bc_grp_type  = NULL;
  egrp->bc_grp_index = NULL;
  egrp->bc_grp_val   = NULL;

  negrp = nelem = 0;
  for (p = _egrp; p; p = p->next) {
    nelem += HECMW_set_int_nval(p->elem);
    negrp++;
  }
  egrp->n_grp = negrp;

  if (egrp->n_grp <= 0) {
    mesh->elem_group = egrp;
    return 0;
  }

  /* grp_name */
  size           = sizeof(*egrp->grp_name) * egrp->n_grp;
  egrp->grp_name = HECMW_malloc(size);
  if (egrp->grp_name == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* grp_index */
  size            = sizeof(*egrp->grp_index) * (egrp->n_grp + 1);
  egrp->grp_index = HECMW_malloc(size);
  if (egrp->grp_index == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* grp_item */
  size           = sizeof(*egrp->grp_item) * nelem;
  egrp->grp_item = HECMW_malloc(size);
  if (egrp->grp_item == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* set */
  egrp->grp_index[0] = 0;
  for (i = 0, p = _egrp; p; p = p->next, i++) {
    int eid;

    HECMW_set_int_iter_init(p->elem);
    for (j = 0; HECMW_set_int_iter_next(p->elem, &eid); j++) {
      int start = egrp->grp_index[i];

      HECMW_assert(start + j < nelem);

      egrp->grp_item[start + j] = get_gid2lid_elem(eid);

      HECMW_assert(egrp->grp_item[start + j] > 0);
    }

    egrp->grp_index[i + 1] = egrp->grp_index[i] + j;
    egrp->grp_name[i]      = HECMW_strdup(p->name);

    if (egrp->grp_name[i] == NULL) {
      set_err(errno, "");
      return -1;
    }
  }

  HECMW_assert(egrp->grp_index[egrp->n_grp] == nelem);

  mesh->elem_group = egrp;

  return 0;
}

static int setup_sgrp(struct hecmwST_local_mesh *mesh) {
  int i, j, nsgrp, nelem;
  size_t size;
  struct hecmwST_surf_grp *sgrp;
  struct hecmw_io_sgrp *p;

  HECMW_assert(mesh);

  sgrp = HECMW_malloc(sizeof(*sgrp));
  if (sgrp == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* initialize */
  sgrp->n_grp        = 0;
  sgrp->n_bc         = 0;
  sgrp->grp_name     = NULL;
  sgrp->grp_index    = NULL;
  sgrp->grp_item     = NULL;
  sgrp->bc_grp_ID    = NULL;
  sgrp->bc_grp_type  = NULL;
  sgrp->bc_grp_index = NULL;
  sgrp->bc_grp_val   = NULL;

  nsgrp = nelem = 0;
  for (p = _sgrp; p; p = p->next) {
    nelem += HECMW_set_int_nval(p->item);
    nsgrp++;
  }
  sgrp->n_grp = nsgrp;

  if (sgrp->n_grp <= 0) {
    mesh->surf_group = sgrp;
    return 0;
  }

  /* grp_name */
  size           = sizeof(*sgrp->grp_name) * sgrp->n_grp;
  sgrp->grp_name = HECMW_malloc(size);
  if (sgrp->grp_name == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* grp_index */
  size            = sizeof(*sgrp->grp_index) * (sgrp->n_grp + 1);
  sgrp->grp_index = HECMW_malloc(size);
  if (sgrp->grp_index == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* grp_item */
  size           = sizeof(*sgrp->grp_item) * nelem * 2;
  sgrp->grp_item = HECMW_malloc(size);
  if (sgrp->grp_item == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* set */
  sgrp->grp_index[0] = 0;
  for (i = 0, p = _sgrp; p; p = p->next, i++) {
    int start = sgrp->grp_index[i] * 2, id;

    HECMW_set_int_iter_init(p->item);
    for (j = 0; HECMW_set_int_iter_next(p->item, &id); j++) {
      int eid, sid;

      decode_surf_key(id, &eid, &sid);

      sgrp->grp_item[start + j * 2]     = get_gid2lid_elem(eid);
      sgrp->grp_item[start + j * 2 + 1] = sid;
    }

    sgrp->grp_index[i + 1] = sgrp->grp_index[i] + j;
    sgrp->grp_name[i]      = HECMW_strdup(p->name);

    if (sgrp->grp_name[i] == NULL) {
      set_err(errno, "");
      return -1;
    }
  }

  HECMW_assert(sgrp->grp_index[sgrp->n_grp] == nelem);

  mesh->surf_group = sgrp;

  return 0;
}

static int setup_mpc(struct hecmwST_local_mesh *mesh) {
  int i, j, nmpc, nneq, start;
  size_t size;
  struct hecmwST_mpc *mpc;
  struct hecmw_io_mpc *p;

  HECMW_assert(mesh);

  mpc = HECMW_malloc(sizeof(*mpc));
  if (mpc == NULL) {
    set_err(errno, "");
    goto error;
  }

  mpc->n_mpc     = 0;
  mpc->mpc_index = NULL;
  mpc->mpc_item  = NULL;
  mpc->mpc_dof   = NULL;
  mpc->mpc_val   = NULL;
  mpc->mpc_const = NULL;

  if (_mpc == NULL) {
    mesh->mpc = mpc;
    return 0;
  }

  /* count total # of mpc, neq */
  nmpc = 0;
  nneq = 0;
  for (p = _mpc; p; p = p->next) {
    nmpc++;
    nneq += p->neq;
  }
  HECMW_assert(nmpc > 0);
  HECMW_assert(nneq > 0);
  mpc->n_mpc = nmpc;

  /* mpc_index */
  size           = sizeof(*mpc->mpc_index) * (mpc->n_mpc + 1);
  mpc->mpc_index = HECMW_malloc(size);
  if (mpc->mpc_index == NULL) {
    set_err(errno, "");
    goto error;
  }

  /* mpc_item */
  size          = sizeof(*mpc->mpc_item) * nneq;
  mpc->mpc_item = HECMW_malloc(size);
  if (mpc->mpc_item == NULL) {
    set_err(errno, "");
    goto error;
  }

  /* mpc_dof */
  size         = sizeof(*mpc->mpc_dof) * nneq;
  mpc->mpc_dof = HECMW_malloc(size);
  if (mpc->mpc_dof == NULL) {
    set_err(errno, "");
    goto error;
  }

  /* mpc_val */
  size         = sizeof(*mpc->mpc_val) * nneq;
  mpc->mpc_val = HECMW_malloc(size);
  if (mpc->mpc_val == NULL) {
    set_err(errno, "");
    goto error;
  }

  /* mpc_const */
  size           = sizeof(*mpc->mpc_const) * nmpc;
  mpc->mpc_const = HECMW_malloc(size);
  if (mpc->mpc_const == NULL) {
    set_err(errno, "");
    goto error;
  }

  /* set */
  i                 = 0;
  mpc->mpc_index[0] = 0;
  for (p = _mpc; p; p = p->next) {
    HECMW_assert(i + 1 <= mpc->n_mpc);
    mpc->mpc_index[i + 1] = mpc->mpc_index[i] + p->neq;
    HECMW_assert(mpc->mpc_index[i + 1] <= nneq);
    start = mpc->mpc_index[i];
    for (j = 0; j < p->neq; j++) {
      HECMW_assert(start + j < nneq);
      mpc->mpc_item[start + j] = get_gid2lid_node(p->item[j].node);
      HECMW_assert(mpc->mpc_item[start + j]);
      mpc->mpc_dof[start + j] = p->item[j].dof;
      mpc->mpc_val[start + j] = p->item[j].a;
    }
    mpc->mpc_const[i] = p->cnst;
    i++;
  }
  HECMW_assert(i == mpc->n_mpc);
  HECMW_assert(mpc->mpc_index[mpc->n_mpc] == nneq);

  mesh->mpc = mpc;

  return 0;
error:
  if (mpc) {
    HECMW_free(mpc->mpc_index);
    HECMW_free(mpc->mpc_item);
    HECMW_free(mpc->mpc_dof);
    HECMW_free(mpc->mpc_val);
    HECMW_free(mpc->mpc_const);
    HECMW_free(mpc);
  }
  return -1;
}

static int setup_amp(struct hecmwST_local_mesh *mesh) {
  int i, j, namp, nitem, start;
  size_t size;
  struct hecmwST_amplitude *amp;
  struct hecmw_io_amplitude *p;

  HECMW_assert(mesh);

  amp = HECMW_malloc(sizeof(*amp));
  if (amp == NULL) {
    set_err(errno, "");
    return -1;
  }

  amp->n_amp               = 0;
  amp->amp_name            = NULL;
  amp->amp_type_definition = NULL;
  amp->amp_type_time       = NULL;
  amp->amp_type_value      = NULL;
  amp->amp_index           = NULL;
  amp->amp_val             = NULL;
  amp->amp_table           = NULL;

  if (_amp == NULL) {
    mesh->amp = amp;
    return 0;
  }

  /* count total # of amplitude,item */
  namp  = 0;
  nitem = 0;
  for (p = _amp; p; p = p->next) {
    struct hecmw_io_amplitude_item *item;
    for (item = p->item; item; item = item->next) {
      nitem++;
    }
    namp++;
  }
  HECMW_assert(namp > 0);
  HECMW_assert(nitem > 0);
  amp->n_amp = namp;

  /* amp_name */
  size          = sizeof(*amp->amp_name) * amp->n_amp;
  amp->amp_name = HECMW_malloc(size);
  if (amp->amp_name == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* amp_type_definition */
  size                     = sizeof(*amp->amp_type_definition) * amp->n_amp;
  amp->amp_type_definition = HECMW_malloc(size);
  if (amp->amp_type_definition == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* amp_type_time */
  size               = sizeof(*amp->amp_type_time) * amp->n_amp;
  amp->amp_type_time = HECMW_malloc(size);
  if (amp->amp_type_time == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* amp_type_val */
  size                = sizeof(*amp->amp_type_value) * amp->n_amp;
  amp->amp_type_value = HECMW_malloc(size);
  if (amp->amp_type_value == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* amp_index */
  size           = sizeof(*amp->amp_index) * (amp->n_amp + 1);
  amp->amp_index = HECMW_malloc(size);
  if (amp->amp_index == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* amp_val */
  size         = sizeof(*amp->amp_val) * nitem;
  amp->amp_val = HECMW_malloc(size);
  if (amp->amp_val == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* amp_table */
  size           = sizeof(*amp->amp_table) * nitem;
  amp->amp_table = HECMW_malloc(size);
  if (amp->amp_table == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* set */
  i                 = 0;
  amp->amp_index[0] = 0;
  for (p = _amp; p; p = p->next) {
    struct hecmw_io_amplitude_item *item;
    int n = 0;
    for (item = p->item; item; item = item->next) {
      n++;
    }
    HECMW_assert(i + 1 <= namp);
    amp->amp_index[i + 1] = amp->amp_index[i] + n;
    HECMW_assert(amp->amp_index[i + 1] <= nitem);
    start = amp->amp_index[i];
    for (item = p->item, j = 0; item; item = item->next, j++) {
      amp->amp_val[start + j]   = item->val;
      amp->amp_table[start + j] = item->table;
    }
    HECMW_assert(strlen(p->name) < HECMW_NAME_LEN);
    amp->amp_name[i] = HECMW_strdup(p->name);
    if (amp->amp_name[i] == NULL) {
      set_err(errno, "");
      return -1;
    }
    amp->amp_type_definition[i] = p->type_def;
    amp->amp_type_time[i]       = p->type_time;
    amp->amp_type_value[i]      = p->type_val;
    i++;
  }
  HECMW_assert(i == amp->n_amp);
  HECMW_assert(amp->amp_index[amp->n_amp] == nitem);

  mesh->amp = amp;

  return 0;
}

static int setup_adapt(struct hecmwST_local_mesh *mesh) {
  HECMW_assert(mesh);

  /* clear */
  mesh->coarse_grid_level       = 0;
  mesh->n_adapt                 = 0;
  mesh->when_i_was_refined_node = NULL;
  mesh->when_i_was_refined_elem = NULL;
  mesh->adapt_parent_type       = NULL;
  mesh->adapt_type              = NULL;
  mesh->adapt_level             = NULL;
  mesh->adapt_parent            = NULL;
  mesh->adapt_children_index    = NULL;
  mesh->adapt_children_item     = NULL;

  return 0;
}

static int setup_refine(struct hecmwST_local_mesh *mesh) {
  HECMW_assert(mesh);

  /* clear */
  mesh->n_refine           = 0;
  mesh->node_old2new       = NULL;
  mesh->node_new2old       = NULL;
  mesh->elem_old2new       = NULL;
  mesh->elem_new2old       = NULL;
  mesh->n_node_refine_hist = NULL;

  return 0;
}

static int setup_pe(struct hecmwST_local_mesh *mesh) {
  HECMW_assert(mesh);

  mesh->PETOT         = HECMW_comm_get_size();
  mesh->my_rank       = HECMW_comm_get_rank();
  mesh->HECMW_COMM    = HECMW_comm_get_comm();
  mesh->PEsmpTOT      = 1;
  mesh->n_subdomain   = 1;
  mesh->errnof        = 0;
  mesh->n_neighbor_pe = 0;
  mesh->neighbor_pe   = NULL;
  mesh->import_index  = NULL;
  mesh->import_item   = NULL;
  mesh->export_index  = NULL;
  mesh->export_item   = NULL;
  mesh->shared_index  = NULL;
  mesh->shared_item   = NULL;

  if (mesh->my_rank == 0) {
    mesh->zero = 1;
  } else {
    mesh->zero = 0;
  }

  return 0;
}

static int setup_elem_check_sectid(struct hecmwST_local_mesh *mesh) {
  int i;

  HECMW_assert(mesh);

  for (i = 0; i < mesh->n_elem; i++) {
    if (mesh->section_ID[i] == -1) {
      set_err(HECMW_IO_E1012, "Element %d", mesh->global_elem_ID[i]);
      return -1;
    }
  }
  return 0;
}

static int setup_sect_set_sectid(struct hecmwST_local_mesh *mesh,
                                 const char *egrp_name, int sectid) {
  int i, eid, egid, start, end;
  struct hecmwST_elem_grp *egrp;

  HECMW_assert(mesh);
  HECMW_assert(mesh->elem_group);
  HECMW_assert(egrp_name);

  egid = HECMW_dist_get_egrp_id(mesh->elem_group, egrp_name);
  HECMW_assert(egid > 0);

  egrp = mesh->elem_group;

  start = egrp->grp_index[egid - 1];
  end   = egrp->grp_index[egid] - 1;

  HECMW_assert(mesh->section_ID);
  for (i = start; i <= end; i++) {
    eid = egrp->grp_item[i];
    if (mesh->section_ID[eid - 1] != -1) {
      set_err(HECMW_IO_E1012, "Element %d has already had section %d",
              mesh->global_elem_ID[eid - 1], mesh->section_ID[eid - 1]);
      return -1;
    }
    mesh->section_ID[eid - 1] = sectid;
  }
  return 0;
}

static int setup_sect(struct hecmwST_local_mesh *mesh) {
  int i, nsect, nint, nreal, nmat;
  size_t size;
  struct hecmwST_section *sect;
  struct hecmw_io_section *p;

  HECMW_assert(mesh);

  /* mesh->section */
  size = sizeof(*sect);
  sect = HECMW_malloc(size);
  if (sect == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* section_ID */
  size             = sizeof(*mesh->section_ID) * mesh->n_elem;
  mesh->section_ID = HECMW_malloc(size);
  if (mesh->section_ID == NULL) {
    set_err(errno, "");
    return -1;
  }
  memset(mesh->section_ID, -1, size); /* initialize */

  nsect = nint = nreal = nmat = 0;
  for (p = _sect; p; p = p->next) {
    nsect++;
    if (p->type == HECMW_SECT_TYPE_SOLID) {
      nreal++; /* thickness */
    } else if (p->type == HECMW_SECT_TYPE_SHELL) {
      nreal++; /* thickness */
      nint++;  /* integpoints */
    } else if (p->type == HECMW_SECT_TYPE_BEAM) {
      nreal += 7; /* vxyz3, Iyy, Izz, Jx */
    } else if (p->type == HECMW_SECT_TYPE_INTERFACE) {
      nreal += 4; /* thickness, gapcon, gaprad1, gaprad2 */
    } else {
      return -1;
    }
    /*
    if(p->composite > 0) {
            nmat += composite;
    } else {
            nmat++;
    }
    */
    nmat++;
  }
  sect->n_sect = nsect;

  sect->sect_type         = NULL;
  sect->sect_opt          = NULL;
  sect->sect_mat_ID_index = NULL;
  sect->sect_mat_ID_item  = NULL;
  sect->sect_I_index      = NULL;
  sect->sect_I_item       = NULL;
  sect->sect_R_index      = NULL;
  sect->sect_R_item       = NULL;

  if (sect->n_sect <= 0) {
    mesh->section = sect;
    return 0;
  }

  /* sect_type */
  size            = sizeof(*sect->sect_type) * sect->n_sect;
  sect->sect_type = HECMW_malloc(size);
  if (sect->sect_type == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* sect_opt */
  size           = sizeof(*sect->sect_opt) * sect->n_sect;
  sect->sect_opt = HECMW_malloc(size);
  if (sect->sect_opt == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* sect_mat_ID_index */
  size = sizeof(*sect->sect_mat_ID_index) * (sect->n_sect + 1);
  sect->sect_mat_ID_index = HECMW_malloc(size);
  if (sect->sect_mat_ID_index == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* sect_mat_ID_item */
  HECMW_assert(nmat > 0);
  size                   = sizeof(*sect->sect_mat_ID_item) * nmat;
  sect->sect_mat_ID_item = HECMW_malloc(size);
  if (sect->sect_mat_ID_item == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* sect_I_index */
  size               = sizeof(*sect->sect_I_index) * (sect->n_sect + 1);
  sect->sect_I_index = HECMW_malloc(size);
  if (sect->sect_I_index == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* sect_I_item */
  sect->sect_I_item = NULL;
  if (nint > 0) {
    size              = sizeof(*sect->sect_I_item) * nint;
    sect->sect_I_item = HECMW_malloc(size);
    if (sect->sect_I_item == NULL) {
      set_err(errno, "");
      return -1;
    }
  }

  /* sect_R_index */
  size               = sizeof(*sect->sect_R_index) * (sect->n_sect + 1);
  sect->sect_R_index = HECMW_malloc(size);
  if (sect->sect_R_index == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* sect_R_item */
  sect->sect_R_item = NULL;
  if (nreal > 0) {
    size              = sizeof(*sect->sect_R_item) * nreal;
    sect->sect_R_item = HECMW_malloc(size);
    if (sect->sect_R_item == NULL) {
      set_err(errno, "");
      return -1;
    }
  }

  /* set */
  sect->sect_I_index[0]      = 0;
  sect->sect_R_index[0]      = 0;
  sect->sect_mat_ID_index[0] = 0;
  for (i = 0, p = _sect; p; p = p->next, i++) {
    int iidx = sect->sect_I_index[i];
    int ridx = sect->sect_R_index[i];
    int midx;
    if (p->type == HECMW_SECT_TYPE_SOLID) {
      sect->sect_I_index[i + 1] = sect->sect_I_index[i] + 0;
      sect->sect_R_index[i + 1] = sect->sect_R_index[i] + 1;
      HECMW_assert(ridx <= nreal);
      sect->sect_R_item[ridx] = p->sect.solid.thickness;
    } else if (p->type == HECMW_SECT_TYPE_SHELL) {
      sect->sect_I_index[i + 1] = sect->sect_I_index[i] + 1;
      sect->sect_R_index[i + 1] = sect->sect_R_index[i] + 1;
      HECMW_assert(iidx <= nint);
      HECMW_assert(ridx <= nreal);
      sect->sect_I_item[iidx] = p->sect.shell.integpoints;
      sect->sect_R_item[ridx] = p->sect.shell.thickness;
    } else if (p->type == HECMW_SECT_TYPE_BEAM) {
      sect->sect_I_index[i + 1] = sect->sect_I_index[i] + 0;
      sect->sect_R_index[i + 1] = sect->sect_R_index[i] + 7;
      HECMW_assert(ridx + 6 <= nreal);
      sect->sect_R_item[ridx]     = p->sect.beam.vxyz[0];
      sect->sect_R_item[ridx + 1] = p->sect.beam.vxyz[1];
      sect->sect_R_item[ridx + 2] = p->sect.beam.vxyz[2];
      sect->sect_R_item[ridx + 3] = p->sect.beam.area;
      sect->sect_R_item[ridx + 4] = p->sect.beam.Iyy;
      sect->sect_R_item[ridx + 5] = p->sect.beam.Izz;
      sect->sect_R_item[ridx + 6] = p->sect.beam.Jx;
    } else if (p->type == HECMW_SECT_TYPE_INTERFACE) {
      sect->sect_I_index[i + 1] = sect->sect_I_index[i] + 0;
      sect->sect_R_index[i + 1] = sect->sect_R_index[i] + 4;
      HECMW_assert(ridx + 3 <= nreal);
      sect->sect_R_item[ridx]     = p->sect.interface.thickness;
      sect->sect_R_item[ridx + 1] = p->sect.interface.gapcon;
      sect->sect_R_item[ridx + 2] = p->sect.interface.gaprad1;
      sect->sect_R_item[ridx + 3] = p->sect.interface.gaprad2;
    } else {
      return -1;
    }
    sect->sect_type[i]             = p->type;
    sect->sect_opt[i]              = p->secopt;
    sect->sect_mat_ID_index[i + 1] = sect->sect_mat_ID_index[i] + 1;
    midx                           = sect->sect_mat_ID_index[i];
    /* must be called setup_mat() before */
    HECMW_assert(mesh->material);
    sect->sect_mat_ID_item[midx] =
        HECMW_dist_get_mat_id(mesh->material, p->material);
    HECMW_assert(sect->sect_mat_ID_item[midx] > 0);

    /* set mesh->section_id */
    /* depends on setup_egrp() */
    if (setup_sect_set_sectid(mesh, p->egrp, i + 1)) return -1;
  }

  mesh->section = sect;

  return 0;
}

static int setup_mpc_sectid(struct hecmwST_local_mesh *mesh) {
  int i;
  struct hecmw_io_element *elem;

  HECMW_assert(mesh);

  for (i = 0; i < mesh->n_elem; i++) {
    /* depends on setup_elem() */
    if (mesh->elem_type[i] < 900) continue;
    if (mesh->elem_type[i] >= 1000) continue;
    elem = HECMW_io_get_elem(mesh->global_elem_ID[i]);
    HECMW_assert(elem);
    HECMW_assert(elem->mpc_sectid != -1);
    mesh->section_ID[i] = elem->mpc_sectid;
  }
  return 0;
}

static int setup_mpc_reorder(struct hecmwST_local_mesh *mesh) {
  if (HECMW_reorder(mesh)) {
    return -1;
  }
  return 0;
}

static int setup_mat(struct hecmwST_local_mesh *mesh) {
  int i, j, k, l, nmat, nmatitem, nmatsubitem, nmattable;
  size_t size;
  struct hecmwST_material *mat;
  struct hecmw_io_material *p;

  HECMW_assert(mesh);

  /* mesh->material */
  size = sizeof(*mat);
  mat  = HECMW_malloc(size);
  if (mat == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* n_mat, n_mat_item, n_mat_subitem, n_mat_table */
  nmat = nmatitem = nmatsubitem = nmattable = 0;
  for (p = _mat; p; p = p->next) {
    nmat++;
    nmatitem += p->nitem;
    for (i = 0; i < p->nitem; i++) {
      struct hecmw_io_matsubitem *msi = p->item[i].subitem;
      nmatsubitem += p->item[i].nval;
      for (msi = p->item[i].subitem; msi; msi = msi->next) {
        nmattable += p->item[i].nval;
      }
    }
  }
  mat->n_mat         = nmat;
  mat->n_mat_item    = nmatitem;
  mat->n_mat_subitem = nmatsubitem;
  mat->n_mat_table   = nmattable;

  mat->mat_name          = NULL;
  mat->mat_item_index    = NULL;
  mat->mat_subitem_index = NULL;
  mat->mat_table_index   = NULL;
  mat->mat_val           = NULL;
  mat->mat_temp          = NULL;

  if (mat->n_mat <= 0) {
    mesh->material = mat;
    return 0;
  }

  /* mat_name */
  size          = sizeof(*mat->mat_name) * mat->n_mat;
  mat->mat_name = HECMW_malloc(size);
  if (mat->mat_name == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* mat_item_index */
  size                = sizeof(*mat->mat_item_index) * (mat->n_mat + 1);
  mat->mat_item_index = HECMW_malloc(size);
  if (mat->mat_item_index == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* mat_subitem_index */
  size = sizeof(*mat->mat_subitem_index) * (mat->n_mat_item + 1);
  mat->mat_subitem_index = HECMW_malloc(size);
  if (mat->mat_subitem_index == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* mat_table_index */
  size = sizeof(*mat->mat_table_index) * (mat->n_mat_subitem + 1);
  mat->mat_table_index = HECMW_malloc(size);
  if (mat->mat_table_index == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* mat_val */
  size         = sizeof(*mat->mat_val) * nmattable;
  mat->mat_val = HECMW_malloc(size);
  if (mat->mat_val == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* mat_temp */
  size          = sizeof(*mat->mat_temp) * nmattable;
  mat->mat_temp = HECMW_malloc(size);
  if (mat->mat_temp == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* set */
  mat->mat_item_index[0]    = 0;
  mat->mat_subitem_index[0] = 0;
  mat->mat_table_index[0]   = 0;
  for (i = 0, p = _mat; p; p = p->next, i++) {
    HECMW_assert(i + 1 <= mat->n_mat);
    mat->mat_item_index[i + 1] = mat->mat_item_index[i] + p->nitem;
    mat->mat_name[i]           = HECMW_strdup(p->name);
    if (mat->mat_name[i] == NULL) {
      set_err(errno, "");
      return -1;
    }
    for (j = 0; j < p->nitem; j++) {
      int ntable                          = 0;
      struct hecmw_io_matitem *item       = &p->item[j];
      struct hecmw_io_matsubitem *subitem = item->subitem;
      int idx                             = mat->mat_item_index[i] + j;
      HECMW_assert(idx + 1 <= mat->n_mat_item);
      mat->mat_subitem_index[idx + 1] =
          mat->mat_subitem_index[idx] + item->nval;
      for (subitem = item->subitem; subitem; subitem = subitem->next) {
        ntable++;
      }
      for (k = 0; k < item->nval; k++) {
        HECMW_assert(mat->mat_item_index[i] + j <= mat->n_mat_item);
        idx = mat->mat_subitem_index[mat->mat_item_index[i] + j] + k;
        HECMW_assert(idx + 1 <= mat->n_mat_subitem);
        mat->mat_table_index[idx + 1] = mat->mat_table_index[idx] + ntable;
      }
      for (k = 0, subitem = item->subitem; subitem;
           subitem = subitem->next, k++) {
        for (l = 0; l < item->nval; l++) {
          int imat   = mat->mat_item_index[i];
          int ismat  = mat->mat_subitem_index[imat + j];
          int itable = mat->mat_table_index[ismat + l];
          idx        = itable + k;
          HECMW_assert(idx < mat->n_mat_table);
          mat->mat_val[idx]  = subitem->val[l];
          mat->mat_temp[idx] = subitem->temp;
        }
      }
    }
  }

  HECMW_assert(mat->mat_item_index[mat->n_mat] == mat->n_mat_item);
  HECMW_assert(mat->mat_subitem_index[mat->n_mat_item] == mat->n_mat_subitem);
  HECMW_assert(mat->mat_table_index[mat->n_mat_subitem] == mat->n_mat_table);

  mesh->material = mat;

  return 0;
}

static int setup_elem_mat(struct hecmwST_local_mesh *mesh) {
  int i, j, n, id, idx, *start, sectid, nmat, *matid;
  struct mat_table {
    int n;
    int *matid;
  } * mat;

  HECMW_assert(mesh);

  mat = HECMW_malloc(sizeof(*mat) * mesh->n_elem);
  if (mat == NULL) {
    set_err(errno, "");
    return -1;
  }

  for (i = 0; i < mesh->n_elem; i++) {
    mat[i].n = 0;
  }

  nmat = 0;
  for (i = 0; i < mesh->n_elem; i++) {
    struct hecmw_io_element *elem = HECMW_io_get_elem(mesh->global_elem_ID[i]);
    HECMW_assert(elem);
    if (mesh->elem_type[i] >= 900 && mesh->elem_type[i] < 1000) {
      n = 1;
      HECMW_assert(elem->mpc_matid != -1);
      start = &elem->mpc_matid;
    } else if (mesh->elem_type[i] >= 1000 && mesh->elem_type[i] < 1100) {
      n = 1;
      HECMW_assert(elem->mpc_matid != -1);
      start = &elem->mpc_matid;
    } else {
      if (elem->nmatitem > 0) {
        HECMW_assert(mesh->material);
        n  = 1;
        id = HECMW_dist_get_mat_id(mesh->material, elem->matname);
        HECMW_assert(id > 0);
        start = &id;
      } else {
        HECMW_assert(mesh->section);
        sectid = mesh->section_ID[i];
        HECMW_assert(sectid > 0);
        idx = mesh->section->sect_mat_ID_index[sectid - 1];
        n   = mesh->section->sect_mat_ID_index[sectid] - idx;
        HECMW_assert(n > 0);
        start = &mesh->section->sect_mat_ID_item[idx];
      }
    }
    matid = HECMW_malloc(sizeof(matid) * n);
    if (matid == NULL) {
      set_err(errno, "");
      return -1;
    }
    for (j = 0; j < n; j++) {
      matid[j] = start[j];
    }
    mat[i].matid = matid;
    mat[i].n     = n;
    nmat += n;
  }

  mesh->n_elem_mat_ID = nmat;
  if (mesh->n_elem_mat_ID > 0) {
    size_t size;
    size = sizeof(*mesh->elem_mat_ID_index) * (mesh->n_elem + 1);
    mesh->elem_mat_ID_index = HECMW_malloc(size);
    if (mesh->elem_mat_ID_index == NULL) {
      set_err(errno, "");
      return -1;
    }

    size                   = sizeof(*mesh->elem_mat_ID_item) * nmat;
    mesh->elem_mat_ID_item = HECMW_malloc(size);
    if (mesh->elem_mat_ID_item == NULL) {
      set_err(errno, "");
      return -1;
    }

    mesh->elem_mat_ID_index[0] = 0;
    for (i = 0; i < mesh->n_elem; i++) {
      mesh->elem_mat_ID_index[i + 1] = mesh->elem_mat_ID_index[i] + mat[i].n;
      for (j = 0; j < mat[i].n; j++) {
        int sidx                         = mesh->elem_mat_ID_index[i];
        mesh->elem_mat_ID_item[sidx + j] = mat[i].matid[j];
      }
    }
  }

  for (i = 0; i < mesh->n_elem; i++) {
    HECMW_free(mat[i].matid);
  }
  HECMW_free(mat);

  return 0;
}

static int setup_contact(struct hecmwST_local_mesh *mesh) {
  int i, npair, slave_gid, master_gid, orislave_sgid;
  size_t size;
  struct hecmwST_contact_pair *cpair;
  struct hecmw_io_contact *p;
  orislave_sgid = 0;
  slave_gid = 0;

  HECMW_assert(mesh);

  cpair = HECMW_malloc(sizeof(*cpair));
  if (cpair == NULL) {
    set_err(errno, "");
    goto error;
  }

  cpair->n_pair        = 0;
  cpair->type          = NULL;
  cpair->name          = NULL;
  cpair->slave_grp_id  = NULL;
  cpair->slave_orisgrp_id  = NULL;
  cpair->master_grp_id = NULL;

  if (_contact == NULL) {
    mesh->contact_pair = cpair;
    return 0;
  }

  /* count total # of contact pairs */
  npair = 0;
  for (p = _contact; p; p = p->next) {
    npair++;
  }
  HECMW_assert(npair > 0);
  cpair->n_pair = npair;

  /* name */
  size        = sizeof(*cpair->name) * (cpair->n_pair);
  cpair->name = HECMW_malloc(size);
  if (cpair->name == NULL) {
    set_err(errno, "");
    return -1;
  }

  /* type */
  size        = sizeof(*cpair->type) * (cpair->n_pair);
  cpair->type = HECMW_malloc(size);
  if (cpair->type == NULL) {
    set_err(errno, "");
    goto error;
  }

  /* slave_grp_id */
  size                = sizeof(*cpair->slave_grp_id) * (cpair->n_pair);
  cpair->slave_grp_id = HECMW_malloc(size);
  if (cpair->slave_grp_id == NULL) {
    set_err(errno, "");
    goto error;
  }

  /* slave_orisgrp_id */
  size                = sizeof(*cpair->slave_orisgrp_id) * (cpair->n_pair);
  cpair->slave_orisgrp_id = HECMW_malloc(size);
  if (cpair->slave_orisgrp_id == NULL) {
    set_err(errno, "");
    goto error;
  }

  /* master_grp_id */
  size                 = sizeof(*cpair->master_grp_id) * (cpair->n_pair);
  cpair->master_grp_id = HECMW_malloc(size);
  if (cpair->master_grp_id == NULL) {
    set_err(errno, "");
    goto error;
  }

  /* set */
  for (p = _contact, i = 0; p; p = p->next, i++) {
    HECMW_assert(i + 1 <= cpair->n_pair);

    HECMW_assert(strlen(p->name) < HECMW_NAME_LEN);
    cpair->name[i] = HECMW_strdup(p->name);
    if (cpair->name[i] == NULL) {
      set_err(errno, "");
      return -1;
    }

    cpair->type[i] = p->type;

    if (p->type == HECMW_CONTACT_TYPE_NODE_SURF) {
      slave_gid = HECMW_dist_get_ngrp_id(mesh->node_group, p->slave_grp);
      orislave_sgid = -1;
    } else if (p->type == HECMW_CONTACT_TYPE_SURF_SURF) {
      slave_gid = HECMW_dist_get_sgrp_id(mesh->surf_group, p->slave_grp);
      orislave_sgid = HECMW_dist_get_sgrp_id(mesh->surf_group, p->slave_orisgrp);
    } else if (p->type == HECMW_CONTACT_TYPE_NODE_ELEM) {
      slave_gid = HECMW_dist_get_ngrp_id(mesh->node_group, p->slave_grp);
      orislave_sgid = -1;
    } else {
      HECMW_assert(0);
    }
    HECMW_assert(slave_gid > 0);

    cpair->slave_grp_id[i] = slave_gid;
    cpair->slave_orisgrp_id[i] = orislave_sgid;

    if (p->type == HECMW_CONTACT_TYPE_NODE_ELEM) {
      master_gid = HECMW_dist_get_egrp_id(mesh->elem_group, p->master_grp);
    } else {
      master_gid = HECMW_dist_get_sgrp_id(mesh->surf_group, p->master_grp);
    }

    HECMW_assert(master_gid > 0);

    cpair->master_grp_id[i] = master_gid;
  }
  HECMW_assert(i == cpair->n_pair);

  mesh->contact_pair = cpair;

  return 0;

error:
  return -1;
}

static int setup_contact_sectid(struct hecmwST_local_mesh *mesh) {
  int i;
  struct hecmw_io_element *elem;

  HECMW_assert(mesh);

  for (i = 0; i < mesh->n_elem; i++) {
    /* depends on setup_elem() */
    if (mesh->elem_type[i] < 1000) continue;
    if (mesh->elem_type[i] >= 1100) continue;
    elem = HECMW_io_get_elem(mesh->global_elem_ID[i]);
    HECMW_assert(elem);
    HECMW_assert(elem->mpc_sectid != -1);
    mesh->section_ID[i] = elem->mpc_sectid;
  }
  return 0;
}

/*----------------------------------------------------------------------------*/

static int post_remove_unused_node(void) {
  int id;

  HECMW_assert(_node);

  /* used nodes have been marked in post_elem_check_node_existence() */

  HECMW_map_int_iter_init(_node);
  while (HECMW_map_int_iter_next_unmarked(_node, &id, NULL)) {
    /* remove from NODE GROUP */
    if (HECMW_io_remove_node_in_ngrp(id) < 0) {
      return -1;
    }
  }

  HECMW_map_int_del_unmarked(_node);

  return 0;
}

static int post_node(void) {
  int n_dup;

  if (_node == NULL || HECMW_map_int_nval(_node) == 0) {
    set_err(HECMW_IO_E1014, "");
    return -1;
  }

  n_dup = HECMW_map_int_check_dup(_node);
  if (n_dup > 0) {
    set_warn(HECMW_IO_W1004, "%d node(s) updated", n_dup);
  }

  return 0;
}

static int post_elem_check_node_existence(void) {
  int i, j, ncon, id;
  struct hecmw_io_element *p;

  HECMW_assert(global_node_ID_max > 0);

  if (HECMW_map_int_mark_init(_node)) {
    return -1;
  }

  HECMW_map_int_iter_init(_elem);
  for (i = 0; HECMW_map_int_iter_next(_elem, &id, (void **)&p); i++) {
    ncon = HECMW_get_max_node(p->type);

    HECMW_assert(ncon > 0);

    for (j = 0; j < ncon; j++) {
      HECMW_assert(p->node[j] > 0);
      HECMW_assert(p->node[j] <= global_node_ID_max);

      if (HECMW_map_int_mark(_node, p->node[j])) {
        set_err(HECMW_IO_E1027, "Node %d does not exist", p->node[j]);
        return -1;
      }
    }
  }

  return 0;
}

static char *post_elem_make_matname(int id, char *buf, int bufsize) {
  const char *matname = "HECMW-MAT";

  HECMW_assert(buf);
  HECMW_assert(bufsize > 0);

  sprintf(buf, "%s%d", matname, id);

  return buf;
}

static int post_elem_make_mat(void) {
  int i, j, id;
  char name[HECMW_NAME_LEN + 1];
  struct hecmw_io_element *p;
  struct hecmw_io_material *mat;
  struct hecmw_io_matitem *matitem;
  struct hecmw_io_matsubitem *matsubitem;

  HECMW_map_int_iter_init(_elem);
  for (i = 0; HECMW_map_int_iter_next(_elem, &id, (void **)&p); i++) {
    if (p->nmatitem <= 0) continue;

    mat = HECMW_malloc(sizeof(*mat));
    if (mat == NULL) {
      set_err(errno, "");
      return -1;
    }

    matitem = HECMW_malloc(sizeof(*matitem));
    if (matitem == NULL) {
      set_err(errno, "");
      return -1;
    }

    matsubitem = HECMW_malloc(sizeof(*matsubitem));
    if (matsubitem == NULL) {
      set_err(errno, "");
      return -1;
    }

    matsubitem->val = HECMW_malloc(sizeof(*matsubitem->val) * p->nmatitem);
    if (matsubitem->val == NULL) {
      set_err(errno, "");
      return -1;
    }

    for (j = 0; j < p->nmatitem; j++) {
      matsubitem->val[j] = p->matitem[j];
    }
    matsubitem->temp = 0.0;
    matsubitem->next = NULL;

    matitem->item    = 1;
    matitem->nval    = p->nmatitem;
    matitem->subitem = matsubitem;

    mat->nitem = 1;
    mat->item  = matitem;
    post_elem_make_matname(id, name, sizeof(name));
    strcpy(p->matname, name);
    strcpy(mat->name, name);
    mat->next = NULL;

    if (HECMW_io_add_mat(name, mat) == NULL) return -1;
  }
  return 0;
}

static int post_elem(void) {
  int n_dup;

  if (_elem == NULL) {
    set_err(HECMW_IO_E1015, "");
    return -1;
  }

  n_dup = HECMW_map_int_check_dup(_elem);
  if (n_dup > 0) {
    set_warn(HECMW_IO_W1001, "%d element(s) updated", n_dup);
  }

  if (post_elem_check_node_existence()) return -1;

  if (post_elem_make_mat()) return -1;

  return 0;
}

static int post_ngrp(void) {
  struct hecmw_io_ngrp *p;

  for (p = _ngrp; p; p = p->next) {
    int n_dup, id, i;

    n_dup = HECMW_set_int_check_dup(p->node);
    if (n_dup > 0) set_warn(HECMW_IO_W1006, "%d node(s) in %s", n_dup, p->name);

    HECMW_set_int_iter_init(p->node);
    for (i = 0; HECMW_set_int_iter_next(p->node, &id); i++) {
      if (HECMW_io_get_node(id) == NULL) {
        set_warn(HECMW_IO_W1005, "Node %d doesn't exist", id);
        HECMW_set_int_del(p->node, id);
      }
    }
  }
  return 0;
}

static int post_egrp(void) {
  struct hecmw_io_egrp *p;

  for (p = _egrp; p; p = p->next) {
    int n_dup, id, i;

    n_dup = HECMW_set_int_check_dup(p->elem);
    if (n_dup > 0)
      set_warn(HECMW_IO_W1003, "%d element(s) in %s", n_dup, p->name);

    HECMW_set_int_iter_init(p->elem);
    for (i = 0; HECMW_set_int_iter_next(p->elem, &id); i++) {
      if (HECMW_io_get_elem(id) == NULL) {
        set_warn(HECMW_IO_W1002, "Element %d doesn't exist", id);
        HECMW_set_int_del(p->elem, id);
      }
    }
  }
  return 0;
}

static int post_sgrp(void) {
  struct hecmw_io_sgrp *p;

  for (p = _sgrp; p; p = p->next) {
    int n_dup, id, i;

    n_dup = HECMW_set_int_check_dup(p->item);
    if (n_dup > 0)
      set_warn(HECMW_IO_W1009, "%d surface(s) in %s", n_dup, p->name);

    HECMW_set_int_iter_init(p->item);
    for (i = 0; HECMW_set_int_iter_next(p->item, &id); i++) {
      int eid, sid;
      struct hecmw_io_element *element;

      decode_surf_key(id, &eid, &sid);

      /* check element */
      element = HECMW_io_get_elem(eid);
      if (element == NULL) {
        set_warn(HECMW_IO_W1007, "Element %d doesn't exist", eid);
        HECMW_set_int_del(p->item, id);
        continue;
      }

      /* check surface */
      if (HECMW_get_max_surf(element->type) < sid) {
        set_warn(HECMW_IO_W1008, "Element %d, surface %d", eid, sid);
        HECMW_set_int_del(p->item, id);
      }
    }
  }
  return 0;
}

static int post_initial_check_node_exists(void) {
  int ignore;
  struct hecmw_io_initial *p, *prev, *next;

  if (_init == NULL) return 0;

  /* check node existence */
  prev = NULL;
  for (p = _init; p; p = next) {
    next = p->next;
    if (p->node == -1) {
      if (prev) {
        prev->next = p;
      }
      prev = p;
      continue;
    }

    ignore = 0;
    if (HECMW_io_get_node(p->node) == NULL) {
      set_warn(HECMW_IO_W1016, "Node %d does not eixist", p->node);
      ignore = 1;
    }
    if (ignore) {
      HECMW_free(p);
      if (prev == NULL) {
        _init = next;
      } else {
        prev->next = next;
      }
    } else {
      if (prev == NULL) {
        _init = p;
      } else {
        prev->next = p;
      }
      prev = p;
    }
  }
  return 0;
}

static int post_initial_ngrp_to_node(void) {
  int i, nnode, ignore, *node;
  struct hecmw_io_initial *p, *prev, *next, *new_init;
  struct hecmw_io_id_array *id;

  if (_init == NULL) return 0;

  /* change ngrp to node */
  prev = NULL;
  for (p = _init; p; p = next) {
    next = p->next;
    if (p->node != -1) {
      if (prev) {
        prev->next = p;
      }
      prev = p;
      continue;
    }

    /* check existence */
    ignore = 0;
    if (HECMW_io_get_ngrp(p->ngrp) == NULL) {
      set_warn(HECMW_IO_W1017, "Node group %s does not eixist", p->ngrp);
      ignore = 1;
    }
    if (ignore) {
      HECMW_free(p);
      if (prev == NULL) {
        _init = next;
      } else {
        prev->next = next;
      }
      continue;
    }

    /* check # of node in node group */
    ignore = 0;
    nnode  = HECMW_io_get_nnode_in_ngrp(p->ngrp);
    HECMW_assert(nnode > 0);

    /* replace by node */
    id = HECMW_io_get_node_in_ngrp(p->ngrp);
    HECMW_assert(id);
    HECMW_assert(id->n == nnode);
    node = id->id;
    HECMW_free(id);

    for (i = 0; i < nnode; i++) {
      new_init = HECMW_malloc(sizeof(*new_init));
      if (new_init == NULL) {
        set_err(errno, "");
        return -1;
      }
      memcpy(new_init, p, sizeof(*new_init));
      new_init->next = NULL;
      new_init->node = node[i];

      if (prev == NULL) {
        _init = new_init;
      } else {
        prev->next = new_init;
      }
      prev = new_init;
    }

    HECMW_free(node);
    HECMW_free(p);
  }
  return 0;
}

static int post_initial_check_dup(void) {
  struct hecmw_io_initial *p;
  struct hecmw_set_int set;
  int ndup;

  if (_init == NULL) return 0;

  HECMW_set_int_init(&set);

  /* check duplication */
  for (p = _init; p; p = p->next) {
    HECMW_set_int_add(&set, p->node);
  }
  ndup = HECMW_set_int_check_dup(&set);

  HECMW_set_int_finalize(&set);

  if (ndup > 0) {
    set_err(HECMW_IO_E1018, "Some nodes are initialized more than once");
    return -1;
  }
  return 0;
}

static int post_initial(void) {
  if (_init == NULL) return 0;

  if (post_initial_check_node_exists()) return -1;
  HECMW_log(HECMW_LOG_DEBUG, "post_initial_check_node_exists done");
  if (post_initial_ngrp_to_node()) return -1;
  HECMW_log(HECMW_LOG_DEBUG, "post_initial_ngrp_to_node done");
  if (post_initial_check_dup()) return -1;
  HECMW_log(HECMW_LOG_DEBUG, "post_initial_check_dup done");

  return 0;
}

static int post_equation_check_node_exists(void) {
  int i, ignore;
  struct hecmw_io_mpc *p, *prev, *next;

  if (_mpc == NULL) return 0;

  /* check node existence */
  prev = NULL;
  for (p = _mpc; p; p = next) {
    next = p->next;
    if (p->item[0].node == -1) {
      if (prev) {
        prev->next = p;
      }
      prev = p;
      continue;
    }

    ignore = 0;
    HECMW_assert(p->neq >= 2);
    for (i = 0; i < p->neq; i++) {
      struct hecmw_io_mpcitem *item = &p->item[i];
      if (HECMW_io_get_node(item->node) == NULL) {
        set_warn(HECMW_IO_W1019, "Node %d not found", item->node);
        ignore = 1;
        break;
      }
    }
    if (ignore) {
      HECMW_free(p->item);
      HECMW_free(p);
      if (prev == NULL) {
        _mpc = next;
      } else {
        prev->next = next;
      }
      continue;
    }
  }
  return 0;
}

static int post_equation_ngrp_to_node(void) {
  int i, j, ignore, **node;
  struct hecmw_io_mpc *p, *prev, *next, *new_mpc;

  if (_mpc == NULL) return 0;

  /* change ngrp to node */
  prev = NULL;
  for (p = _mpc; p; p = next) {
    int nnode;
    next = p->next;
    if (p->item[0].node != -1) {
      if (prev) {
        prev->next = p;
      }
      prev = p;
      continue;
    }

    /* check existence */
    ignore = 0;
    HECMW_assert(p->neq >= 2);
    for (i = 0; i < p->neq; i++) {
      struct hecmw_io_mpcitem *item = &p->item[i];
      HECMW_assert(item->node == -1);
      if (HECMW_io_get_ngrp(item->ngrp) == NULL) {
        set_warn(HECMW_IO_W1020, "Node group %s not found", item->ngrp);
        ignore = 1;
        break;
      }
    }
    if (ignore) {
      HECMW_free(p->item);
      HECMW_free(p);
      if (prev == NULL) {
        _mpc = next;
      } else {
        prev->next = next;
      }
      continue;
    }

    /* check # of node in node group */
    ignore = 0;
    nnode  = HECMW_io_get_nnode_in_ngrp(p->item[0].ngrp);
    HECMW_assert(nnode > 0);
    for (i = 1; i < p->neq; i++) {
      struct hecmw_io_mpcitem *item = &p->item[i];
      int n                         = HECMW_io_get_nnode_in_ngrp(item->ngrp);
      if (n != nnode) {
        set_err(HECMW_IO_E1021, "%d node%s in %s, %d node%s in %s", nnode,
                (nnode != 0) ? "s" : "", p->item[0].ngrp, n,
                (n != 0) ? "s" : "", p->item[i].ngrp);
        return -1;
      }
    }

    /* replace by node */
    node = HECMW_malloc(sizeof(node) * p->neq);
    if (node == NULL) {
      set_err(errno, "");
      return -1;
    }

    for (i = 0; i < p->neq; i++) {
      struct hecmw_io_id_array *id = HECMW_io_get_node_in_ngrp(p->item[i].ngrp);
      HECMW_assert(id);
      HECMW_assert(id->n == nnode);
      node[i] = id->id;
      HECMW_free(id);
    }

    for (i = 0; i < nnode; i++) {
      new_mpc = HECMW_malloc(sizeof(*new_mpc));
      if (new_mpc == NULL) {
        set_err(errno, "");
        return -1;
      }
      memcpy(new_mpc, p, sizeof(*new_mpc));
      new_mpc->next = NULL;
      new_mpc->item = NULL;

      new_mpc->item = HECMW_malloc(sizeof(*new_mpc->item) * (p->neq));
      if (new_mpc == NULL) {
        set_err(errno, "");
        return -1;
      }

      for (j = 0; j < new_mpc->neq; j++) {
        struct hecmw_io_mpcitem *item = &new_mpc->item[j];
        item->node                    = node[j][i];
        item->dof                     = p->item[j].dof;
        item->a                       = p->item[j].a;
      }

      if (prev == NULL) {
        _mpc = new_mpc;
      } else {
        prev->next = new_mpc;
      }
      prev = new_mpc;
    }

    for (i = 0; i < p->neq; i++) {
      HECMW_free(node[i]);
    }
    HECMW_free(node);

    HECMW_free(p->item);
    HECMW_free(p);
  }
  return 0;
}

/*
 * must be node(not allow ngrp)
 */
static int post_equation_check_dup(void) {
  int i;
  struct hecmw_io_mpc *p, *q;

  if (_mpc == NULL) return 0;

  /* check duplication */
  for (p = _mpc; p; p = p->next) {
    int nod = p->item[0].node;
    int dof = p->item[0].dof;
    for (q = _mpc; q; q = q->next) {
      HECMW_assert(q->neq >= 2);
      for (i = 1; i < q->neq; i++) {
        HECMW_assert(q->item[i].node != -1);
        if (q->item[i].node == nod && q->item[i].dof == dof) {
          set_err(HECMW_IO_E1022, "Node:%d and DOF:%d", nod, dof);
          return -1;
        }
      }
    }
  }
  return 0;
}

static int post_equation_add_elem(void) {
  int i, j, mpc_id, elem_id, dof1, dof2, type;
  int node[2];
  struct hecmw_io_mpc *p;
  struct hecmw_io_element *elem;

  if (_mpc == NULL) return 0;

  /* max element ID */
  elem_id = HECMW_io_get_elem_max_id();
  elem_id++;

  /* add element */
  for (p = _mpc, mpc_id = 1; p; p = p->next, mpc_id++) {
    HECMW_assert(p->neq >= 2);
    for (j = 0; j < p->neq - 1; j++) {
      dof1 = p->item[j].dof;
      for (i = j + 1; i < p->neq; i++) {
        dof2 = p->item[i].dof;
        /* make element type */
        type = 900 + dof1 * 10 + dof2;
        HECMW_assert(HECMW_get_max_node(type) == 2);

        /* set node */
        node[0] = p->item[j].node;
        node[1] = p->item[i].node;

        /* add */
        elem = HECMW_io_add_elem(elem_id, type, node, 0, NULL);
        if (elem == NULL) {
          return -1;
        }

        elem->mpc_matid  = (j + 1) * 10 + (i + 1);
        elem->mpc_sectid = mpc_id;

        if (HECMW_io_add_egrp("ALL", 1, &elem_id) < 0) {
          return -1;
        }
        elem_id++;
      }
    }
  }
  return 0;
}

static int post_equation(void) {
  if (_mpc == NULL) return 0;

  if (post_equation_check_node_exists()) return -1;
  if (post_equation_ngrp_to_node()) return -1;
  /* Delete because performance grow worse at large number of equations
          if(post_equation_check_dup()) return -1;
  */
  if (post_equation_add_elem()) return -1;

  return 0;
}

static int post_section_check_exists(void) {
  if (_sect == NULL) {
    set_err(HECMW_IO_E1023, "");
    return -1;
  }
  return 0;
}

static int post_section_check_egrp(void) {
  int i, eid;
  struct hecmw_io_section *p;

  for (p = _sect; p; p = p->next) {
    struct hecmw_io_egrp *egrp = HECMW_io_get_egrp(p->egrp);

    if (egrp == NULL) {
      set_err(HECMW_IO_E1024, "Element group %s not found", p->egrp);
      goto error;
    }

    HECMW_set_int_iter_init(egrp->elem);
    for (i = 0; HECMW_set_int_iter_next(egrp->elem, &eid); i++) {
      struct hecmw_io_element *elem = HECMW_io_get_elem(eid);

      HECMW_assert(elem);

      if (HECMW_is_etype_link(elem->type)) continue;

      switch (p->type) {
        case HECMW_SECT_TYPE_SHELL:
          if (!HECMW_is_etype_shell(elem->type)) {
            set_err(HECMW_IO_E1026, "Only shell element allowed in EGRP %s",
                    p->egrp);
            goto error;
          }
          break;
        case HECMW_SECT_TYPE_SOLID:
          if (!HECMW_is_etype_solid(elem->type)) {
            set_err(HECMW_IO_E1026, "Only solid element allowed in EGRP %s",
                    p->egrp);
            goto error;
          }
          break;
        case HECMW_SECT_TYPE_BEAM:
          if (!HECMW_is_etype_beam(elem->type)) {
            set_err(HECMW_IO_E1026, "Only beam element allowed in EGRP %s",
                    p->egrp);
            goto error;
          }
          break;
        case HECMW_SECT_TYPE_INTERFACE:
          if (!HECMW_is_etype_interface(elem->type)) {
            set_err(HECMW_IO_E1026, "Only interface element allowed in EGRP %s",
                    p->egrp);
            goto error;
          }
          break;
        default:
          HECMW_assert(0);
      }
    }
  }
  return 0;
error:
  return -1;
}

static int post_section_check_mat_exists(void) {
  int found;
  struct hecmw_io_section *p;
  struct hecmw_io_material *mat;
  extern hecmw_hash_p *hash_mat;

  for (p = _sect; p; p = p->next) {
    found = 0;
    /* for(mat=_mat; mat; mat=mat->next) {
            if(strcmp(p->material, mat->name) == 0) {
                    found = 1;
                    break;
            }
    }*/
    if ((struct hecmw_io_material *)hecmw_hash_p_get(hash_mat, p->material) !=
        NULL) {
      found = 1;
    }
    if (p->type != HECMW_SECT_TYPE_INTERFACE && !found) {
      set_err(HECMW_IO_E1025, "MATERIAL %s not found", p->material);
      return -1;
    }
  }
  return 0;
}

static int post_section(void) {
  if (post_section_check_exists()) return -1;
  if (post_section_check_egrp()) return -1;
  if (post_section_check_mat_exists()) return -1;

  return 0;
}

static int post_contact_check_grp(void) {
  int i;
  struct hecmw_io_contact *p;

  for (p = _contact; p; p = p->next) {
    struct hecmw_io_ngrp *ngrp;
    struct hecmw_io_sgrp *sgrp;
    struct hecmw_io_egrp *egrp;

    if (p->type == HECMW_CONTACT_TYPE_NODE_SURF) {
      ngrp = HECMW_io_get_ngrp(p->slave_grp);
      if (ngrp == NULL) {
        set_err(HECMW_IO_E1029, "Node group %s not found", p->slave_grp);
        goto error;
      }
      sgrp = HECMW_io_get_sgrp(p->master_grp);
      if (sgrp == NULL) {
        set_err(HECMW_IO_E1028, "Surface group %s not found", p->master_grp);
        goto error;
      }
    } else if (p->type == HECMW_CONTACT_TYPE_SURF_SURF) {
      sgrp = HECMW_io_get_sgrp(p->slave_grp);
      if (sgrp == NULL) {
        set_err(HECMW_IO_E1028, "Surface group %s not found", p->slave_grp);
        goto error;
      }
      sgrp = HECMW_io_get_sgrp(p->master_grp);
      if (sgrp == NULL) {
        set_err(HECMW_IO_E1028, "Surface group %s not found", p->master_grp);
        goto error;
      }
    } else if (p->type == HECMW_CONTACT_TYPE_NODE_ELEM) {
      ngrp = HECMW_io_get_ngrp(p->slave_grp);
      if (ngrp == NULL) {
        set_err(HECMW_IO_E1029, "Node group %s not found", p->slave_grp);
        goto error;
      }
      egrp = HECMW_io_get_egrp(p->master_grp);
      if (egrp == NULL) {
        set_err(HECMW_IO_E1028, "Element group %s not found", p->master_grp);
        goto error;
      }
    } else {
      fprintf(stderr, "ERROR: CONTACT_PAIR: TYPE=%d\n", p->type);
      HECMW_assert(0);
    }

  }
  return 0;
error:
  return -1;
}

static int post_contact_convert_sgroup(void)
{
  struct hecmw_io_contact *p;
  int elem_id, contact_id;

  elem_id = HECMW_io_get_elem_max_id();
  elem_id++;

  for (p = _contact, contact_id = 1; p; p = p->next, contact_id++) {
    struct hecmw_io_sgrp *sgrp;
    int n_item, i, id, ret;
    int *elem, *surf;
    char new_sgrp_name[HECMW_NAME_LEN+1];

    if (p->type != HECMW_CONTACT_TYPE_SURF_SURF) continue;

    sgrp = HECMW_io_get_sgrp(p->slave_grp);
    HECMW_assert(sgrp);

    n_item = HECMW_set_int_nval(sgrp->item);
    if (n_item == 0) continue;

    elem = (int *) malloc(sizeof(int) * n_item);
    surf = (int *) malloc(sizeof(int) * n_item);
    if (!elem || !surf) {
      set_err(errno, "");
      return -1;
    }

    HECMW_set_int_iter_init(sgrp->item);
    for (i = 0; HECMW_set_int_iter_next(sgrp->item, &id); i++) {
      int eid, sid, etype, surf_etype, surf_nnode, j;
      struct hecmw_io_element *element, *ptelem;
      const int *surf_nodes;
      int nodes[8];

      decode_surf_key(id, &eid, &sid);

      element = HECMW_io_get_elem(eid);
      HECMW_assert(element);
      etype = element->type;

      /* extract surface */
      surf_nodes = HECMW_get_surf_nodes(etype, sid, &surf_etype);
      HECMW_assert( HECMW_is_etype_patch(surf_etype) );

      surf_nnode = HECMW_get_max_node(surf_etype);

      for (j = 0; j < surf_nnode; j++) {
        nodes[j] = element->node[surf_nodes[j] - 1];
      }

      /* add surface patch elem */
      ptelem = HECMW_io_add_elem(elem_id, surf_etype, nodes, 0, NULL);
      if (ptelem == NULL) {
        return -1;
      }

      ptelem->mpc_matid = surf_etype % 100;
      ptelem->mpc_sectid = contact_id;

      elem[i] = elem_id;
      surf[i] = 1;

      elem_id++;
    }

    /* add newly added patch elems to egrp "ALL" */
    if (HECMW_io_add_egrp("ALL", n_item, elem) < 0)
      return -1;

    /* generate name for new sgrp with patch elems */
    ret = snprintf(new_sgrp_name, sizeof(new_sgrp_name), "_PT_%s", sgrp->name);
    if (ret >= sizeof(new_sgrp_name)) {
      set_err(HECMW_IO_E0001, "Surface group name: %s", sgrp->name);
      return -1;
    } else if (HECMW_io_get_sgrp(new_sgrp_name) != NULL) {
      set_err(HECMW_IO_E0003, "Surface group name: %s", new_sgrp_name);
      return -1;
    }

    /* add sgrp with patch elems */
    if (HECMW_io_add_sgrp(new_sgrp_name, n_item, elem, surf) < 0)
      return -1;

    free(elem);
    free(surf);

    /* replace slave group by newly added sgrp with patch elems */
    strcpy(p->slave_grp, new_sgrp_name);
  }
  return 0;
}

static int post_contact(void) {
  if (post_contact_check_grp()) return -1;
  if (post_contact_convert_sgroup()) return -1;

  return 0;
}

/*----------------------------------------------------------------------------*/

int HECMW_io_check_mpc_dof(int dof) {
  if (dof < 1 || dof > 6) return -1;
  return 0;
}

int HECMW_io_is_reserved_name(const char *name) {
  if (name == NULL) return 0;
  if (strncmp("HECMW", name, 5) == 0) return 1;
  return 0;
}

int HECMW_io_get_version(void) { return HECMW_FLAG_VERSION; }

int HECMW_hash_init(void) {
  extern hecmw_hash_p *hash_ng;
  extern hecmw_hash_p *hash_eg;
  extern hecmw_hash_p *hash_sg;
  extern hecmw_hash_p *hash_mat;

  hash_ng = hecmw_hash_p_new(1);
  if (hash_ng == NULL) return 1;
  hash_eg = hecmw_hash_p_new(2);
  if (hash_eg == NULL) return 1;
  hash_sg = hecmw_hash_p_new(3);
  if (hash_sg == NULL) return 1;
  hash_mat = hecmw_hash_p_new(4);
  if (hash_mat == NULL) return 1;

  return 0;
}

int HECMW_hash_finalize(void) {
  extern hecmw_hash_p *hash_ng;
  extern hecmw_hash_p *hash_eg;
  extern hecmw_hash_p *hash_sg;
  extern hecmw_hash_p *hash_mat;

  hecmw_hash_p_delete(hash_ng);
  hecmw_hash_p_delete(hash_eg);
  hecmw_hash_p_delete(hash_sg);
  hecmw_hash_p_delete(hash_mat);

  return 0;
}

int HECMW_io_init(void) {
  HECMW_log(HECMW_LOG_DEBUG, "Initializing IO process...");

  if (HECMW_hash_init()) {
    printf("ERROE:HECMW_HASHTABLE INIT \n");
    return -1;
  }
  if (clear()) {
    return -1;
  }

  return 0;
}

int HECMW_io_finalize(void) {
  HECMW_log(HECMW_LOG_DEBUG, "Finalizing IO process...");

  if (HECMW_hash_finalize()) {
    printf("ERROE:HECMW_HASHTABLE FINALIZE \n");
    return -1;
  }
  if (clear()) {
    return -1;
  }

  return 0;
}

int HECMW_io_pre_process(void) { return 0; }

int HECMW_io_post_process(void) {
  HECMW_log(HECMW_LOG_DEBUG, "Running post process...");

  if (post_node()) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "post_node done");
  if (post_elem()) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "post_elem done");
  if (post_ngrp()) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "post_ngrp done");
  if (post_egrp()) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "post_egrp done");
  if (post_sgrp()) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "post_sgrp done");
  if (post_remove_unused_node()) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "post_remove_unused_node done");
  if (post_initial()) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "post_initial done");
  if (post_equation()) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "post_equation done");
  if (post_section()) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "post_section done");
  if (post_contact()) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "post_contact done");
  return 0;
error:
  return -1;
}

struct hecmwST_local_mesh *HECMW_io_make_local_mesh(void) {
  struct hecmwST_local_mesh *mesh;

  HECMW_log(HECMW_LOG_DEBUG, "Creating hecmwST_local_mesh...");

  mesh = HECMW_calloc(1, sizeof(*mesh));
  if (mesh == NULL) {
    set_err(errno, "");
    goto error;
  }

  if (setup_flags(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_flags done");
  if (setup_gridfile(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_gridfile done");
  if (setup_files(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_files done");
  if (setup_header(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_header done");
  if (setup_zero(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_zero done");
  if (setup_node(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_node done");
  if (setup_init(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_init done");
  if (setup_elem(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_elem done");
  if (setup_ngrp(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_ngrp done");
  if (setup_egrp(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_egrp done");
  if (setup_sgrp(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_sgrp done");
  if (setup_pe(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_pe done");
  if (setup_adapt(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_adapt done");
  if (setup_refine(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_refine done");
  if (setup_mpc(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_mpc done");
  if (setup_amp(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_amp done");
  if (setup_mat(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_mat done");
  if (setup_sect(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_sect done");
  if (setup_mpc_sectid(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_mpc_sectid done");
  if (setup_contact_sectid(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_contact_sectid done");
  if (setup_elem_check_sectid(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_elem_check_sectid done");
  if (setup_elem_mat(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_elem_mat done");
  if (setup_mpc_reorder(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_mpc_reorder done");
  if (setup_contact(mesh)) goto error;
  HECMW_log(HECMW_LOG_DEBUG, "setup_contact done");
  return mesh;
error:
  return NULL;
}
