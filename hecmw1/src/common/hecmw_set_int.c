/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.1                                               *
 *                                                                     *
 *     Last Update : 2007/06/29                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuya Goto (AdvanceSoft)                     *
 *                                                                     *
 *     Contact address :  IIS, The University of Tokyo RSS21 project   *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "hecmw_util.h"
#include "hecmw_malloc.h"
#include "hecmw_config.h"
#include "hecmw_set_int.h"

enum { SET_MAX_VAL_INIT = 64, SET_MAX_VAL_GROW = 2 };


int HECMW_set_int_init(struct hecmw_set_int *set)
{
  HECMW_assert(set);

  set->n_val = 0;
  set->max_val = 0;

  set->vals = NULL;

  set->checked = 1;
  set->sorted = 1;

  set->iter = -1;

  return HECMW_SUCCESS;
}

void HECMW_set_int_finalize(struct hecmw_set_int *set)
{
  HECMW_assert(set);

  if (set->max_val == 0) {
    HECMW_assert(set->n_val == 0);
    return;
  }

  HECMW_free(set->vals);
  return;
}


int HECMW_set_int_nval(const struct hecmw_set_int *set)
{
  HECMW_assert(set);

  return set->n_val;
}

static int set_resize(struct hecmw_set_int *set, int new_max_val)
{
  int *new_vals;

  HECMW_assert(set);
  HECMW_assert(set->n_val <= new_max_val);

  if (set->max_val == new_max_val) return HECMW_SUCCESS;

  if (new_max_val == 0) {
    HECMW_assert(set->vals);

    free(set->vals);
    set->vals = NULL;
    set->max_val = 0;

    return HECMW_SUCCESS;
  }

  new_vals = (int *) HECMW_realloc(set->vals, sizeof(int) * new_max_val);
  if (new_vals == NULL) {
  	return HECMW_ERROR;
  }

  set->vals = new_vals;
  set->max_val = new_max_val;

  return HECMW_SUCCESS;
}

static int set_grow(struct hecmw_set_int *set)
{
  int new_max_val;

  HECMW_assert(set);

  if (set->max_val == 0)
    new_max_val = SET_MAX_VAL_INIT;
  else
    new_max_val = set->max_val * SET_MAX_VAL_GROW;

  return set_resize(set, new_max_val); 
}

int HECMW_set_int_add(struct hecmw_set_int *set, int value)
{
  HECMW_assert(set);

  if (set->n_val == set->max_val)
    if (set_grow(set) != HECMW_SUCCESS)
      return HECMW_ERROR;

  set->vals[set->n_val] = value;

  if (set->n_val > 0 && set->sorted) {
    int val_prev = set->vals[set->n_val - 1];

    if (val_prev > value)
      set->sorted = set->checked = 0;

    if (set->checked && val_prev == value)
      set->checked = 0;
  }

  set->n_val++;

  return HECMW_SUCCESS;
}

static int int_cmp(const void *v1, const void *v2)
{
  const int *i1, *i2;

  i1 = (const int *) v1;
  i2 = (const int *) v2;

  if (*i1 < *i2) return -1;
  if (*i1 > *i2) return 1;
  return 0;
}

int HECMW_set_int_check_dup(struct hecmw_set_int *set)
{
  int i, n_dup = 0;

  HECMW_assert(set);

  if (set->checked) return 0;

  if (!set->sorted) {
    qsort(set->vals, set->n_val, sizeof(int), int_cmp);
    set->sorted = 1;
  }

  for (i = 1; i < set->n_val; i++) {
    if (set->vals[i-1] == set->vals[i]) {
      n_dup++;
    } else {
      if (n_dup > 0) {
        set->vals[i-n_dup] = set->vals[i];
      }
    }
  }

  set->n_val -= n_dup;
  set->checked = 1;

  set_resize(set, set->n_val); /* reduce memory usage */

  return n_dup;
}

static int set_del_index(struct hecmw_set_int *set, int index)
{
  HECMW_assert(set);
  HECMW_assert(0 <= index && index < set->n_val);

  if (index < set->n_val - 1)
    memmove(set->vals + index,
            set->vals + index + 1, sizeof(int) * (set->n_val - index - 1));

  set->n_val--;

  if (index < set->iter) set->iter--;

  return HECMW_SUCCESS;
}

int HECMW_set_int_del(struct hecmw_set_int *set, int value)
{
  int *p;

  HECMW_assert(set);
  HECMW_assert(set->checked);

  p = bsearch(&value, set->vals, set->n_val, sizeof(int), int_cmp);
  if (p == NULL) return HECMW_ERROR;

  return set_del_index(set, p - set->vals);
}


void HECMW_set_int_iter_init(struct hecmw_set_int *set)
{
  HECMW_assert(set);
  HECMW_assert(set->checked);

  set->iter = 0;
  return;
}

int HECMW_set_int_iter_next(struct hecmw_set_int *set, int *value)
{
  HECMW_assert(set);
  HECMW_assert(0 <= set->iter && set->iter <= set->n_val);

  if (set->iter == set->n_val) {
  	set->iter = -1;
  	return 0;
  }

  *value = set->vals[set->iter];
  set->iter++;

  return 1;
}
