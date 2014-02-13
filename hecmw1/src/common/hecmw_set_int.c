/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2013/12/18                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuya Goto (PExProCS)                        *
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
#include "hecmw_varray_int.h"
#include "hecmw_set_int.h"


int HECMW_set_int_init(struct hecmw_set_int *set)
{
  HECMW_assert(set);

  set->vals = (struct hecmw_varray_int *) HECMW_malloc(sizeof(struct hecmw_varray_int));
  if (set->vals == NULL) {
    return HECMW_ERROR;
  }

  if (HECMW_varray_int_init(set->vals) != HECMW_SUCCESS)
    return HECMW_ERROR;

  set->checked = 1;
  set->sorted = 1;

  set->in_iter = 0;
  set->iter = 0;

  return HECMW_SUCCESS;
}

void HECMW_set_int_finalize(struct hecmw_set_int *set)
{
  HECMW_assert(set);

  HECMW_varray_int_finalize(set->vals);
  HECMW_free(set->vals);

  return;
}


size_t HECMW_set_int_nval(struct hecmw_set_int *set)
{
  HECMW_assert(set);

  if (!set->checked) {
    HECMW_set_int_check_dup(set);
  }
  return HECMW_varray_int_nval(set->vals);
}

int HECMW_set_int_is_empty(const struct hecmw_set_int *set)
{
  HECMW_assert(set);

  return HECMW_varray_int_nval(set->vals) == 0 ? 1 : 0;
}

int HECMW_set_int_add(struct hecmw_set_int *set, int value)
{
  HECMW_assert(set);

  size_t nval = HECMW_varray_int_nval(set->vals);
  if (nval > 0 && set->sorted) {
    int val_prev = HECMW_varray_int_get(set->vals, nval-1);

    if (val_prev > value)
      set->sorted = set->checked = 0;

    if (set->checked && val_prev == value)
      set->checked = 0;
  }

  if (HECMW_varray_int_append(set->vals, value) != HECMW_SUCCESS)
    return HECMW_ERROR;

  return HECMW_SUCCESS;
}

size_t HECMW_set_int_check_dup(struct hecmw_set_int *set)
{
  size_t i, n_dup = 0;

  HECMW_assert(set);

  if (set->checked) return 0;

  if (!set->sorted) {
    HECMW_varray_int_sort(set->vals);
    set->sorted = 1;
  }

  n_dup = HECMW_varray_int_uniq(set->vals);
  set->checked = 1;

  return n_dup;
}

int HECMW_set_int_del(struct hecmw_set_int *set, int value)
{
  size_t index;

  HECMW_assert(set);

  if (!set->checked) {
    HECMW_assert(!set->in_iter);
    HECMW_set_int_check_dup(set);
  }

  if (HECMW_varray_int_search(set->vals, value, &index) != HECMW_SUCCESS)
    return HECMW_ERROR;

  HECMW_varray_int_delete(set->vals, index);

  if (index < set->iter) set->iter--;

  return HECMW_SUCCESS;
}


void HECMW_set_int_iter_init(struct hecmw_set_int *set)
{
  HECMW_assert(set);

  if (!set->checked) {
    HECMW_set_int_check_dup(set);
  }

  set->in_iter = 1;
  set->iter = 0;
  return;
}

int HECMW_set_int_iter_next(struct hecmw_set_int *set, int *value)
{
  size_t nval = HECMW_varray_int_nval(set->vals);
  HECMW_assert(set);
  HECMW_assert(set->in_iter);
  HECMW_assert(set->iter <= nval);

  if (set->iter == nval) {
    set->in_iter = 0;
    set->iter = 0;
    return 0;
  }

  *value = HECMW_varray_int_get(set->vals, set->iter);
  set->iter++;

  return 1;
}
