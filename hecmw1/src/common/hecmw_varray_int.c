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
#include "hecmw_bit_array.h"
#include "hecmw_varray_int.h"

enum { VARRAY_MAX_VAL_INIT = 64, VARRAY_MAX_VAL_GROW = 2 };


int HECMW_varray_int_init(struct hecmw_varray_int *varray)
{
  HECMW_assert(varray);

  varray->n_val = 0;
  varray->max_val = 0;

  varray->vals = NULL;

  return HECMW_SUCCESS;
}

void HECMW_varray_int_finalize(struct hecmw_varray_int *varray)
{
  HECMW_assert(varray);

  if (varray->max_val == 0) {
    HECMW_assert(varray->n_val == 0);
    return;
  }

  HECMW_free(varray->vals);
  return;
}


size_t HECMW_varray_int_nval(const struct hecmw_varray_int *varray)
{
  HECMW_assert(varray);

  return varray->n_val;
}

static int varray_resize(struct hecmw_varray_int *varray, size_t new_max_val)
{
  int *new_vals;

  HECMW_assert(varray);
  HECMW_assert(varray->n_val <= new_max_val);

  if (varray->max_val == new_max_val) return HECMW_SUCCESS;

  if (new_max_val == 0) {
    HECMW_assert(varray->vals);

    free(varray->vals);
    varray->vals = NULL;
    varray->max_val = 0;

    return HECMW_SUCCESS;
  }

  new_vals = (int *) HECMW_realloc(varray->vals, sizeof(int) * new_max_val);
  if (new_vals == NULL) {
    return HECMW_ERROR;
  }

  varray->vals = new_vals;
  varray->max_val = new_max_val;

  return HECMW_SUCCESS;
}

static int varray_grow(struct hecmw_varray_int *varray)
{
  size_t new_max_val;

  HECMW_assert(varray);

  if (varray->max_val == 0)
    new_max_val = VARRAY_MAX_VAL_INIT;
  else
    new_max_val = varray->max_val * VARRAY_MAX_VAL_GROW;

  return varray_resize(varray, new_max_val);
}

int HECMW_varray_int_append(struct hecmw_varray_int *varray, int value)
{
  HECMW_assert(varray);

  if (varray->n_val == varray->max_val)
    if (varray_grow(varray) != HECMW_SUCCESS)
      return HECMW_ERROR;

  varray->vals[varray->n_val] = value;
  varray->n_val++;

  return HECMW_SUCCESS;
}

int HECMW_varray_int_get(const struct hecmw_varray_int *varray, size_t index)
{
  HECMW_assert(varray);
  HECMW_assert(0 <= index && index < varray->n_val);

  return varray->vals[index];
}

int HECMW_varray_int_cat(struct hecmw_varray_int *varray,
                         const struct hecmw_varray_int *varray2)
{
  size_t i;

  HECMW_assert(varray);
  HECMW_assert(varray2);

  while (varray->n_val + varray2->n_val > varray->max_val) {
    if (varray_grow(varray) != HECMW_SUCCESS)
      return HECMW_ERROR;
  }
  for (i = 0; i < varray2->n_val; i++) {
    varray->vals[varray->n_val] = varray2->vals[i];
    varray->n_val++;
  }
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

void HECMW_varray_int_sort(struct hecmw_varray_int *varray)
{
  HECMW_assert(varray);
  qsort(varray->vals, varray->n_val, sizeof(int), int_cmp);
}

int HECMW_varray_int_search(struct hecmw_varray_int *varray, int value, size_t *index)
{
  int *p;
  p = bsearch(&value, varray->vals, varray->n_val, sizeof(int), int_cmp);
  if (p == NULL) return HECMW_ERROR;
  *index = p - varray->vals;
  return HECMW_SUCCESS;
}

size_t HECMW_varray_int_uniq(struct hecmw_varray_int *varray)
{
  size_t i, n_dup = 0;

  HECMW_assert(varray);

  for (i = 1; i < varray->n_val; i++) {
    if (varray->vals[i-1] == varray->vals[i]) {
      n_dup++;
    } else {
      if (n_dup > 0) {
        varray->vals[i-n_dup] = varray->vals[i];
      }
    }
  }

  varray->n_val -= n_dup;

  if (varray->n_val * 2 < varray->max_val)
    varray_resize(varray, varray->n_val); /* reduce memory usage */

  return n_dup;
}

int HECMW_varray_int_resize(struct hecmw_varray_int *varray, size_t len)
{
  HECMW_assert(varray);

  if (varray->max_val < len) {
    if (varray_resize(varray, len) != HECMW_SUCCESS)
      return HECMW_ERROR;
  }
  varray->n_val = len;
  return HECMW_SUCCESS;
}

int *HECMW_varray_int_get_v(struct hecmw_varray_int *varray)
{
  HECMW_assert(varray);
  return varray->vals;
}

const int *HECMW_varray_int_get_cv(const struct hecmw_varray_int *varray)
{
  HECMW_assert(varray);
  return varray->vals;
}

int HECMW_varray_int_copy(const struct hecmw_varray_int *varray,
                          struct hecmw_varray_int *varray2)
{
  size_t i;

  HECMW_assert(varray);

  if (HECMW_varray_int_resize(varray2, varray->n_val) != HECMW_SUCCESS)
    return HECMW_ERROR;

  for (i = 0; i < varray->n_val; i++)
    varray2->vals[i] = varray->vals[i];

  return HECMW_SUCCESS;
}

int HECMW_varray_int_rmdup(struct hecmw_varray_int *varray)
{
  struct hecmw_varray_int tmp_array;

  HECMW_assert(varray);

  if (HECMW_varray_int_init(&tmp_array) != HECMW_SUCCESS) {
    return HECMW_ERROR;
  }
  if (HECMW_varray_int_copy(varray, &tmp_array) != HECMW_SUCCESS) {
    return HECMW_ERROR;
  }
  HECMW_varray_int_sort(&tmp_array);

  if (HECMW_varray_int_uniq(&tmp_array) != 0) {
    struct hecmw_bit_array ba;
    size_t i, n_dup = 0;

    HECMW_bit_array_init(&ba, tmp_array.n_val);
    for (i = 0; i < varray->n_val; i++) {
      int *key = varray->vals + i;
      int *res = bsearch(key, tmp_array.vals, tmp_array.n_val, sizeof(int), int_cmp);
      size_t idx;

      HECMW_assert(res != NULL);

      idx = res - tmp_array.vals;

      if (HECMW_bit_array_get(&ba, idx)) {
        n_dup++;
      } else {
        HECMW_bit_array_set(&ba, idx);
        varray->vals[i-n_dup] = varray->vals[i];
      }
    }
    varray->n_val -= n_dup;
    HECMW_bit_array_finalize(&ba);

    HECMW_assert(varray->n_val == tmp_array.n_val);
  }
  HECMW_varray_int_finalize(&tmp_array);
  return HECMW_SUCCESS;
}

int HECMW_varray_int_assign(struct hecmw_varray_int *varray,
                            size_t begin, size_t end, int val)
{
  size_t i;

  HECMW_assert(varray);
  HECMW_assert(0 <= begin);
  HECMW_assert(end <= varray->n_val);

  for (i = begin; i < end; i++) {
    varray->vals[i] = val;
  }
  return HECMW_SUCCESS;
}

int HECMW_varray_int_insert(struct hecmw_varray_int *varray,
                            size_t index, int val)
{
  size_t i;

  HECMW_assert(varray);
  HECMW_assert(0 <= index && index <= varray->n_val);

  if (varray->n_val == varray->max_val)
    if (varray_grow(varray) != HECMW_SUCCESS)
      return HECMW_ERROR;

  /* for (i = varray->n_val; i > index; i--) { */
  /*   varray->vals[i] = varray->vals[i-1]; */
  /* } */
  memmove(varray->vals + index + 1, varray->vals + index,
          sizeof(int) * (varray->n_val - index));

  varray->vals[index] = val;
  varray->n_val++;

  return HECMW_SUCCESS;
}

int HECMW_varray_int_delete(struct hecmw_varray_int *varray,
                            size_t index)
{
  size_t i;

  HECMW_assert(varray);
  HECMW_assert(0 <= index && index <= varray->n_val);

  /* for (i = index+1; i < varray->n_val; i++) { */
  /*   varray->vals[i-1] = varray->vals[i]; */
  /* } */
  memmove(varray->vals + index, varray->vals + index + 1,
          sizeof(int) * (varray->n_val - index - 1));

  varray->n_val--;

  return HECMW_SUCCESS;
}
