/*****************************************************************************
 * Copyright (c) 2026 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "hecmw_util.h"
#include "hecmw_malloc.h"
#include "hecmw_config.h"
#include "hecmw_bit_array.h"
#include "hecmw_varray_idx.h"

enum { VARRAY_MAX_VAL_INIT = 64, VARRAY_MAX_VAL_GROW = 2 };

int HECMW_varray_idx_init(struct hecmw_varray_idx *varray) {
  HECMW_assert(varray);

  varray->n_val   = 0;
  varray->max_val = 0;

  varray->vals = NULL;

  return HECMW_SUCCESS;
}

void HECMW_varray_idx_finalize(struct hecmw_varray_idx *varray) {
  HECMW_assert(varray);

  if (varray->max_val == 0) {
    HECMW_assert(varray->n_val == 0);
    return;
  }

  HECMW_free(varray->vals);
  return;
}

size_t HECMW_varray_idx_nval(const struct hecmw_varray_idx *varray) {
  HECMW_assert(varray);

  return varray->n_val;
}

static int varray_resize(struct hecmw_varray_idx *varray, size_t new_max_val) {
  idx_t *new_vals;

  HECMW_assert(varray);
  HECMW_assert(varray->n_val <= new_max_val);

  if (varray->max_val == new_max_val) return HECMW_SUCCESS;

  if (new_max_val == 0) {
    HECMW_assert(varray->vals);

    HECMW_free(varray->vals);
    varray->vals    = NULL;
    varray->max_val = 0;

    return HECMW_SUCCESS;
  }

  new_vals = (idx_t *)HECMW_realloc(varray->vals, sizeof(idx_t) * new_max_val);
  if (new_vals == NULL) {
    return HECMW_ERROR;
  }

  varray->vals    = new_vals;
  varray->max_val = new_max_val;

  return HECMW_SUCCESS;
}

static int varray_grow(struct hecmw_varray_idx *varray) {
  size_t new_max_val;

  HECMW_assert(varray);

  if (varray->max_val == 0)
    new_max_val = VARRAY_MAX_VAL_INIT;
  else
    new_max_val = varray->max_val * VARRAY_MAX_VAL_GROW;

  return varray_resize(varray, new_max_val);
}

int HECMW_varray_idx_append(struct hecmw_varray_idx *varray, idx_t value) {
  HECMW_assert(varray);

  if (varray->n_val == varray->max_val)
    if (varray_grow(varray) != HECMW_SUCCESS) return HECMW_ERROR;

  varray->vals[varray->n_val] = value;
  varray->n_val++;

  return HECMW_SUCCESS;
}

idx_t HECMW_varray_idx_get(const struct hecmw_varray_idx *varray,
                           size_t index) {
  HECMW_assert(varray);
  HECMW_assert(0 <= index && index < varray->n_val);

  return varray->vals[index];
}

int HECMW_varray_idx_cat(struct hecmw_varray_idx *varray,
                         const struct hecmw_varray_idx *varray2) {
  size_t i;

  HECMW_assert(varray);
  HECMW_assert(varray2);

  while (varray->n_val + varray2->n_val > varray->max_val) {
    if (varray_grow(varray) != HECMW_SUCCESS) return HECMW_ERROR;
  }
  for (i = 0; i < varray2->n_val; i++) {
    varray->vals[varray->n_val] = varray2->vals[i];
    varray->n_val++;
  }
  return HECMW_SUCCESS;
}

static int idx_cmp(const void *v1, const void *v2) {
  const idx_t *i1, *i2;

  i1 = (const idx_t *)v1;
  i2 = (const idx_t *)v2;

  if (*i1 < *i2) return -1;
  if (*i1 > *i2) return 1;
  return 0;
}

void HECMW_varray_idx_sort(struct hecmw_varray_idx *varray) {
  HECMW_assert(varray);
  qsort(varray->vals, varray->n_val, sizeof(idx_t), idx_cmp);
}

int HECMW_varray_idx_search(struct hecmw_varray_idx *varray, idx_t value,
                            size_t *index) {
  idx_t *p;
  p = bsearch(&value, varray->vals, varray->n_val, sizeof(idx_t), idx_cmp);
  if (p == NULL) return HECMW_ERROR;
  *index = p - varray->vals;
  return HECMW_SUCCESS;
}

size_t HECMW_varray_idx_uniq(struct hecmw_varray_idx *varray) {
  size_t i, n_dup = 0;

  HECMW_assert(varray);

  for (i = 1; i < varray->n_val; i++) {
    if (varray->vals[i - 1] == varray->vals[i]) {
      n_dup++;
    } else {
      if (n_dup > 0) {
        varray->vals[i - n_dup] = varray->vals[i];
      }
    }
  }

  varray->n_val -= n_dup;

  if (varray->n_val * 2 < varray->max_val)
    varray_resize(varray, varray->n_val); /* reduce memory usage */

  return n_dup;
}

int HECMW_varray_idx_resize(struct hecmw_varray_idx *varray, size_t len) {
  HECMW_assert(varray);

  if (varray->max_val < len) {
    if (varray_resize(varray, len) != HECMW_SUCCESS) return HECMW_ERROR;
  }
  varray->n_val = len;
  return HECMW_SUCCESS;
}

idx_t *HECMW_varray_idx_get_v(struct hecmw_varray_idx *varray) {
  HECMW_assert(varray);
  return varray->vals;
}

const idx_t *HECMW_varray_idx_get_cv(const struct hecmw_varray_idx *varray) {
  HECMW_assert(varray);
  return varray->vals;
}

int HECMW_varray_idx_copy(const struct hecmw_varray_idx *varray,
                          struct hecmw_varray_idx *varray2) {
  size_t i;

  HECMW_assert(varray);

  if (HECMW_varray_idx_resize(varray2, varray->n_val) != HECMW_SUCCESS)
    return HECMW_ERROR;

  for (i = 0; i < varray->n_val; i++) varray2->vals[i] = varray->vals[i];

  return HECMW_SUCCESS;
}

int HECMW_varray_idx_rmdup(struct hecmw_varray_idx *varray) {
  struct hecmw_varray_idx tmp_array;

  HECMW_assert(varray);

  if (HECMW_varray_idx_init(&tmp_array) != HECMW_SUCCESS) {
    return HECMW_ERROR;
  }
  if (HECMW_varray_idx_copy(varray, &tmp_array) != HECMW_SUCCESS) {
    return HECMW_ERROR;
  }
  HECMW_varray_idx_sort(&tmp_array);

  if (HECMW_varray_idx_uniq(&tmp_array) != 0) {
    struct hecmw_bit_array ba;
    size_t i, n_dup = 0;

    HECMW_bit_array_init(&ba, tmp_array.n_val);
    for (i = 0; i < varray->n_val; i++) {
      idx_t *key = varray->vals + i;
      idx_t *res =
          bsearch(key, tmp_array.vals, tmp_array.n_val, sizeof(idx_t), idx_cmp);
      size_t idx;

      HECMW_assert(res != NULL);

      idx = res - tmp_array.vals;

      if (HECMW_bit_array_get(&ba, idx)) {
        n_dup++;
      } else {
        HECMW_bit_array_set(&ba, idx);
        varray->vals[i - n_dup] = varray->vals[i];
      }
    }
    varray->n_val -= n_dup;
    HECMW_bit_array_finalize(&ba);

    HECMW_assert(varray->n_val == tmp_array.n_val);
  }
  HECMW_varray_idx_finalize(&tmp_array);
  return HECMW_SUCCESS;
}

int HECMW_varray_idx_assign(struct hecmw_varray_idx *varray, size_t begin,
                            size_t end, idx_t val) {
  size_t i;

  HECMW_assert(varray);
  HECMW_assert(0 <= begin);
  HECMW_assert(end <= varray->n_val);

  for (i = begin; i < end; i++) {
    varray->vals[i] = val;
  }
  return HECMW_SUCCESS;
}

int HECMW_varray_idx_insert(struct hecmw_varray_idx *varray, size_t index,
                            idx_t val) {
  HECMW_assert(varray);
  HECMW_assert(0 <= index && index <= varray->n_val);

  if (varray->n_val == varray->max_val)
    if (varray_grow(varray) != HECMW_SUCCESS) return HECMW_ERROR;

  memmove(varray->vals + index + 1, varray->vals + index,
          sizeof(idx_t) * (varray->n_val - index));

  varray->vals[index] = val;
  varray->n_val++;

  return HECMW_SUCCESS;
}

int HECMW_varray_idx_delete(struct hecmw_varray_idx *varray, size_t index) {
  HECMW_assert(varray);
  HECMW_assert(0 <= index && index <= varray->n_val);

  memmove(varray->vals + index, varray->vals + index + 1,
          sizeof(idx_t) * (varray->n_val - index - 1));

  varray->n_val--;

  return HECMW_SUCCESS;
}
