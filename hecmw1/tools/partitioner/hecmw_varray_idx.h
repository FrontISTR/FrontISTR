/*****************************************************************************
 * Copyright (c) 2026 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_VARRAY_IDX_INCLUDED
#define HECMW_VARRAY_IDX_INCLUDED

#include <stddef.h>

#ifdef HECMW_PART_WITH_METIS
#include "metis.h"
#ifndef METIS_VER_MAJOR
/* METIS 4 names its index type idxtype and has no idx_t */
typedef idxtype idx_t;
#endif
#else
typedef long long idx_t;
#endif

struct hecmw_varray_idx {
  size_t n_val;
  size_t max_val;

  idx_t *vals;
};

extern int HECMW_varray_idx_init(struct hecmw_varray_idx *varray);

extern void HECMW_varray_idx_finalize(struct hecmw_varray_idx *varray);

extern size_t HECMW_varray_idx_nval(const struct hecmw_varray_idx *varray);

extern int HECMW_varray_idx_append(struct hecmw_varray_idx *varray, idx_t value);

extern idx_t HECMW_varray_idx_get(const struct hecmw_varray_idx *varray,
                                  size_t index);

extern int HECMW_varray_idx_cat(struct hecmw_varray_idx *varray,
                                const struct hecmw_varray_idx *varray2);

extern void HECMW_varray_idx_sort(struct hecmw_varray_idx *varray);

extern int HECMW_varray_idx_search(struct hecmw_varray_idx *varray, idx_t value,
                                   size_t *index);

extern size_t HECMW_varray_idx_uniq(struct hecmw_varray_idx *varray);

extern int HECMW_varray_idx_resize(struct hecmw_varray_idx *varray, size_t len);

extern idx_t *HECMW_varray_idx_get_v(struct hecmw_varray_idx *varray);

extern const idx_t *HECMW_varray_idx_get_cv(
    const struct hecmw_varray_idx *varray);

extern int HECMW_varray_idx_copy(const struct hecmw_varray_idx *varray,
                                 struct hecmw_varray_idx *varray2);

extern int HECMW_varray_idx_rmdup(struct hecmw_varray_idx *varray);

extern int HECMW_varray_idx_assign(struct hecmw_varray_idx *varray,
                                   size_t begin, size_t end, idx_t val);

extern int HECMW_varray_idx_insert(struct hecmw_varray_idx *varray,
                                   size_t index, idx_t val);

extern int HECMW_varray_idx_delete(struct hecmw_varray_idx *varray,
                                   size_t index);

#endif /* HECMW_VARRAY_IDX_INCLUDED */
