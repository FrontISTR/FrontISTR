/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_COUPLE_INTRA_IFTABLE
#define INC_HECMW_COUPLE_INTRA_IFTABLE

#include "hecmw_struct.h"
#include "hecmw_couple_comm.h"
#include "hecmw_couple_boundary_info.h"

struct hecmw_couple_intra_iftable {
  int n_neighbor_pe;
  int *neighbor_pe;
  int *import_index;
  int *import_item;
  int *export_index;
  int *export_item;
};

extern void HECMW_couple_free_intra_iftable(
    struct hecmw_couple_intra_iftable *intra_tbl);

extern struct hecmw_couple_intra_iftable *HECMW_couple_alloc_intra_iftable(
    void);

extern struct hecmw_couple_intra_iftable *HECMW_couple_make_intra_iftable(
    const struct hecmwST_local_mesh *mesh,
    const struct hecmw_couple_boundary *boundary,
    const struct hecmw_couple_comm *intracomm);

#endif /* INC_HECMW_COUPLE_INTRA_IFTABLE */
