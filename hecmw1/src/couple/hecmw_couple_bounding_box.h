/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_COUPLE_BOUNDING_BOX
#define INC_HECMW_COUPLE_BOUNDING_BOX

#include "hecmw_struct.h"
#include "hecmw_couple_boundary_info.h"

struct hecmw_couple_box {
  double min_x;
  double min_y;
  double min_z;
  double max_x;
  double max_y;
  double max_z;
};

struct hecmw_couple_bounding_box {
  double tolerance;
  double coef;
  struct hecmw_couple_box *just;
  struct hecmw_couple_box *enlarged;
};

extern void HECMW_couple_free_bounding_box(
    struct hecmw_couple_bounding_box *bounding_box);

extern struct hecmw_couple_bounding_box *HECMW_couple_set_bounding_box(
    const char *boundary_id, const struct hecmwST_local_mesh *mesh,
    const struct hecmw_couple_boundary *boundary);

#endif /* INC_HECMW_COUPLE_BOUNDING_BOX */
