/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_COUPLE_S2N_WITH_AREA
#define INC_HECMW_COUPLE_S2N_WITH_AREA

#include "hecmw_struct.h"
#include "hecmw_couple_weight.h"

extern struct hecmw_couple_weight_list *HECMW_couple_s2n_with_area(
    const struct hecmwST_local_mesh *mesh,
    const struct hecmw_couple_boundary *boundary);

#endif /* INC_HECMW_COUPLE_S2N_WITH_AREA */
