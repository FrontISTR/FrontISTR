/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_COUPLE_S2N_AVERAGE
#define INC_HECMW_COUPLE_S2N_AVERAGE

#include "hecmw_struct.h"
#include "hecmw_couple_weight.h"

extern struct hecmw_couple_weight_list *HECMW_couple_s2n_average(
    const struct hecmwST_local_mesh *mesh,
    const struct hecmw_couple_boundary *boundary);

#endif /* INC_HECMW_COUPLE_S2N_AVERAGE */
