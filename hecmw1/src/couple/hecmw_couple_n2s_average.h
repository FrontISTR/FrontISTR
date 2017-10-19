/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_COUPLE_N2S_AVERAGE
#define INC_HECMW_COUPLE_N2S_AVERAGE

#include "hecmw_struct.h"
#include "hecmw_couple_weight.h"
#include "hecmw_couple_struct.h"
#include "hecmw_couple_intra_iftable.h"

extern struct hecmw_couple_weight_list *HECMW_couple_n2s_average(
    const struct hecmwST_local_mesh *mesh,
    const struct hecmw_couple_boundary *boundary,
    const struct hecmw_couple_comm *intracomm,
    const struct hecmw_couple_intra_iftable *intra_tbl);

#endif /* INC_HECMW_COUPLE_N2S_AVERAGE */
