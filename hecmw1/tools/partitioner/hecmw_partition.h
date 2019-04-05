/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_PARTITION
#define INC_HECMW_PARTITION

#include "hecmw_struct.h"
#include "hecmw_part_struct.h"

extern struct hecmwST_local_mesh *HECMW_partition_inner(
    struct hecmwST_local_mesh *global_mesh,
    struct hecmw_part_cont_data *cont_data);

extern struct hecmwST_local_mesh *HECMW_partition(
    struct hecmwST_local_mesh *local_mesh);

#endif /* INC_HPWMC_PARTITION */
