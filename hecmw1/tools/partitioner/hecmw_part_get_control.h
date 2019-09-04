/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_PART_GET_CONTROL
#define INC_HECMW_PART_GET_CONTROL

#include "hecmw_part_struct.h"

extern int HECMW_part_set_ctrl_file_name(char *fname);

extern struct hecmw_part_cont_data *HECMW_part_get_control();

extern void HECMW_part_free_control(struct hecmw_part_cont_data *cont_data);

#endif /* INC_HECMW_PART_GET_CONTROL */
