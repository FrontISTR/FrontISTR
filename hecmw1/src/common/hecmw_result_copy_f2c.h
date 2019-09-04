/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_RESULT_COPY_F2C
#define INC_HECMW_RESULT_COPY_F2C

#include "hecmw_result.h"

extern int HECMW_result_copy_f2c_init(struct hecmwST_result_data *result_data,
                                      int n_node, int n_elem);
extern int HECMW_result_copy_f2c_finalize(void);

#endif
