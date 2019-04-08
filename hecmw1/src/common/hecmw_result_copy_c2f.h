/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_RESULT_COPY_C2F_INCLUDED
#define HECMW_RESULT_COPY_C2F_INCLUDED

#include "hecmw_result.h"

extern int HECMW_result_copy_c2f_init(struct hecmwST_result_data *result_data,
                                      int n_node, int n_elem);
extern int HECMW_result_copy_c2f_finalize(void);

#endif
