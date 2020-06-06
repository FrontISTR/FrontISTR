/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_RESULT_IO_BIN_INCLUDED
#define HECMW_RESULT_IO_BIN_INCLUDED

#include "hecmw_result.h"

extern int HECMW_judge_result_bin_file(char *filename);
extern int HECMW_result_write_bin_by_fname(char *filename);
extern int HECMW_result_write_bin_ST_by_fname(
    char *filename, struct hecmwST_result_data *result, int n_node, int n_elem,
    char *header, char *comment);
extern struct hecmwST_result_data *HECMW_result_read_bin_by_fname(char *filename);

#endif
