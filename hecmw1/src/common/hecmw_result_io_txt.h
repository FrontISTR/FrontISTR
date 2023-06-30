/*****************************************************************************
 * Copyright (c) 2023 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_RESULT_IO_TXT_INCLUDED
#define HECMW_RESULT_IO_TXT_INCLUDED

#include "hecmw_result.h"

extern int HECMW_result_io_txt_write_by_fname(char *filename);
extern int HECMW_result_io_txt_write_ST_by_fname(
    char *filename, struct hecmwST_result_data *result, int n_node, int n_elem,
    char *header, char *comment);
extern struct hecmwST_result_data *HECMW_result_io_txt_read_by_fname(char *filename);

#endif
