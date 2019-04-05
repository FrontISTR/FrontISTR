/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_FSTR_OUTPUT_FEMAP_H_INCLUDED
#define HECMW_FSTR_OUTPUT_FEMAP_H_INCLUDED

#include "hecmw_struct.h"
#include "hecmw_result.h"
#include "hecmw_util.h"

void HECMW_fstr_output_femap(struct hecmwST_local_mesh *mesh,
                             struct hecmwST_result_data *data, char *outfile,
                             HECMW_Comm VIS_COMM);
void HECMW_avs_output(struct hecmwST_local_mesh *mesh,
                      struct hecmwST_result_data *data, char *outfile,
                      HECMW_Comm VIS_COMM);
void HECMW_reorder_avs_output(struct hecmwST_local_mesh *mesh,
                              struct hecmwST_result_data *data, char *outfile,
                              HECMW_Comm VIS_COMM);
void HECMW_microavs_output(struct hecmwST_local_mesh *mesh,
                           struct hecmwST_result_data *data, char *outfile,
                           HECMW_Comm VIS_COMM);
void HECMW_bin_avs_output(struct hecmwST_local_mesh *mesh,
                          struct hecmwST_result_data *data, char *outfile,
                          HECMW_Comm VIS_COMM);

void HECMW_separate_avs_output(struct hecmwST_local_mesh *mesh,
                               struct hecmwST_result_data *data, char *outfile);

#endif /* HECMW_FSTR_OUTPUT_FEMAP_H_INCLUDED */
