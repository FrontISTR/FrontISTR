/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *
 * Exodus II output via NetCDF API (without libexodus dependency)
 *****************************************************************************/
#ifndef HECMW_FSTR_OUTPUT_EXO_H_INCLUDED
#define HECMW_FSTR_OUTPUT_EXO_H_INCLUDED

#include "hecmw_struct.h"
#include "hecmw_result.h"

/* output_type=18 (EXODUS): All timesteps in one file */
void HECMW_exodus_output(struct hecmwST_local_mesh *mesh,
                         struct hecmwST_result_data *data,
                         char *outfile, char *outfile1,
                         HECMW_Comm VIS_COMM);

/* output_type=19 (STEP_EXODUS): One file per timestep */
void HECMW_exodus_step_output(struct hecmwST_local_mesh *mesh,
                              struct hecmwST_result_data *data,
                              char *outfile, char *outfile1,
                              HECMW_Comm VIS_COMM);

#endif /* HECMW_FSTR_OUTPUT_EXO_H_INCLUDED */
