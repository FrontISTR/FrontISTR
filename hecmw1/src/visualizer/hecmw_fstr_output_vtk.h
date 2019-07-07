/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_FSTR_OUTPUT_VKT_H_INCLUDED
#define HECMW_FSTR_OUTPUT_VKT_H_INCLUDED

#include "hecmw_struct.h"
#include "hecmw_result.h"
#include "hecmw_util.h"

void
HECMW_vtk_output (struct hecmwST_local_mesh *mesh,
		struct hecmwST_result_data *data, char *outfile, char *outfile1, HECMW_Comm VIS_COMM);

void
HECMW_bin_vtk_output (struct hecmwST_local_mesh *mesh,
		struct hecmwST_result_data *data, char *outfile, char *outfile1, HECMW_Comm VIS_COMM);

#endif /* HECMW_FSTR_OUTPUT_VKT_H_INCLUDED */
