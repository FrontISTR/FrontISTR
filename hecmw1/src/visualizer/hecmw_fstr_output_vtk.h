#ifndef HECMW_FSTR_OUTPUT_VKT_H_INCLUDED
#define HECMW_FSTR_OUTPUT_VKT_H_INCLUDED

#include "hecmw_struct.h"
#include "hecmw_result.h"
#include "hecmw_util.h"

void
HECMW_vkt_output (struct hecmwST_local_mesh *mesh,
		struct hecmwST_result_data *data, char *outfile, char *outfile1, int *max_timestep, HECMW_Comm VIS_COMM);

void
HECMW_bin_vkt_output (struct hecmwST_local_mesh *mesh,
		struct hecmwST_result_data *data, char *outfile, char *outfile1, int *max_timestep, HECMW_Comm VIS_COMM);

#endif /* HECMW_FSTR_OUTPUT_VKT_H_INCLUDED */
