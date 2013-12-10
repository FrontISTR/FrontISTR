#ifndef HECMW_VIS_PVR_MAIN_H_INCLUDED
#define HECMW_VIS_PVR_MAIN_H_INCLUDED

#include "hecmw_struct.h"
#include "hecmw_result.h"
#include "hecmw_util.h"
#include "hecmw_vis_ray_trace.h"

void HECMW_vis_pvr_rendering(struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data, int *timestep, int *init_flag,
		int num_of_pvr, Parameter_vr *vr,
		int stat_para[NUM_CONTROL_PVR], char *outfile,  HECMW_Comm VIS_COMM);

#endif /* HECMW_VIS_PVR_MAIN_H_INCLUDED */







