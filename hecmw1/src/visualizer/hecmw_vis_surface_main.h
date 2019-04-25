/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_VIS_SURFACE_MAIN_H_INCLUDED
#define HECMW_VIS_SURFACE_MAIN_H_INCLUDED

#include "hecmw_struct.h"
#include "hecmw_result.h"
#include "hecmw_util.h"
#include "hecmw_vis_SF_geom.h"
#include "hecmw_vis_psf_rendering.h"

void HECMW_vis_psf_rendering(struct hecmwST_local_mesh *mesh,
                             struct hecmwST_result_data *data, int *timestep,
                             struct surface_module *sf, Parameter_rendering *sr,
                             int stat_para[NUM_CONTROL_PSF], char *outfile1, char *body,
                             HECMW_Comm VIS_COMM);

#endif /* HECMW_VIS_SURFACE_MAIN_H_INCLUDED */
