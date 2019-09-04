/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_VIS_RENDERING_H_INCLUDED
#define HECMW_VIS_RENDERING_H_INCLUDED

#include "hecmw_struct.h"
#include "hecmw_result.h"
#include "hecmw_util.h"
#include "hecmw_vis_SF_geom.h"
#include "hecmw_vis_psf_rendering.h"

void HECMW_vis_rendering_surface(struct surface_module *sf,
                                 Parameter_rendering *sr,
                                 struct hecmwST_local_mesh *mesh,
                                 struct hecmwST_result_data *data, int tvertex,
                                 int tpatch, int *color_list, double *minvalue,
                                 double *maxvalue, Result *result,
                                 char *outfile, int stat_para[NUM_CONTROL_PSF],
                                 HECMW_Comm VIS_COMM, int timestep);

#endif /* HECMW_VIS_RENDERING_H_INCLUDED */
