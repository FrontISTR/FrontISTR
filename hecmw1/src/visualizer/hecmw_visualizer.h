/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_VISUALIZER_H_INCLUDED
#define HECMW_VISUALIZER_H_INCLUDED

#include "hecmw_struct.h"
#include "hecmw_result.h"
#include "hecmw_util.h"
#include "hecmw_vis_SF_geom.h"
#include "hecmw_vis_ray_trace.h"

extern PSF_link *psf;
extern PVR_link *pvr;

int HECMW_visualize_init(void);

int HECMW_visualize_init_by_comm(HECMW_Comm comm);

int HECMW_visualize(struct hecmwST_local_mesh *mesh,
                    struct hecmwST_result_data *result, int timestep);

int HECMW_visualize_finalize(void);

#endif /* HECMW_VISUALIZER_H_INCLUDED */
