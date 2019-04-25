/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_DIST_REFINE_INCLUDED
#define HECMW_DIST_REFINE_INCLUDED

#include "hecmw_struct.h"

extern int HECMW_dist_refine(struct hecmwST_local_mesh **mesh, int refine,
                             const char *cad_filename,
                             const char *part_filename);

#endif
