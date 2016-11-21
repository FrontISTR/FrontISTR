/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_COUPLE_JUDGE
#define INC_HECMW_COUPLE_JUDGE

#include "hecmw_struct.h"

extern int
HECMW_couple_judge_tet1(const struct hecmwST_local_mesh *local_mesh, int elem, int surf_id,
		double coord_px, double coord_py, double coord_pz, double *dot_product, double *distance);
extern int
HECMW_couple_judge_hex1(const struct hecmwST_local_mesh *local_mesh, int elem, int surf_id,
		double coord_px, double coord_py, double coord_pz, double *dot_product, double *distance);

#endif	/* INC_HECMW_COUPLE_JUDGE */
