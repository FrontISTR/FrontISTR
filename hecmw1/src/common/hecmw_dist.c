/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "hecmw_dist.h"
#include "hecmw_util.h"
#include "hecmw_struct.h"


int
HECMW_dist_get_mat_id(const struct hecmwST_material *mat, const char *name)
{
	static int i;

	if(mat == NULL) return -1;
	if(name == NULL) return -1;

	if(i < mat->n_mat) {
		if(strcmp(mat->mat_name[i], name) == 0) {
			i++;
			return i;
		}
	}
	for(i=0; i < mat->n_mat; i++) {
		if(strcmp(mat->mat_name[i], name) == 0) {
			i++;
			return i;
		}
	}
	return -1;
}


int
HECMW_dist_get_ngrp_id(const struct hecmwST_node_grp *ngrp, const char *name)
{
	static int i;

	if(ngrp == NULL) return -1;
	if(name == NULL) return -1;

	if(i < ngrp->n_grp) {
		if(strcmp(ngrp->grp_name[i], name) == 0) {
			i++;
			return i;
		}
	}
	for(i=0; i < ngrp->n_grp; i++) {
		if(strcmp(ngrp->grp_name[i], name) == 0) {
			i++;
			return i;
		}
	}
	return -1;
}


int
HECMW_dist_get_egrp_id(const struct hecmwST_elem_grp *egrp, const char *name)
{
	static int i=0;

	if(egrp == NULL) return -1;
	if(name == NULL) return -1;

	if(i < egrp->n_grp) {
		if(strcmp(egrp->grp_name[i], name) == 0) {
			i++;
			return i;
		}
	}
	for(i=0; i < egrp->n_grp; i++) {
		if(strcmp(egrp->grp_name[i], name) == 0) {
			i++;
			return i;
		}
	}
	return -1;
}


int
HECMW_dist_get_sgrp_id(const struct hecmwST_surf_grp *sgrp, const char *name)
{
	static int i;

	if(sgrp == NULL) return -1;
	if(name == NULL) return -1;

	if(i < sgrp->n_grp) {
		if(strcmp(sgrp->grp_name[i], name) == 0) {
			i++;
			return i;
		}
	}
	for(i=0; i < sgrp->n_grp; i++) {
		if(strcmp(sgrp->grp_name[i], name) == 0) {
			i++;
			return i;
		}
	}
	return -1;
}


int
HECMW_dist_gid2lid_node(const struct hecmwST_local_mesh *mesh, int gid)
{
	static int i;

	if(mesh == NULL) return -1;

	if(i < mesh->n_node) {
		if(mesh->global_node_ID[i] == gid) {
			i++;
			return i;
		}
	}
	for(i=0; i < mesh->n_node; i++) {
		if(mesh->global_node_ID[i] == gid) {
			i++;
			return i;
		}
	}
	return -1;
}


int
HECMW_dist_gid2lid_elem(const struct hecmwST_local_mesh *mesh, int gid)
{
	static int i;

	if(mesh == NULL) return -1;

	if(i < mesh->n_elem) {
		if(mesh->global_elem_ID[i] == gid) {
			i++;
			return i;
		}
	}
	for(i=0; i < mesh->n_elem; i++) {
		if(mesh->global_elem_ID[i] == gid) {
			i++;
			return i;
		}
	}
	return -1;
}


