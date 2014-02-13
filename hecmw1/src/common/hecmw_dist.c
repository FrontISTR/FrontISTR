/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuaki Sakane (RIST)                         *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/



#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "hecmw_dist.h"
#include "hecmw_util.h"
#include "hecmw_struct.h"


int
HECMW_dist_get_mat_id(const struct hecmwST_material *mat, const char *name)
{
	int i;

	if(mat == NULL) return -1;
	if(name == NULL) return -1;

	for(i=0; i < mat->n_mat; i++) {
		if(strcmp(mat->mat_name[i], name) == 0) {
			return i+1;
		}
	}
	return -1;
}


int
HECMW_dist_get_ngrp_id(const struct hecmwST_node_grp *ngrp, const char *name)
{
	int i;

	if(ngrp == NULL) return -1;
	if(name == NULL) return -1;

	for(i=0; i < ngrp->n_grp; i++) {
		if(strcmp(ngrp->grp_name[i], name) == 0) {
			return i+1;
		}
	}
	return -1;
}


int
HECMW_dist_get_egrp_id(const struct hecmwST_elem_grp *egrp, const char *name)
{
	int i;

	if(egrp == NULL) return -1;
	if(name == NULL) return -1;

	for(i=0; i < egrp->n_grp; i++) {
		if(strcmp(egrp->grp_name[i], name) == 0) {
			return i+1;
		}
	}
	return -1;
}


int
HECMW_dist_get_sgrp_id(const struct hecmwST_surf_grp *sgrp, const char *name)
{
	int i;

	if(sgrp == NULL) return -1;
	if(name == NULL) return -1;

	for(i=0; i < sgrp->n_grp; i++) {
		if(strcmp(sgrp->grp_name[i], name) == 0) {
			return i+1;
		}
	}
	return -1;
}


int
HECMW_dist_gid2lid_node(const struct hecmwST_local_mesh *mesh, int gid)
{
	int i;

	if(mesh == NULL) return -1;

	for(i=0; i < mesh->n_node; i++) {
		if(mesh->global_node_ID[i] == gid) return i+1;
	}
	return -1;
}


int
HECMW_dist_gid2lid_elem(const struct hecmwST_local_mesh *mesh, int gid)
{
	int i;

	if(mesh == NULL) return -1;

	for(i=0; i < mesh->n_elem; i++) {
		if(mesh->global_elem_ID[i] == gid) return i+1;
	}
	return -1;
}


