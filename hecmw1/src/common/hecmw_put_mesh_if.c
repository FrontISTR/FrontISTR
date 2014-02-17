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
#include "hecmw_struct.h"
#include "hecmw_util.h"
#include "hecmw_dist_copy_f2c.h"
#include "hecmw_io_put_mesh.h"
#include "hecmw_dist_free.h"

static struct hecmwST_local_mesh *mesh;

static int
alloc_struct(void)
{
	mesh = HECMW_malloc(sizeof(*mesh));
	if(mesh == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}

	mesh->section = HECMW_malloc(sizeof(*mesh->section));
	if(mesh->section == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}

	mesh->material = HECMW_malloc(sizeof(*mesh->material));
	if(mesh->material == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}

	mesh->mpc = HECMW_malloc(sizeof(*mesh->mpc));
	if(mesh->mpc == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}

	mesh->amp = HECMW_malloc(sizeof(*mesh->amp));
	if(mesh->amp == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}

	mesh->node_group = HECMW_malloc(sizeof(*mesh->node_group));
	if(mesh->node_group == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}

	mesh->elem_group = HECMW_malloc(sizeof(*mesh->elem_group));
	if(mesh->elem_group == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}

	mesh->surf_group = HECMW_malloc(sizeof(*mesh->surf_group));
	if(mesh->surf_group == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}

	mesh->contact_pair = HECMW_malloc(sizeof(*mesh->contact_pair));
	if(mesh->contact_pair == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}

	return 0;
}


/*----------------------------------------------------------------------------*/


void
hecmw_put_mesh_if(char *name_ID, int *err, int len)
{
	char cname[HECMW_NAME_LEN+1];

	*err = 1;

	if(HECMW_strcpy_f2c_r(name_ID, len, cname, sizeof(cname)) == NULL) {
		return;
	}

	if(HECMW_put_mesh(mesh, cname)) {
		return;
	}

	*err = 0; 
}



void
hecmw_put_mesh_if_(char *name_ID, int *err, int len)
{
	hecmw_put_mesh_if(name_ID, err, len);
}



void
hecmw_put_mesh_if__(char *name_ID, int *err, int len)
{
	hecmw_put_mesh_if(name_ID, err, len);
}



void
HECMW_PUT_MESH_IF(char *name_ID, int *err, int len)
{
	hecmw_put_mesh_if(name_ID, err, len);
}


/*----------------------------------------------------------------------------*/


void
hecmw_put_mesh_init_if(int *err)
{
	*err = 1;

	if(alloc_struct()) {
		return;
	}

	if(HECMW_dist_copy_f2c_init(mesh)) {
		return;
	}

	*err = 0;
}



void
hecmw_put_mesh_init_if_(int *err)
{
	hecmw_put_mesh_init_if(err);
}



void
hecmw_put_mesh_init_if__(int *err)
{
	hecmw_put_mesh_init_if(err);
}



void
HECMW_PUT_MESH_INIT_IF(int *err)
{
	hecmw_put_mesh_init_if(err);
}


/*----------------------------------------------------------------------------*/


void
hecmw_put_mesh_finalize_if(int *err)
{
	*err = 1;

	if(HECMW_dist_copy_f2c_finalize()) {
		return;
	}
	HECMW_dist_free(mesh);
	mesh = NULL;

	*err = 0;
}



void
hecmw_put_mesh_finalize_if_(int *err)
{
	hecmw_put_mesh_finalize_if(err);
}



void
hecmw_put_mesh_finalize_if__(int *err)
{
	hecmw_put_mesh_finalize_if(err);
}



void
HECMW_PUT_MESH_FINALIZE_IF(int *err)
{
	hecmw_put_mesh_finalize_if(err);
}
