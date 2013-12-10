/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : Coupling Interface                                *
 *                                                                     *
 *            Written by Shin'ichi Ezure (RIST)                        *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using Hight End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "hecmw_struct.h"
#include "hecmw_dist_copy_f2c.h"
#include "hecmw_couple_define.h"
#include "hecmw_couple_init.h"

static struct hecmwST_local_mesh *mesh_unit1 = NULL, *mesh_unit2 = NULL;

/*================================================================================================*/

static struct hecmwST_local_mesh *
alloc_struct_local_mesh(void)
{
	struct hecmwST_local_mesh *mesh = NULL;

	mesh = (struct hecmwST_local_mesh *)HECMW_malloc(sizeof(struct hecmwST_local_mesh));
	if(mesh == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}
	mesh->section = NULL;
	mesh->material = NULL;
	mesh->mpc = NULL;
	mesh->amp = NULL;
	mesh->node_group = NULL;
	mesh->elem_group = NULL;
	mesh->surf_group = NULL;

	mesh->section = (struct hecmwST_section *)HECMW_malloc(sizeof(struct hecmwST_section));
	if(mesh->section == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}

	mesh->material = (struct hecmwST_material *)HECMW_malloc(sizeof(struct hecmwST_material));
	if(mesh->material == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}

	mesh->mpc = (struct hecmwST_mpc *)HECMW_malloc(sizeof(struct hecmwST_mpc));
	if(mesh->mpc == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}

	mesh->amp = (struct hecmwST_amplitude *)HECMW_malloc(sizeof(struct hecmwST_amplitude));
	if(mesh->amp == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}

	mesh->node_group = (struct hecmwST_node_grp *)HECMW_malloc(sizeof(struct hecmwST_node_grp));
	if(mesh->node_group == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}

	mesh->elem_group = (struct hecmwST_elem_grp *)HECMW_malloc(sizeof(struct hecmwST_elem_grp));
	if(mesh->elem_group == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}

	mesh->surf_group = (struct hecmwST_surf_grp *)HECMW_malloc(sizeof(struct hecmwST_surf_grp));
	if(mesh->surf_group == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}

	return mesh;

error:
	HECMW_dist_free(mesh);
	return NULL;
}


/*================================================================================================*/

extern void
hecmw_couple_init_if(char *boundary_id, int *err, int len)
{
	char c_id[HECMW_NAME_LEN+1];

	*err = 1;

	if(HECMW_strcpy_f2c_r(boundary_id, len, c_id, sizeof(c_id)) == NULL) return;
	if(HECMW_couple_init(c_id, mesh_unit1, mesh_unit2) != HECMW_SUCCESS) return;

	*err = 0;
}



extern void
hecmw_couple_init_if_(char *boundary_id, int *err, int len)
{
	hecmw_couple_init_if(boundary_id, err, len);
}



extern void
hecmw_couple_init_if__(char *boundary_id, int *err, int len)
{
	hecmw_couple_init_if(boundary_id, err, len);
}



extern void
HECMW_COUPLE_INIT_IF(char *boundary_id, int *err, int len)
{
	hecmw_couple_init_if(boundary_id, err, len);
}


/*------------------------------------------------------------------------------------------------*/

extern void
hecmw_couple_init_init_if(int *unit_specifier, int *err)
{
	*err = 1;

	if(*unit_specifier == HECMW_COUPLE_UNIT1) {
		if((mesh_unit1 = alloc_struct_local_mesh()) == NULL) return;
		if(HECMW_dist_copy_f2c_init(mesh_unit1)) return;
	} else if(*unit_specifier == HECMW_COUPLE_UNIT2) {
		if((mesh_unit2 = alloc_struct_local_mesh()) == NULL) return;
		if(HECMW_dist_copy_f2c_init(mesh_unit2)) return;
	} else {
		return;
	}

	*err = 0;
}



extern void
hecmw_couple_init_init_if_(int *unit_specifier, int *err)
{
	hecmw_couple_init_init_if(unit_specifier, err);
}



extern void
hecmw_couple_init_init_if__(int *unit_specifier, int *err)
{
	hecmw_couple_init_init_if(unit_specifier, err);
}



extern void
HECMW_COUPLE_INIT_INIT_IF(int *unit_specifier, int *err)
{
	hecmw_couple_init_init_if(unit_specifier, err);
}


/*------------------------------------------------------------------------------------------------*/

extern void
hecmw_couple_init_final_if(int *err)
{
	*err = 1;

	if(HECMW_dist_copy_f2c_finalize()) return;
	HECMW_dist_free(mesh_unit1);
	HECMW_dist_free(mesh_unit2);

	mesh_unit1 = NULL;
	mesh_unit2 = NULL;

	*err = 0;
}



extern void
hecmw_couple_init_final_if_(int *err)
{
	hecmw_couple_init_final_if(err);
}



extern void
hecmw_couple_init_final_if__(int *err)
{
	hecmw_couple_init_final_if(err);
}



extern void
HECMW_COUPLE_INIT_FINAL_IF(int *err)
{
	hecmw_couple_init_final_if(err);
}
