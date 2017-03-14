/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hecmw_struct.h"
#include "hecmw_io_get_mesh.h"
#include "hecmw_dist_free.h"

#include "hecmw_couple_struct.h"
#include "hecmw_couple_info.h"



extern struct hecmwST_local_mesh *
HECMW_couple_get_mesh(char *name_ID, char *unit_ID)
{
	struct hecmwST_local_mesh *mesh = NULL;
	struct hecmw_couple_comm *intracomm = NULL;

	if((mesh = HECMW_get_mesh(name_ID)) == NULL) return NULL;
	if((intracomm = HECMW_couple_get_intracomm_u(unit_ID)) == NULL) goto error;

	mesh->HECMW_COMM = intracomm->comm;
	mesh->PETOT      = intracomm->psize;
	mesh->my_rank    = intracomm->rank;

	HECMW_couple_free_comm(intracomm);

	return mesh;

error:
	HECMW_dist_free(mesh);
	return NULL;
}
