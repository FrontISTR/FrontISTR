/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "hecmw_struct.h"
#include "hecmw_util.h"
#include "hecmw_io_put_mesh.h"
#include "hecmw_io_dist.h"


int
HECMW_put_mesh(struct hecmwST_local_mesh *mesh, char *name_ID)
{
	struct hecmw_ctrl_meshfile *file;
	struct hecmw_ctrl_meshfiles *files;
	char filename[HECMW_FILENAME_LEN+1];

	files = HECMW_ctrl_get_meshfiles(name_ID);
	if(files == NULL) return -1;

	file = &files->meshfiles[0];

	strcpy(filename, files->meshfiles[0].filename);
	if(HECMW_put_dist_mesh(mesh, filename)) return -1;

	HECMW_ctrl_free_meshfiles(files);

	return 0;
}
