/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuaki Sakane (RIST)                         *
 *                       Shin'ichi Ezure (RIST)                        *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "hecmw_struct.h"
#include "hecmw_util.h"
#include "hecmw_io_get_mesh.h"
#include "hecmw_io_mesh.h"
#include "hecmw_io_hec.h"
#include "hecmw_io_geofem.h"
#include "hecmw_io_abaqus.h"
#include "hecmw_io_dist.h"
#include "hecmw_dist_refine.h"


static struct hecmwST_local_mesh *
get_entire_mesh(struct hecmw_ctrl_meshfiles *files)
{
	int i;
	struct hecmwST_local_mesh *mesh;

	if(HECMW_io_init()) return NULL;
	HECMW_log(HECMW_LOG_DEBUG, "io_init done");

	if(HECMW_io_pre_process()) return NULL;
	HECMW_log(HECMW_LOG_DEBUG, "io_pre_process done");

	for(i=0; i < files->n_mesh; i++) {
		struct hecmw_ctrl_meshfile *file = &files->meshfiles[i];

		switch(file->type) {
			case HECMW_CTRL_FTYPE_HECMW_ENTIRE:
				if(HECMW_read_entire_mesh(file->filename)) return NULL;
				break;
			case HECMW_CTRL_FTYPE_GEOFEM:
				if(HECMW_read_geofem_mesh(file->filename)) return NULL;
				break;
			case HECMW_CTRL_FTYPE_ABAQUS:
				if(HECMW_read_abaqus_mesh(file->filename)) return NULL;
				break;
			default:
				HECMW_assert(0);
		}
	}
	HECMW_log(HECMW_LOG_DEBUG, "reading mesh done\n");

	if(HECMW_io_post_process()) return NULL;
	HECMW_log(HECMW_LOG_DEBUG, "post_process done\n");

	mesh = HECMW_io_make_local_mesh();
	if(mesh == NULL) return NULL;
	HECMW_log(HECMW_LOG_DEBUG, "converting mesh done\n");

	if(HECMW_io_finalize()) return NULL;
	HECMW_log(HECMW_LOG_DEBUG, "io_finalize done\n");

	return mesh;
}


struct hecmwST_local_mesh *
HECMW_get_mesh(char *name_ID)
{
	struct hecmw_ctrl_meshfiles *files;
	struct hecmwST_local_mesh *mesh;
	char filename[HECMW_FILENAME_LEN+1];
	char *cad_filename;
	FILE* fp;

	files = HECMW_ctrl_get_meshfiles(name_ID);
	if(files == NULL) return NULL;

	if(files->n_mesh == 1 && files->meshfiles[0].type == HECMW_CTRL_FTYPE_HECMW_DIST) {
		strcpy(filename, files->meshfiles[0].filename);
		mesh = HECMW_get_dist_mesh(filename);
	} else {
		mesh = get_entire_mesh(files);
	}

	strcpy(filename, files->meshfiles[0].filename);
	strtok(filename, ".");
	strcat(filename, ".rnf");
	if((fp = fopen(filename, "r")) == NULL) {
		cad_filename = NULL;
	} else {
		fclose(fp);
		cad_filename = filename;
	}

	if(HECMW_dist_refine(&mesh, files->meshfiles[0].refine, cad_filename, NULL) != HECMW_SUCCESS) {
		HECMW_dist_free(mesh);
		return NULL;
	}

	HECMW_ctrl_free_meshfiles(files);

	return mesh;
}
