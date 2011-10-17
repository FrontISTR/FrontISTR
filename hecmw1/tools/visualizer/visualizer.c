/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.3                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : HEC-MW Utility                                    *
 *                                                                     *
 *            Written by Li Chen (Univ. of Tokyo)                      *
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
#include "hecmw_result.h"
#include "hecmw_util.h"
#include "hecmw_io.h"
#include "hecmw_visualizer.h"

extern PSF_link *psf;
extern PVR_link *pvr;


int
main(int argc , char **argv)
{
	int i;
	char *name_ID = NULL;
	struct hecmwST_local_mesh *mesh;
	struct hecmwST_result_data *data;
	char  *resultfile, buf2[HECMW_FILENAME_LEN], resultfile1[HECMW_FILENAME_LEN];
	int min_step, max_step;
	PSF_link *tp1;
	PVR_link *tv1;
	int  timestep;
	int mynode;

	if(HECMW_init(&argc, &argv)) {
		abort();
	}

	if(argc == 2) {
		name_ID = argv[1];
	}
	mesh = HECMW_get_mesh("mesh");
	if(mesh == NULL) {
		HECMW_abort(HECMW_comm_get_comm());
	}
	HECMW_DEBUG(("hecmw_get_mesh OK"));

	HECMW_Comm_rank(mesh->HECMW_COMM, &mynode);

	HECMW_visualize_init();

	min_step=100000000;
	max_step=-10000000;
	tp1=psf->next_psf;
	tv1=pvr->next_pvr;
	for(i=0;i<psf->num_of_psf;i++) {
		if((tp1->visual_start_step!=-1) && (tp1->visual_start_step<min_step))
			min_step=tp1->visual_start_step;
		if((tp1->visual_end_step!=-1) && (tp1->visual_end_step>max_step))
			max_step=tp1->visual_end_step;
		tp1=tp1->next_psf;
	}
	for(i=0;i<pvr->num_of_pvr;i++) {
		if((tv1->visual_start_step!=-1) && (tv1->visual_start_step<min_step))
			min_step=tv1->visual_start_step;
		if((tv1->visual_end_step!=-1) && (tv1->visual_end_step>max_step))
			max_step=tv1->visual_end_step;
		tv1=tv1->next_pvr;
	}
	if((min_step==100000000) && (max_step==-10000000)) {
		min_step = 1;
		max_step = 1;
	}

	for(timestep=min_step;timestep<=max_step;timestep++) {
		int read_flag = 0;
		int visual_id;
		if(psf->num_of_psf>0) {
			tp1=psf->next_psf;
			for(visual_id=0;visual_id<psf->num_of_psf;visual_id++) {
				if((((timestep-tp1->visual_start_step) % tp1->visual_interval_step) == 0) || (timestep==max_step)) {
					read_flag = 1;
					break;
				}
				tp1=tp1->next_psf;
			}
		}
		if(read_flag == 0 && pvr->num_of_pvr>0) {
			tv1=pvr->next_pvr;
			for(visual_id=0;visual_id<pvr->num_of_pvr;visual_id++) {
				if((((timestep-tv1->visual_start_step) % tv1->visual_interval_step) == 0) || (timestep==max_step)) {
					read_flag = 1;
					break;
				}
			}
		}
		if (read_flag == 0) continue;

		resultfile=HECMW_ctrl_get_result_fileheader("result", buf2, HECMW_FILENAME_LEN);
		sprintf(resultfile1, "%s.%d.%d", resultfile, mynode, timestep);
		data=HECMW_result_read_by_fname(resultfile1);
		if(data == NULL) {
			HECMW_abort(HECMW_comm_get_comm());
		}

		HECMW_visualize(mesh, data, timestep, max_step, 0);

		HECMW_result_free(data);
	}

	HECMW_visualize_finalize();

	HECMW_dist_free(mesh);

	HECMW_finalize();

	return 0;
}
