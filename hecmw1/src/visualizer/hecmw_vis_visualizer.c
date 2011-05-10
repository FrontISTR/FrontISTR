/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.1                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : Visualization                                     *
 *                                                                     *
 *            Written by Li Chen (Univ. of Tokyo)                      *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/


#include "hecmw_vis_SF_geom.h"
#include "hecmw_vis_ray_trace_c.h"
extern void HECMW_vis_read_control(FILE *fp, int pesize, int mynode);
extern void HECMW_vis_psf_rendering(struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data, int *timestep, int *init_flag,
		int num_of_psf, struct surface_module *sf, Parameter_rendering *sr, int stat_para[NUM_CONTROL_PSF],
		char *outfile1,  HECMW_Comm VIS_COMM);
extern void HECMW_vis_pvr_rendering(struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data, int *timestep, int *init_flag,
		int num_of_pvr, Parameter_vr *vr,
		int stat_para[NUM_CONTROL_PVR], char *outfile,  HECMW_Comm VIS_COMM);


PSF_link *psf;
PVR_link *pvr;

int
main(int argc , char **argv)
{
	int i, j, ii;
	char *name_ID = NULL;
	struct hecmwST_local_mesh *mesh;
	struct hecmwST_result_data *data;
	char  *contfile, *outfile, *resultfile, buf[HECMW_FILENAME_LEN], buf1[HECMW_FILENAME_LEN], buf2[HECMW_FILENAME_LEN],
	resultfile1[HECMW_FILENAME_LEN], outfile1[HECMW_FILENAME_LEN];
	FILE  *contfp;
	int min_step, max_step;
	PSF_link *tp1, *tp2;
	PVR_link *tv1, *tv2;
	int  visual_id, init_flag, timestep;
	Parameter_rendering   *sr;
	struct surface_module	*sf;
	Parameter_vr   *vr;
	int  stat_para_sf[NUM_CONTROL_PSF], stat_para_vr[NUM_CONTROL_PVR];
	HECMW_Comm  VIS_COMM;
	int pesize, mynode;
	int flag_read;

	if(HECMW_init(&argc, &argv)) {
		abort();
	}

	if(argc == 2) {
		name_ID = argv[1];
	}
	mesh = HECMW_get_mesh("test");
	if(mesh == NULL) {
		HECMW_abort(HECMW_comm_get_comm());
	}
	HECMW_DEBUG(("hecmw_get_mesh OK"));
	HECMW_Comm_dup(mesh->HECMW_COMM, &VIS_COMM);
	HECMW_Comm_size(VIS_COMM, &pesize);
	HECMW_Comm_rank(VIS_COMM, &mynode);

	outfile=HECMW_ctrl_get_result_fileheader("vis_out", buf1, HECMW_FILENAME_LEN);

	if((contfp=fopen("hecmw_vis.ini","r"))== NULL){
		contfile=HECMW_ctrl_get_control_file("vis", buf, HECMW_FILENAME_LEN);
		if((contfp=fopen(contfile,"r"))== NULL)
			HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0011: Cannot find control file");
	}
	HECMW_vis_read_control(contfp, pesize, mynode);
	fclose(contfp);
	min_step=100000000;
	max_step=-10000000;
	tp1=psf->next_psf;
	tv1=pvr->next_pvr;
	for(i=0;i<psf->num_of_psf;i++) {
		if((tp1->visual_start_step!=-1) && (tp1->visual_start_step<min_step))
			min_step=tp1->visual_start_step;
		if((tp1->visual_end_step!=-1) && (tp1->visual_end_step>max_step))
			max_step=tp1->visual_end_step;
	}
	for(i=0;i<pvr->num_of_pvr;i++) {
		if((tv1->visual_start_step!=-1) && (tv1->visual_start_step<min_step))
			min_step=tv1->visual_start_step;
		if((tv1->visual_end_step!=-1) && (tv1->visual_end_step>max_step))
			max_step=tv1->visual_end_step;
	}
	if((min_step==100000000) && (max_step==-10000000)) {
		resultfile=HECMW_ctrl_get_result_fileheader("result", buf2, HECMW_FILENAME_LEN);
		sprintf(resultfile1, "%s.%d.1", resultfile, mynode);
		/*		if(mynode==0)
        fprintf(stderr, "the result file is %s\n", resultfile1);
		 */
		data=HECMW_result_read_by_fname(resultfile1);
		if(data == NULL) {
			HECMW_abort(HECMW_comm_get_comm());
		}
		HECMW_DEBUG(("hecmw_get_result OK"));
		/*	exchange(mesh, data);

		for(j=0;j<mesh->n_node;j++) {
			if(mesh->global_node_ID[j]==1001) {
          fprintf(stderr, "NODEID: 1001\n");
        for(i=0;i<16;i++)
          fprintf(stderr, "%lf  ", data->node_val_item[16*j+i]);
          fprintf(stderr, "\n");
			}
			if(mesh->global_node_ID[j]==1003) {
          fprintf(stderr, "NODEID: 1003\n");
        for(i=0;i<16;i++)
          fprintf(stderr, "%lf  ", data->node_val_item[16*j+i]);
          fprintf(stderr, "\n");
			}
			if(mesh->global_node_ID[j]==1021) {
          fprintf(stderr, "NODEID: 1021\n");
        for(i=0;i<16;i++)
          fprintf(stderr, "%lf  ", data->node_val_item[16*j+i]);
          fprintf(stderr, "\n");
			}
			if(mesh->global_node_ID[j]==1221) {
          fprintf(stderr, "NODEID: 1221\n");
        for(i=0;i<16;i++)
          fprintf(stderr, "%lf  ", data->node_val_item[16*j+i]);
          fprintf(stderr, "\n");
			}
		}
		 */

		if(psf->num_of_psf>0) {
			timestep=0;
			init_flag=1;
			tp1=psf->next_psf;
			for(visual_id=0;visual_id<psf->num_of_psf;visual_id++) {
				if(mynode==0)
					fprintf(stderr, "Start visualize PSF %d\n", visual_id+1);
				sf=tp1->sf;
				sr=tp1->sr;
				for(ii=0;ii<NUM_CONTROL_PSF;ii++)
					stat_para_sf[ii]=tp1->stat_para[ii];
				tp1=tp1->next_psf;
				if(psf->num_of_psf>1)
					sprintf(outfile1, "%s_psf%d", outfile, visual_id+1);
				else
					sprintf(outfile1, "%s_psf", outfile);
				HECMW_vis_psf_rendering(mesh, data, &timestep, &init_flag, psf->num_of_psf, sf, sr, stat_para_sf, outfile1,  VIS_COMM);
				init_flag=0;
			}
		}

		if(pvr->num_of_pvr>0) {
			tv1=pvr->next_pvr;
			init_flag=1;
			timestep=0;
			for(visual_id=0;visual_id<pvr->num_of_pvr;visual_id++) {
				if(mynode==0)
					fprintf(stderr, "Start visualize PVR %d\n", visual_id+1);
				vr=tv1->vr;
				for(ii=0;ii<NUM_CONTROL_PVR;ii++)
					stat_para_vr[ii]=tv1->stat_para[ii];
				tv1=tv1->next_pvr;
				if(pvr->num_of_pvr>1)
					sprintf(outfile1, "%s_pvr%d", outfile, visual_id+1);
				else
					sprintf(outfile1, "%s_pvr", outfile);

				HECMW_vis_pvr_rendering(mesh, data, &timestep, &init_flag, pvr->num_of_pvr, vr, stat_para_vr, outfile1,  VIS_COMM);
				init_flag=0;
			}
		}
		HECMW_result_free(data);
	}
	else {
		for(timestep=min_step;timestep<=max_step;timestep++) {
			flag_read=1;
			if(psf->num_of_psf>0) {
				init_flag=1;
				tp1=psf->next_psf;
				for(visual_id=0;visual_id<psf->num_of_psf;visual_id++) {
					if(((timestep-tp1->visual_start_step) % tp1->visual_interval_step) == 0) {
						if(mynode==0)
							fprintf(stderr, "Start visualize PSF %d at timestep %d\n", visual_id+1, timestep);
						if(flag_read==1) {
							resultfile=HECMW_ctrl_get_result_fileheader("result", buf2, HECMW_FILENAME_LEN);
							sprintf(resultfile1, "%s.%d.%d", resultfile, mynode, timestep);
							data=HECMW_result_read_by_fname(resultfile1);
							if(data == NULL) {
								HECMW_abort(HECMW_comm_get_comm());
							}
							flag_read=0;
						}

						sf=tp1->sf;
						sr=tp1->sr;
						for(ii=0;ii<NUM_CONTROL_PSF;ii++)
							stat_para_sf[ii]=tp1->stat_para[ii];
						tp1=tp1->next_psf;
						if(psf->num_of_psf>1)
							sprintf(outfile1, "%s_psf%d.%d", outfile, visual_id+1, timestep);
						else
							sprintf(outfile1, "%s_psf.%d", outfile, timestep);
						HECMW_vis_psf_rendering(mesh, data, &timestep, &init_flag, psf->num_of_psf, sf, sr, stat_para_sf, outfile1,  VIS_COMM);
						init_flag=0;
					}
				}
			}
			if(pvr->num_of_pvr>0) {
				tv1=pvr->next_pvr;
				init_flag=1;
				for(visual_id=0;visual_id<pvr->num_of_pvr;visual_id++) {
					if(((timestep-tv1->visual_start_step) % tv1->visual_interval_step) == 0) {
						if(mynode==0)
							fprintf(stderr, "Start visualize PVR %d at timestep %d\n", visual_id+1, timestep);
						if(flag_read==1) {
							resultfile=HECMW_ctrl_get_result_fileheader("result", buf2, HECMW_FILENAME_LEN);
							sprintf(resultfile1, "%s.%d.%d", resultfile, mynode, timestep);
							data=HECMW_result_read_by_fname(resultfile1);
							if(data == NULL) {
								HECMW_abort(HECMW_comm_get_comm());
							}
							flag_read=0;
						}
						vr=tv1->vr;
						for(ii=0;ii<NUM_CONTROL_PVR;ii++)
							stat_para_vr[ii]=tv1->stat_para[ii];
						tv1=tv1->next_pvr;
						if(pvr->num_of_pvr>1)
							sprintf(outfile1, "%s_pvr%d", outfile, visual_id+1);
						else
							sprintf(outfile1, "%s_pvr", outfile);

						HECMW_vis_pvr_rendering(mesh, data, &timestep, &init_flag, pvr->num_of_pvr, vr, stat_para_vr, outfile1,  VIS_COMM);

						init_flag=0;
					}
				}
			}
			if(flag_read==0)
				HECMW_result_free(data);
		}
	}





	HECMW_dist_free(mesh);

	HECMW_finalize();

	return 0;
}
/*
void  exchange(	struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data)
{
	int i, j, k;
	int tn_component;
	double *result1;
	int tmp_count=0;

	tn_component=0;
	for(i=0;i<data->nn_component;i++)
		tn_component+=data->nn_dof[i];
	result1=(double *)calloc(mesh->n_node*tn_component, sizeof(double));
	if(result1==NULL)
		HECMW_vis_print_exit("Not enough memory for result exchange");
	tmp_count=0;
	for(j=0;j<data->nn_component;j++) {
		for(i=0;i<mesh->n_node;i++) {
			for(k=0;k<data->nn_dof[j];k++)
		  result1[tn_component*i+tmp_count+k]=data->node_val_item[tmp_count*mesh->n_node+i*data->nn_dof[j]+k];
		}
		tmp_count+=data->nn_dof[j];
	}
	free(data->node_val_item);
	data->node_val_item=result1;
	return;
}

 */

