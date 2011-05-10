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

#include "hecmw_visualizer.h"

#include <stdio.h>
#include <stdlib.h>
#include "hecmw_vis_mem_util.h"
#include "hecmw_vis_read_control.h"
#include "hecmw_vis_surface_main.h"
#include "hecmw_vis_pvr_main.h"

PSF_link *psf;
PVR_link *pvr;


int
HECMW_visualize_init(void)
{
	return HECMW_visualize_init_by_comm(HECMW_comm_get_comm());
}

int
HECMW_visualize_init_by_comm(HECMW_Comm VIS_COMM)
{
	FILE   *contfp;
	int    pesize, mynode;
	char  *contfile, buf[HECMW_FILENAME_LEN];

	HECMW_Comm_size(VIS_COMM, &pesize);
	HECMW_Comm_rank(VIS_COMM, &mynode);

	if((contfp=fopen("hecmw_vis.ini","r"))== NULL){
		contfile=HECMW_ctrl_get_control_file("vis_ctrl", buf, HECMW_FILENAME_LEN);
		if((contfp=fopen(contfile,"r"))== NULL)
			HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0011: Cannot open control file");
	}

	psf=(PSF_link *)malloc(sizeof(PSF_link));
	if(psf==NULL)
		HECMW_vis_memory_exit("psf");
	psf->next_psf=NULL;
	psf->num_of_psf=0;
	pvr=(PVR_link *)malloc(sizeof(PVR_link));
	if(pvr==NULL)
		HECMW_vis_memory_exit("pvr");
	pvr->next_pvr=NULL;
	pvr->num_of_pvr=0;

	HECMW_vis_read_control(contfp, pesize, mynode, psf, pvr);
	fclose(contfp);

	return 0;
}

int
HECMW_visualize( struct hecmwST_local_mesh *mesh,
		struct hecmwST_result_data *data,
		int timestep, int max_timestep, int is_force )
{
	int ii;
	char   *outfile, buf1[HECMW_FILENAME_LEN], outfile1[HECMW_FILENAME_LEN];
	PSF_link *tp1;
	PVR_link *tv1;
	int  visual_id, init_flag;
	Parameter_rendering   *sr;
	struct surface_module	*sf;
	Parameter_vr   *vr;
	int  stat_para_sf[NUM_CONTROL_PSF], stat_para_vr[NUM_CONTROL_PVR];
	HECMW_Comm  VIS_COMM;
	int pesize, mynode;
	/*
  int nn_item, ne_item;

  fprintf( stderr, "step=%d, max_step is_force=%d\n", timestep, max_timestep, is_force );

  fprintf( stderr, "nn_component=%d\n", result->nn_component );
  for( i=0; i<result->nn_component; i++ ) {
    fprintf( stderr, "nn_dof[%d]=%d\n", i, result->nn_dof[i] );
	}
  for( i=0; i<result->nn_component; i++ ) {
    fprintf( stderr, "node_label[%d]='%s'\n", i, result->node_label[i] );
  }
  for( nn_item=0, i=0; i<result->nn_component; i++ ) {
    nn_item += result->nn_dof[i];
  }
  for( i=0; i<mesh->n_node; i++ ) {
    for( j=0; j<nn_item; j++ ) {
      fprintf( stderr, " %f", result->node_val_item[nn_item*i+j] );
    }
    fprintf( stderr, "\n" );
  }

  fprintf( stderr, "ne_component=%d\n", result->ne_component );
  for( i=0; i<result->ne_component; i++ ) {
    fprintf( stderr, "ne_dof[%d]=%d\n", i, result->ne_dof[i] );
	}
  for( i=0; i<result->ne_component; i++ ) {
    fprintf( stderr, "elem_label[%d]=%s\n", i, result->elem_label[i] );
  }
  for( ne_item=0, i=0; i<result->ne_component; i++ ) {
    ne_item += result->ne_dof[i];
  }
  for( i=0; i<mesh->n_elem; i++ ) {
    for( j=0; j<ne_item; j++ ) {
      fprintf( stderr, " %f", result->elem_val_item[ne_item*i+j] );
    }
    fprintf( stderr, "\n" );
  }
	 */
	HECMW_Comm_dup(mesh->HECMW_COMM, &VIS_COMM);
	HECMW_Comm_size(VIS_COMM, &pesize);
	HECMW_Comm_rank(VIS_COMM, &mynode);
	outfile=HECMW_ctrl_get_result_fileheader("vis_out", buf1, HECMW_FILENAME_LEN);
	if(psf->num_of_psf>0) {
		init_flag=1;
		tp1=psf->next_psf;
		for(visual_id=0;visual_id<psf->num_of_psf;visual_id++) {
			if((((timestep-tp1->visual_start_step) % tp1->visual_interval_step) == 0) || (is_force==1) || (timestep ==max_timestep)){
				if(mynode==0)
					fprintf(stderr, "Start visualize PSF %d at timestep %d\n", visual_id+1, timestep);
				sf=tp1->sf;
				sr=tp1->sr;
				for(ii=0;ii<NUM_CONTROL_PSF;ii++)
					stat_para_sf[ii]=tp1->stat_para[ii];
				tp1=tp1->next_psf;
				if(psf->num_of_psf>1) {
					if(timestep>=1000)
						sprintf(outfile1, "%s_psf%d.%d", outfile, visual_id+1, timestep);
					else if((timestep>=100) && (timestep<=999))
						sprintf(outfile1, "%s_psf%d.0%d", outfile, visual_id+1, timestep);
					else if((timestep>=10) && (timestep<=99))
						sprintf(outfile1, "%s_psf%d.00%d", outfile, visual_id+1, timestep);
					else if(timestep<=9)
						sprintf(outfile1, "%s_psf%d.000%d", outfile, visual_id+1, timestep);
				}
				else {
					if(timestep>=1000)
						sprintf(outfile1, "%s_psf.%d", outfile,timestep);
					else if((timestep>=100) && (timestep<=999))
						sprintf(outfile1, "%s_psf.0%d", outfile, timestep);
					else if((timestep>=10) && (timestep<=99))
						sprintf(outfile1, "%s_psf.00%d", outfile, timestep);
					else if(timestep<=9)
						sprintf(outfile1, "%s_psf.000%d", outfile, timestep);
				}
				/*#ifdef DEBUG
             fprintf(stderr, "start call PSF\n");
             fprintf(stderr, "surface_num=%d surface_style=%d\n", sf[0].surface_style, sf[1].surface_style);
#endif
				 */
				HECMW_vis_psf_rendering(mesh, data, &timestep, sf, sr, stat_para_sf, outfile1,  VIS_COMM);
				init_flag=0;
			}
		}
	}
	if(pvr->num_of_pvr>0) {
		tv1=pvr->next_pvr;
		init_flag=1;
		for(visual_id=0;visual_id<pvr->num_of_pvr;visual_id++) {
			if((((timestep-tv1->visual_start_step) % tv1->visual_interval_step) == 0) || (is_force==1) || (timestep==max_timestep)) {
				if(mynode==0)
					fprintf(stderr, "Start visualize PVR %d at timestep %d\n", visual_id+1, timestep);
				vr=tv1->vr;
				for(ii=0;ii<NUM_CONTROL_PVR;ii++)
					stat_para_vr[ii]=tv1->stat_para[ii];
				tv1=tv1->next_pvr;
				if(pvr->num_of_pvr>1) {
					if(timestep>=1000)
						sprintf(outfile1, "%s_pvr%d.%d", outfile, visual_id+1, timestep);
					else if((timestep>=100) && (timestep<=999))
						sprintf(outfile1, "%s_pvr%d.0%d", outfile, visual_id+1, timestep);
					else if((timestep>=10) && (timestep<=99))
						sprintf(outfile1, "%s_pvr%d.00%d", outfile, visual_id+1, timestep);
					else if(timestep<=9)
						sprintf(outfile1, "%s_pvr%d.000%d", outfile, visual_id+1, timestep);
				}
				else {
					if(timestep>=1000)
						sprintf(outfile1, "%s_pvr.%d", outfile,timestep);
					else if((timestep>=100) && (timestep<=999))
						sprintf(outfile1, "%s_pvr.0%d", outfile, timestep);
					else if((timestep>=10) && (timestep<=99))
						sprintf(outfile1, "%s_pvr.00%d", outfile, timestep);
					else if(timestep<=9)
						sprintf(outfile1, "%s_pvr.000%d", outfile, timestep);
				}/*#ifdef DEBUG
            fprintf(stderr, "start call PVR\n");
            fprintf(stderr, "image_resolution=%d  %d  color_compo=%s\n", vr->xr, vr->yr, vr->color_comp_name);
#endif
				 */
				HECMW_vis_pvr_rendering(mesh, data, &timestep, &init_flag, pvr->num_of_pvr, vr, stat_para_vr, outfile1,  VIS_COMM);
				init_flag=0;
			}
		}
	}

	return 0;
}

extern int
HECMW_visualize_finalize( void )
{
	PSF_link *tp1, *tp2;
	PVR_link *tv1, *tv2;
	int i;
	if(psf->num_of_psf>0) {
		tp1=psf->next_psf;
		for(i=0;i<psf->num_of_psf;i++) {
			tp2=tp1;
			tp1=tp1->next_psf;
			free(tp2);
		}
	}
	free(psf);
	if(pvr->num_of_pvr>0) {
		tv1=pvr->next_pvr;
		for(i=0;i<pvr->num_of_pvr;i++) {
			tv2=tv1;
			tv1=tv1->next_pvr;
			free(tv2);
		}
	}
	free(pvr);


	return 0;
}
