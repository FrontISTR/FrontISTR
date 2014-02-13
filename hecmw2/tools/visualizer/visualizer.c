/*=====================================================================*
!                                                                      !
! Software Name : HEC-MW Ver 4.3                                      !
!                                                                      !
!      Module Name : Visualizer Main                                   !
!                                                                      !
!            Written by Keiji Suemitsu (Advancesoft)                   !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
*=====================================================================*/

#include "hecmw_visualizer.h"
#include "MW3_visualizer.h"
#include "API_C.h"

extern PSF_link *psf;
extern PVR_link *pvr;

int
main( int argc , char **argv )
{
	char controlfile[] = "./hecmw_ctrl.dat";
	int mynode, level, partID;
	struct hecmwST_local_mesh* mesh;
	struct hecmwST_result_data* data;
	char* header;
	char buf[HECMW_FILENAME_LEN], resultfile[HECMW_FILENAME_LEN];
	char ctrlfile[HECMW_FILENAME_LEN];
	int min_step, max_step, timestep, i, fg_text;
	PSF_link *tp1;
	PVR_link *tv1;

	mw_initialize_fstr_( &argc, argv, controlfile );

	mynode = mw_get_rank_();

	mw_file_read_fstr_();
	level = 0;
	partID = -1;
	mesh = MW3_get_mesh( level, partID );

	if( argc > 1 ) {
		strcpy( ctrlfile, argv[1] );
		if( MW3_visualize_init_ex( ctrlfile ) < 0 ) {
			HECMW_abort( HECMW_comm_get_comm() );
		}
	} else {
		if( MW3_visualize_init() < 0 ) {
			HECMW_abort( HECMW_comm_get_comm() );
		}
	}

/* MW2 world begin */

	min_step = 100000000;
	max_step =-100000000;
	tp1 = psf->next_psf;
	tv1 = pvr->next_pvr;
	for( i = 0; i < psf->num_of_psf; i++ ) {
		if( (tp1->visual_start_step!=-1) && (tp1->visual_start_step<min_step) )
			min_step = tp1->visual_start_step;
		if( (tp1->visual_end_step!=-1) && (tp1->visual_end_step>max_step) )
			max_step = tp1->visual_end_step;
		tp1 = tp1->next_psf;
	}
	for( i = 0; i < pvr->num_of_pvr; i++) {
		if( (tv1->visual_start_step!=-1) && (tv1->visual_start_step<min_step) )
			min_step = tv1->visual_start_step;
		if ((tv1->visual_end_step!=-1) && (tv1->visual_end_step>max_step) )
			max_step = tv1->visual_end_step;
		tv1 = tv1->next_pvr;
	}
	if( (min_step==100000000) && (max_step==-100000000) ) {
		min_step = 1;
		max_step = 1;
	}

	for( timestep=min_step; timestep <= max_step; timestep++ ) {
		int read_flag = 0;
		int visual_id;
		if( psf->num_of_psf > 0 ) {
			tp1 = psf->next_psf;
			for( visual_id = 0; visual_id < psf->num_of_psf; visual_id++ ) {
				if( ( ((timestep-tp1->visual_start_step)%tp1->visual_interval_step) == 0 ) ||
					( timestep == max_step ) ) {
					read_flag = 1;
					break;
				}
				tp1 = tp1->next_psf;
			}
		}
		if( read_flag == 0 && pvr->num_of_pvr > 0 ) {
			tv1 = pvr->next_pvr;
			for( visual_id = 0; visual_id < pvr->num_of_pvr; visual_id++ ) {
				if( ( ((timestep-tv1->visual_start_step)%tv1->visual_interval_step) == 0 ) ||
					( timestep == max_step ) ) {
					read_flag = 1;
					break;
				}
			}
		}
		if( read_flag == 0 ) continue;

		header = HECMW_ctrl_get_result_fileheader( "fstrRES", max_step, timestep, &fg_text );
		sprintf( resultfile, "%s.%d.%d", header, mynode, timestep );
		data = HECMW_result_read_by_fname( resultfile );
		if( data == NULL ) HECMW_abort( HECMW_comm_get_comm() );

		HECMW_visualize( mesh, data, timestep, max_step, 0 );

		HECMW_result_free( data );
	}

/* MW2 world end */

	MW3_visualize_finalize();

	MW3_free_mesh( mesh );

	mw_finalize_refine_();
	mw_finalize_();

	return 0;
}
