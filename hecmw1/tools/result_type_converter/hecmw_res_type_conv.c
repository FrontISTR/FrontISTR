/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.3                                               *
 *                                                                     *
 *     Last Update : 2010/08/26                                        *
 *        Category : HEC-MW Utility                                    *
 *                                                                     *
 *            Written by Keiji Suemitsu (AdvanceSoft)                  *
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
#include "hecmw_result.h"
#include "hecmw_util.h"
#include "hecmw_io.h"

int stepid = 1;
char out_file[HECMW_NAME_LEN+1] = "";


void help(void)
{
	printf(" HECMW Result File Type Converter\n");
	printf("usage)  res_conv [options]\n");
	printf("[option]\n");
	printf(" -h             : help\n");
	printf(" -o [file]      : output file name without rank and step number\n");
	printf(" -s [step]      : step number (default:%d)\n", stepid);
}


int set_params( int argc, char** argv )
{
	int i;

	for( i = 1; i < argc; i++ ) {
		if( strcmp(argv[i], "-h" ) == 0 ) {
			help();
			return -1;
		} else if( strcmp(argv[i], "-o" ) == 0 ) {
			if( argc == i+1 ) {
				fprintf( stderr, "Error : paramter required after %s\n", argv[i] );
				return -1;
			}
			i++;
			strcpy( out_file, argv[i] );
		} else if( strcmp(argv[i], "-s" ) == 0 ) {
			if( argc == i+1 ) {
				fprintf( stderr, "Error : paramter required after %s\n", argv[i] );
				return -1;
			}
			i++;
			if(sscanf( argv[i], "%d", &stepid ) != 1 ) {
				fprintf( stderr, "Error : parameter %s cannot be converted to step number\n", argv[i] );
				return -1;
			}
		} else {
			fprintf(stderr, "Error : invalid parameter %s\n", argv[i] );
			help();
			return -1;
		}
	}

	return 0;
}


int
main(int argc , char **argv)
{
	struct hecmwST_result_data *data;
	char *fileheader, buf[HECMW_FILENAME_LEN], resultfile[HECMW_FILENAME_LEN];
	char header[HECMW_HEADER_LEN];
	int n_node, n_elem, rcode;
	int mynode;

	if( HECMW_init( &argc, &argv ) ) {
		HECMW_abort( HECMW_comm_get_comm() );
	}

	if( set_params( argc, argv ) )  {
		HECMW_abort( HECMW_comm_get_comm() );
	}

	mynode = HECMW_comm_get_rank();

	fileheader = HECMW_ctrl_get_result_fileheader( "result", buf, HECMW_FILENAME_LEN );
	sprintf( resultfile, "%s.%d.%d", fileheader, mynode, stepid );
	fprintf(stdout,"Input file : %s\n", resultfile);
	data = HECMW_result_read_by_fname( resultfile );
	if(!data ) {
		HECMW_abort( HECMW_comm_get_comm() );
	}

	if ( out_file[0] ) {
		sprintf( resultfile, "%s.%d.%d", out_file, mynode, stepid );
	}
	fprintf(stdout,"Output file : %s\n", resultfile);
	n_node = HECMW_result_get_nnode();
	n_elem = HECMW_result_get_nelem();
	HECMW_result_get_header( header );
	rcode = HECMW_result_write_txt_ST_by_fname( resultfile, data, n_node, n_elem, header );
	if( rcode ) {
		HECMW_abort( HECMW_comm_get_comm() );
	}

	HECMW_result_free(data);

	HECMW_finalize();

	return 0;
}

