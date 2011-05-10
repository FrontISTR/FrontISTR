/*
	fstr_rmerge Ver.1.3
	2010.09.27 by K.Suemitsu (AdvanceSoft)
	2004.11.10 by N.Imai (RIST)
	--------------------------------------------------
	FISTR result 出力ファイルを１つのファイルに統合する。
	FISTR 中で HECMW_write_result で作成されたランクごと
	のファイルを一つのファイルに統合する独立したツール
	計算に用いたコントロールファイル、メッシュデータが必要
	単一領域にも対応
*/

#include "fstr_rmerge_util.h"

FILE* log_fp;

void error_stop( void )
{
	char* msg;
	HECMW_get_error( &msg );
	fputs( msg, stderr );
	exit( -1 );
}

void set_fname( int argc, char** argv, char* out_fname, int* binary, int* refine )
{
	int i;
	char* fname = NULL;

	*binary = 0;
	*refine = 0;
	for( i=1; i<argc; i++ ) {
		if( strcmp(argv[i], "-h") == 0 ) {
			fprintf( stderr, "[usage] fstr_rmerge [-o binary/text -refine] out_file\n" );
			exit( -1 );
		} else if( strcmp(argv[i], "-o" ) == 0 ) {
			if( argc == i+1 ) {
				fprintf( stderr, "[usage] fstr_rmerge [-o binary/text -refine] out_file\n" );
				exit( -1 );
			}
			i++;
			if( strcmp(argv[i], "binary" ) == 0 ) {
				*binary = 1;
			} else if( strcmp(argv[i], "text" ) == 0 ) {
				*binary = 0;
			} else {
				fprintf( stderr, "Error : text or binary is required after -o\n" );
				exit( -1 );
			}
		} else if( strcmp(argv[i], "-refine" ) == 0 ) {
			*refine = 1;
		} else {
			fname = argv[i];
		}
	}
	if( !fname ) {
		fprintf( stderr, "[usage] fstr_rmerge [-o binary/text -refine] out_file\n" );
		exit( -1 );
	} else {
		strcpy( out_fname, fname );
	}
}

int main(int argc, char** argv )
{
	int area_n, step_n, binary, refine;
	int step, rcode;

	char out_fheader[ HECMW_FILENAME_LEN+1 ];
	char out_fname[ HECMW_FILENAME_LEN+1 ];
	struct hecmwST_local_mesh** mesh;
	struct hecmwST_local_mesh* glmesh;
	struct hecmwST_result_data* data;
	fstr_res_info** res;
	fstr_gl_t* glt;
	int n_node, n_elem;
	char header[ HECMW_HEADER_LEN+1 ];

	log_fp = stderr;

	if( HECMW_init( &argc, &argv )) { error_stop(); }

	set_fname( argc, argv, out_fheader, &binary, &refine );
	fprintf( log_fp, "out file name header is %s\n", out_fheader );

	mesh = fstr_get_all_local_mesh( "fstrMSH", &area_n );
	if( !mesh ) { error_stop(); }

	fprintf( log_fp, "table creating .. \n" );
	glt = fstr_create_glt( mesh, area_n );
	if( !glt ) {
		fprintf( stderr, "ERROR : Cannot create global_local table." );
		fstr_free_mesh( mesh, area_n );
		exit( -1 );
	}

	glmesh = fstr_create_glmesh( glt );
	if( !glmesh ) {
		fprintf( stderr, "ERROR : Cannot create global table." );
		fstr_free_mesh( mesh, area_n );
		fstr_free_gl_t( glt );
		exit( -1 );
	}

	step_n = fstr_get_step_n( "fstrRES" );

	for( step=1; step<=step_n; step++ ){
		fprintf( log_fp, "step:%d .. reading .. ", step );
		res = fstr_get_all_result( "fstrRES", step, area_n );
		if( !res ){
			fprintf( stderr, "ERROR : Cannot create result structure.\n" );
			fstr_free_mesh( mesh, area_n );
			fstr_free_gl_t( glt );
			fstr_free_glmesh( glmesh );
			exit( -1 );
		}
		fprintf( log_fp, "end\n" );

		if( refine ) {
			fprintf( log_fp, "refined element table creating .. \n" );
			glt = fstr_refine_glt( glt, res, area_n );
			if( !glt ) {
				fprintf( stderr, "ERROR : Cannot create refined element table." );
				fstr_free_mesh( mesh, area_n );
				fstr_free_gl_t( glt );
				fstr_free_glmesh( glmesh );
				fstr_free_result( res, area_n );
				exit( -1 );
			}
		}

		fprintf( log_fp, "step:%d .. combining .. ", step );
		data = fstr_all_result( glt, res, &n_node, &n_elem );
		if( !data ){
			fprintf( stderr, "ERROR : Cannot combine result structure.\n" );
			fstr_free_mesh( mesh, area_n );
			fstr_free_gl_t( glt );
			fstr_free_glmesh( glmesh );
			fstr_free_result( res, area_n );
			exit( -1 );
		}
		fprintf( log_fp, "end\n" );

		sprintf( out_fname, "%s.%d", out_fheader, step );
		fprintf( log_fp, "output to %s .. ", out_fname );
		HECMW_result_get_header( header );
		HECMW_result_init( glmesh, step, header );
		if( binary ) {
			rcode = HECMW_result_write_bin_ST_by_fname( out_fname, data, n_node, n_elem, header );
		} else {
			rcode = HECMW_result_write_txt_ST_by_fname( out_fname, data, n_node, n_elem, header );
		}
		if( rcode ) {
			fprintf( stderr, "ERROR : Cannot open/write file %s\n", out_fname );
			fstr_free_mesh( mesh, area_n );
			fstr_free_gl_t( glt );
			fstr_free_glmesh( glmesh );
			fstr_free_result( res, area_n );
			HECMW_result_free( data );
			exit( -1 );
		}
		fprintf( log_fp, "end\n" );

		fstr_free_result( res, area_n );
		HECMW_result_free( data );
	}

	fstr_free_mesh( mesh, area_n );
	fstr_free_gl_t( glt );
	fstr_free_glmesh( glmesh );

	HECMW_finalize();

	exit( 0 );
}

