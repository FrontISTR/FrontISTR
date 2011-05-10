/*
	fstr_res_util Ver.1.3
	2010.09.27 by K.Suemitsu (AdvanceSoft)
	2004.10.26 by N.Imai (RIST)
	--------------------------------------------------------
	分散で計算された結果を読込み処理するためのユーティリティ
*/

#include "fstr_rmerge_util.h"

FILE* log_fp;
int fg_log = 1;

static
void out_log( const char* fmt, ... )
{
	va_list arg;
	if( fg_log ){
		va_start( arg, fmt );
		vfprintf( log_fp, fmt, arg );
		va_end( arg );
	}
}


/*****************************************************************************/
/* 全分散メッシュの読込み                                                    */
/*****************************************************************************/

/* 分散の場合：fg_single==0, fname ファイル名のヘッダ */
/* それ以外　：fg_single==1, fname にはなにもセットされない */
/* 備考）メッシュファイル名にランク番号は付加されない */
static
int get_dist_fname( char* name_ID, char* fname, int* fg_single )
{
	struct hecmw_ctrl_meshfiles *files;

	files = HECMW_ctrl_get_meshfiles_header( name_ID );
	if(files == NULL) return -1;

	if( files->n_mesh == 1 ) {
		strcpy( fname, files->meshfiles[0].filename );
		if( files->meshfiles[0].type == HECMW_CTRL_FTYPE_HECMW_DIST ) {
			*fg_single = 0;
		} else {
			*fg_single = 1;
		}
	} else {
		HECMW_ctrl_free_meshfiles( files );
		return -1;
	}

	HECMW_ctrl_free_meshfiles( files );
	return 0;
}

/* 領域数をファイル数から得る */
static
int get_area_n( char* fname )
{
	FILE* fp;
	char buff[ HECMW_FILENAME_LEN+1 ];
	int area = 0;

	while(1) {
		sprintf( buff, "%s.%d", fname, area );
		out_log( "try open : %s  ... ", buff );
		fp = fopen( buff, "r" );
		if( !fp ) {
			out_log( "fail\n" );
			out_log( "area number is %d\n", area );
			return area;
		} else {
			out_log( "success\n" );
			fclose(fp);
		}
		area++;
	}
}

/* 全分散メッシュの読込み( area_n : 分散領域の数) */
struct hecmwST_local_mesh** fstr_get_all_local_mesh( char* name_ID, int* area_number )
{
	int i;
	int area_n;
	int fg_single;
	char fname[ HECMW_FILENAME_LEN+1 ];
	struct hecmwST_local_mesh** mesh;

	if( get_dist_fname( name_ID, fname, &fg_single )) return NULL;

	if( fg_single ) {
		out_log( "mesh file type is NOT HECMW_DIST.\n" );
		area_n = 1;
		out_log( "area number is %d\n", area_n );
		mesh = HECMW_malloc( area_n * sizeof( struct hecmwST_local_mesh* ) );
		mesh[0] = HECMW_get_mesh( name_ID );
		if( mesh[0] == NULL ) return NULL;
	} else {
		out_log( "mesh file type is HECMW_DIST.\n" );
		area_n = get_area_n( fname );
		if( area_n <= 0 ) return NULL;
		mesh = HECMW_malloc( area_n * sizeof( struct hecmwST_local_mesh* ) );
		for( i=0; i<area_n; i++ ) {
			char fname_a[256];
			sprintf( fname_a, "%s.%d", fname, i );
			out_log( "loading dist mesh from %s\n", fname_a );
			mesh[i] = HECMW_get_dist_mesh( fname_a );
			if( mesh[i] == NULL ) return NULL;
		}
	}
	*area_number = area_n;
	return mesh;
}


/*****************************************************************************/
/* メッシュの削除                                                            */
/*****************************************************************************/

void fstr_free_mesh( struct hecmwST_local_mesh** mesh, int area_n )
{
	int i;

	if( !mesh ) return;

	for( i=0; i<area_n; i++ ){
		HECMW_dist_free( mesh[i] );
	}
	HECMW_free( mesh );
	return;
}


/*****************************************************************************/
/* ステップ数を調べる（ファイルの存在を調べる）                              */
/*****************************************************************************/

int fstr_get_step_n( char* name_ID )
{
	FILE* fp;
	int step;
	char header[ HECMW_FILENAME_LEN+1 ];
	char fname[ HECMW_FILENAME_LEN+1 ];

	if( HECMW_ctrl_get_result_fileheader( name_ID, header, HECMW_FILENAME_LEN-10 ) == NULL )
		return 0;

	step = 1;
	while(1) {
		sprintf( fname, "%s.0.%d", header, step );
		out_log( "try open : %s  ... ", fname );
		fp = fopen( fname, "r" );
		if( !fp ) {
			out_log( "fail\n" );
			out_log( "step number is %d\n", step-1 );
			return step-1;
		} else {
			out_log( "success\n" );
			fclose(fp);
		}
		step++;
	}
}


/*****************************************************************************/
/* ステップ step の全領域のデータを得る                                      */
/*****************************************************************************/

fstr_res_info** fstr_get_all_result( char* name_ID, int step, int area_n )
{
	int i;
	char fname[ HECMW_FILENAME_LEN+1 ];
	fstr_res_info** res;

	if( HECMW_ctrl_get_result_fileheader( name_ID, fname, HECMW_FILENAME_LEN ) == NULL )
		return NULL;

	res = HECMW_malloc( area_n * sizeof(fstr_res_info*) );
	if( res == NULL ) return NULL;

	for( i=0; i<area_n; i++ ) {
		char fname_a[256];
		sprintf( fname_a, "%s.%d.%d", fname, i, step );
		res[i] = HECMW_malloc( sizeof(fstr_res_info) );
		if( res[i] == NULL ) return NULL;		
		res[i]->result = HECMW_result_read_by_fname( fname_a );
		if( res[i]->result == NULL ) return NULL;
		res[i]->nnode_gid = HECMW_result_get_nnode( );
		res[i]->node_gid = HECMW_malloc( res[i]->nnode_gid * sizeof(int) );
		if( res[i]->node_gid == NULL ) return NULL;
		HECMW_result_get_nodeID( res[i]->node_gid );
		res[i]->nelem_gid = HECMW_result_get_nelem( );
		res[i]->elem_gid = HECMW_malloc( res[i]->nelem_gid * sizeof(int) );
		if( res[i]->elem_gid == NULL ) return NULL;
		HECMW_result_get_elemID( res[i]->elem_gid );
	}
	return res;
}


/*****************************************************************************/
/* ステップ step の全領域のデータを結合する                                  */
/*****************************************************************************/

struct hecmwST_result_data* fstr_all_result( fstr_gl_t* glt,
	fstr_res_info** res, int* n_node, int* n_elem )
{
	int i, j, k, g_id, l_id, area;
	int nitem, eitem, count, irec;
	struct hecmwST_result_data* data;

	nitem = 0;
	for( j=0; j<res[0]->result->nn_component; j++ ) {
		nitem += res[0]->result->nn_dof[j];
	}
	eitem = 0;
	for( j=0; j<res[0]->result->ne_component; j++ ) {
		eitem += res[0]->result->ne_dof[j];
	}

	data = HECMW_calloc( 1, sizeof(*data) );
	data->nn_dof = HECMW_malloc( sizeof(*data->nn_dof)*res[0]->result->nn_component );
	data->node_label = HECMW_malloc( sizeof(*data->node_label)*res[0]->result->nn_component );
	data->node_val_item = HECMW_malloc( sizeof(*data->node_val_item)*nitem*glt->node_n );
	data->ne_dof = HECMW_malloc( sizeof(*data->ne_dof)*res[0]->result->ne_component );
	data->elem_label = HECMW_malloc( sizeof(*data->elem_label)*res[0]->result->ne_component );
	data->elem_val_item = HECMW_malloc( sizeof(*data->elem_val_item)*eitem*glt->elem_n );

/* for node */
	count = 0;
	g_id = -1;
	for( i=0; i<glt->node_n; i++ ) {
		if( g_id == glt->nrec[i].global ) {
			continue;
		} else {
			g_id = glt->nrec[i].global;
		}
		l_id = glt->nrec[i].local;
		area = glt->nrec[i].area;
		irec = nitem * l_id;
		for( j=0; j<res[0]->result->nn_component; j++ ) {
			for( k=0; k<res[0]->result->nn_dof[j]; k++ ) {
				data->node_val_item[count++] = res[area]->result->node_val_item[irec++];
			}
		}
	}
	data->nn_component = res[0]->result->nn_component;
	for( j=0; j<res[0]->result->nn_component; j++ ) {
		data->nn_dof[j] = res[0]->result->nn_dof[j];
		data->node_label[j] = HECMW_strdup( res[0]->result->node_label[j] );
	}
	*n_node = count/nitem;

/* for element */
	count = 0;
	g_id = -1;
	for( i=0; i<glt->elem_n; i++ ) {
		if( glt->erec[i].global < 0 ) continue;
		if( g_id == glt->erec[i].global ) {
			continue;
		} else {
			g_id = glt->erec[i].global;
		}
		l_id = glt->erec[i].local;
		area = glt->erec[i].area;
		irec = eitem * l_id;
		for( j=0; j<res[0]->result->ne_component; j++ ) {
			for( k=0; k<res[0]->result->ne_dof[j]; k++ ) {
				data->elem_val_item[count++] = res[area]->result->elem_val_item[irec++];
			}
		}
	}
	data->ne_component = res[0]->result->ne_component;
	for( j=0; j<res[0]->result->ne_component; j++ ) {
		data->ne_dof[j] = res[0]->result->ne_dof[j];
		data->elem_label[j] = HECMW_strdup( res[0]->result->elem_label[j] );
	}
	*n_elem = count/eitem;

	return data;
}


/*****************************************************************************/
/* result の削除                                                             */
/*****************************************************************************/

void fstr_free_result( fstr_res_info** res, int area_n )
{
	int i;

	if( !res ) return;

	for( i=0; i<area_n; i++ ) {
		HECMW_result_free( res[i]->result );
		HECMW_free( res[i]->node_gid );
		HECMW_free( res[i]->elem_gid );
		HECMW_free( res[i] );
	}
	HECMW_free( res );
	return;
}


/*****************************************************************************/
/* グローバルとローカル，所属領域のテーブル fstr_gl_t の作成                 */
/*****************************************************************************/

static int cmp_global_glt( const fstr_gl_rec* g1, const fstr_gl_rec* g2)
{
	return ( g1->global - g2->global );
}

typedef int (*cmp_func)(const void*, const void*);

fstr_gl_t* fstr_create_glt( struct hecmwST_local_mesh** mesh, int area_n )
{
	int i, j, count;
	int all_n, all_e;
	struct hecmwST_local_mesh* mp;
	fstr_gl_rec* nrec;
	fstr_gl_rec* erec;
	fstr_gl_t* glt;

/* for node */
	all_n = 0;
	for( i=0; i<area_n; i++ ) {
		out_log( "area:%d -- nn_internal:%d\n", i, mesh[i]->nn_internal );
		all_n += mesh[i]->nn_internal;
	}
	out_log( "total nn_internal:%d\n", all_n );

	nrec = HECMW_malloc( sizeof(fstr_gl_rec) * all_n );
	if( !nrec ) return NULL;

	count = 0;
	for( i=0; i<area_n; i++ ) {
		mp = mesh[i];
		for( j=0; j<mp->nn_internal; j++ ) {
			nrec[count].global = mp->global_node_ID[j];
			nrec[count].local = j;
			nrec[count].area = i;
			count++;
		}
	}

	qsort( nrec, all_n, sizeof(fstr_gl_rec), (cmp_func)cmp_global_glt );

/* for element */
	all_e = 0;
	for( i=0; i<area_n; i++ ) {
		out_log( "area:%d -- n_elem:%d\n", i, mesh[i]->n_elem );
		all_e += mesh[i]->n_elem;
	}
	out_log( "total n_elem:%d\n", all_e );

	erec = HECMW_malloc( sizeof(fstr_gl_rec) * all_e );
	if( !erec ) return NULL;

	count = 0;
	for( i=0; i<area_n; i++ ) {
		mp = mesh[i];
		for( j=0; j<mp->n_elem; j++ ) {
			erec[count].global = mp->global_elem_ID[j];
			erec[count].local = j;
			erec[count].area = i;
			count++;
		}
	}

	qsort( erec, all_e, sizeof(fstr_gl_rec), (cmp_func)cmp_global_glt );

/* common data */
	glt = HECMW_malloc( sizeof(fstr_gl_t) );
	if( !glt ) return NULL;
	glt->nrec = nrec;
	glt->erec = erec;
	glt->node_n = all_n;
	glt->elem_n = all_e;

	return glt;
}

fstr_gl_t* fstr_refine_glt( fstr_gl_t* glt, fstr_res_info** res, int area_n )
{
	int i, j, count;
	int all_e;
	fstr_gl_rec* erec;

/* for element */
	all_e = 0;
	for( i=0; i<area_n; i++ ) {
		all_e += res[i]->nelem_gid;
	}

	erec = HECMW_malloc( sizeof(fstr_gl_rec) * all_e );
	if( !erec ) return NULL;

	count = 0;
	for( i=0; i<area_n; i++ ) {
		for( j=0; j<res[i]->nelem_gid; j++ ) {
			erec[count].global = res[i]->elem_gid[j];
			erec[count].local = j;
			erec[count].area = i;
			count++;
		}
	}

	qsort( erec, all_e, sizeof(fstr_gl_rec), (cmp_func)cmp_global_glt );

/* common data */
	HECMW_free( glt->erec );
	glt->erec = erec;
	glt->elem_n = all_e;

	return glt;
}

/*****************************************************************************/
/* fstr_gl_t の削除                                                          */
/*****************************************************************************/

void fstr_free_gl_t( fstr_gl_t* glt )
{
	if( !glt ) return;

	HECMW_free( glt->nrec );
	HECMW_free( glt->erec );
	HECMW_free( glt );
	return;
}


/*****************************************************************************/
/* グローバルIDメッシュの作成                                                */
/*****************************************************************************/

struct hecmwST_local_mesh* fstr_create_glmesh( fstr_gl_t* glt )
{
	int i, count, g_id;
	struct hecmwST_local_mesh* mp;

	mp = HECMW_calloc( 1, sizeof(*mp) );
	mp->global_node_ID = HECMW_malloc( sizeof(*mp->global_node_ID)*glt->node_n );
	mp->global_elem_ID = HECMW_malloc( sizeof(*mp->global_elem_ID)*glt->elem_n );

	count = 0;
	g_id = -1;
	for( i=0; i<glt->node_n; i++ ) {
		if( g_id == glt->nrec[i].global ) {
			continue;
		} else {
			g_id = glt->nrec[i].global;
			mp->global_node_ID[count++] = g_id;
		}
	}
	mp->n_node = count;

	count = 0;
	g_id = -1;
	for( i=0; i<glt->elem_n; i++ ) {
		if( g_id == glt->erec[i].global ) {
			continue;
		} else {
			g_id = glt->erec[i].global;
			mp->global_elem_ID[count++] = g_id;
		}
	}
	mp->n_elem = count;

	return mp;
}


/*****************************************************************************/
/* グローバルIDメッシュの削除                                                */
/*****************************************************************************/

void fstr_free_glmesh( struct hecmwST_local_mesh* mp )
{
	if( !mp ) return;

	HECMW_free( mp->global_node_ID );
	HECMW_free( mp->global_elem_ID );
	HECMW_free( mp );
	return;
}

