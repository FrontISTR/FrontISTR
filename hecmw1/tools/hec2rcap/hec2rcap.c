/*
	hec2rcap ver.1.0
	2007.01.18 by N.Imai
	-------------------------------------
	HEC-MW Mesh -> REVOCAP Mesh Converter
*/


#include <stdio.h>
#include <string.h>
#include "hecmw.h"

#define DEF_MESH_ID  "fstrMSH"
#define DEF_RCAP_ID  "revocap"

static
char* rcap_fname_header = NULL;


/*=============================================================================
 * convert
 *===========================================================================*/


static
int conv_node( struct hecmwST_local_mesh* mesh, FILE* fp )
{
	int i, j;
	int node_n = mesh->n_node;

	fprintf( fp, "Number_of_Node_Ids %d\n",  node_n );
	j = 0;
	for( i=0; i<node_n; i++, j+=3) {
		int nid = mesh->global_node_ID[i];
		double x = mesh->node[j];
		double y = mesh->node[j+1];
		double z = mesh->node[j+2];
		fprintf( fp,"%d  %lf %lf %lf\n", nid, x, y, z );
	}
	return 0;
}

static
int conv_elem( struct hecmwST_local_mesh* mesh, FILE* fp )
{
	int i, j;
	int is,ie, icel;
	int etype;
	char* rcap_etype;
	int nn, elem_n;
	int etype_count;

	int tbl341[] = { 0,1,2,3 };
	int tbl342[] = { 0,1,2,3,6,5,7,4,9,8 };
	int tbl361[] = { 0,1,2,3,4,5,6,7 };
	int tbl362[] = { 0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15 };
	int tbl351[] = { 0,1,2,3,4,5 };
	int tbl371[] = { 4,0,3,2,1 };
	int* tbl;

	for( i= 0, etype_count= 0, elem_n= 0; i<mesh->n_elem_type; i++ ) {
		is = mesh->elem_type_index[i];
		ie = mesh->elem_type_index[i+1];
		etype = mesh->elem_type[is];
		if( etype >= 900 ) {
			/* no structural element */
			continue;
		}
		etype_count++;
		elem_n += ie - is;
	}
	if( etype_count == 1 ) {
		switch( etype ) {
		case 341: nn =  4; rcap_etype =  "Tet4"; tbl=tbl341; break;
		case 342: nn = 10; rcap_etype = "Tet10"; tbl=tbl342; break;
		case 361: nn =  8; rcap_etype =  "Hex8"; tbl=tbl361; break;
		case 362: nn = 20; rcap_etype = "Hex20"; tbl=tbl362; break;
		case 351: nn =  6; rcap_etype =  "Wed6"; tbl=tbl351; break;
		case 371: nn =  5; rcap_etype =  "Pyr5"; tbl=tbl371; break;
		default:
			fprintf( stderr, "#Error : not supported element type %d\n", etype );
			return 1;
		}
		fprintf( fp, "Element_Type %s\n", rcap_etype );
	} else {
		fprintf( fp, "Element_Type Multi\n" );
	}
	fprintf( fp, "Number_of_Elemen_Ids %d\n",  elem_n );

	for( i= 0; i<mesh->n_elem_type; i++ ) {
		is = mesh->elem_type_index[i];
		ie = mesh->elem_type_index[i+1];
		etype = mesh->elem_type[is];
		if( etype >= 900 ) {
			/* no structural element */
			continue;
		}
		for( icel = is; icel<ie; icel++ ) {
			int eid = mesh->global_elem_ID[icel];
			int n0 = mesh->elem_node_index[icel];
			fprintf( fp, "%d", eid );
			if( etype_count > 1 ) {
				switch( etype ) {
				case 341: nn =  4; rcap_etype =  "Tet4"; tbl=tbl341; break;
				case 342: nn = 10; rcap_etype = "Tet10"; tbl=tbl342; break;
				case 361: nn =  8; rcap_etype =  "Hex8"; tbl=tbl361; break;
				case 362: nn = 20; rcap_etype = "Hex20"; tbl=tbl362; break;
				case 351: nn =  6; rcap_etype =  "Wed6"; tbl=tbl351; break;
				case 371: nn =  5; rcap_etype =  "Pyr5"; tbl=tbl371; break;
				default:
					fprintf( stderr, "#Error : not supported element type %d\n", etype );
					return 1;
				}
				fprintf( fp, " %s", rcap_etype);
			}
			for( j=0; j<nn; j++) {
				int id = mesh->elem_node_item[n0+tbl[j]];
				int nid = mesh->global_node_ID[id-1];
				fprintf( fp, " %d", nid );
			}
			fprintf( fp, "\n" );
		}
	}
	return 0;
}


static
int convert( struct hecmwST_local_mesh* mesh, FILE* fp )
{
	fprintf( fp, "Solid_PartID %d\n", mesh->my_rank );
	if( conv_elem( mesh, fp ) ) return 1;
	if( conv_node( mesh, fp ) ) return 1;
	return 0;
}


/*=============================================================================
 * main & etc
 *===========================================================================*/

/*
static
void help(void)
{
}
*/

static
int set_param(  int argc, char** argv )
{
	if( argc < 2 ) {
		fprintf( stderr, "#Error : revocap mesh file name is required\n" );
		return 1;
	}
	rcap_fname_header = argv[1];
	return 0;
}


int main( int argc, char** argv )
{
	/* input mesh */
	char* name_ID = DEF_MESH_ID;
	struct hecmwST_local_mesh* mesh = NULL;

	/* output mesh */
	char rcap_fname[HECMW_FILENAME_LEN +1];
	struct hecmw_ctrl_meshfiles *files;
	char dirname[HECMW_HEADER_LEN+1];
	char buff[HECMW_HEADER_LEN+1];
	char *ptoken, *ntoken;
	FILE* fp;

	if( HECMW_init(&argc, &argv) ) return 1;
	if( set_param(argc,argv) ) {
		HECMW_abort(HECMW_comm_get_comm());
	}

	mesh = HECMW_get_mesh( name_ID );
	if( mesh == NULL ) {
		HECMW_abort(HECMW_comm_get_comm());
	}
	printf( "domain #%d : get_mesh ok\n", mesh->my_rank );

	files = HECMW_ctrl_get_meshfiles_header( name_ID );
	if( files == NULL ) {
		HECMW_abort(HECMW_comm_get_comm());
	}
	strcpy( buff, files->meshfiles[0].filename );
	strcpy( dirname, "" );
	ptoken = strtok (buff, "/" );
	ntoken = strtok( NULL, "/" );
	while( ntoken ) {
		strcat( dirname, ptoken );
		strcat( dirname, "/" );
		ptoken = ntoken;
		ntoken = strtok( NULL, "/" );
	}
	sprintf( rcap_fname, "%s%s.%d", dirname, rcap_fname_header, mesh->my_rank );

	fp = fopen( rcap_fname, "w" );
	if( !fp ) {
		fprintf( stderr, "#Error at domain #%d : Cannot open revocap file %s\n",
				 mesh->my_rank, rcap_fname );
		HECMW_abort(HECMW_comm_get_comm());
	}

	if( convert( mesh, fp ) ) {
		HECMW_abort(HECMW_comm_get_comm());
	}

	fclose( fp );
	printf( "completed\n" );
	HECMW_finalize();
	return 0;
}
