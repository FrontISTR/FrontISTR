/*=====================================================================*
!                                                                      !
! Software Name : HEC-MW Ver 4.3                                      !
!                                                                      !
!      Module Name : Visualizer Utility                                !
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

static struct hecmwST_local_mesh *static_mesh;
static struct hecmwST_result_data *static_data;
static int *external_node;

static void memory_error(char *var)
{
	fprintf(stderr, "#### HEC-MW-VIS-E0001:There is no enough memory allocated for variable %s\n", var);
	mw_finalize_();
	abort();
}

struct hecmwST_local_mesh* MW3_get_mesh( int level, int partID )
{
	struct hecmwST_local_mesh *mesh;
	int *part, *nnode, *nelem, *neighbor, *import, *export, *importid, *reorderid, *sendbuff, *recvbuff;
	int npart, total_nnode, total_nelem, nproc, iproc, inode, etype, mynode;
	int i, j, k, n, m, count, ncount, ecount;
	HECMW_Request request;

	mynode = mw_get_rank_();

	mw_select_assemble_model_( &level );

	nproc = mw_get_num_of_process_();
	neighbor = calloc( nproc, sizeof(*neighbor) );
	if( !neighbor ) memory_error( "neighbor" );
	if( partID < 0 ) {
		npart = mw_get_num_of_mesh_part_();
		part = malloc( sizeof(*part) * npart );
		if( !part ) memory_error( "part" );
		nnode = malloc( sizeof(*nnode) * (npart+1) );
		if( !nnode ) memory_error( "nnode" );
		nelem = malloc( sizeof(*nelem) * (npart+1) );
		if( !nelem ) memory_error( "nelem" );
		nnode[0] = 0;
		nelem[0] = 0;
		for( i = 0; i < npart; i++ ) {
			part[i] = i;
			mw_select_mesh_part_( &part[i] );
			nnode[i+1] = nnode[i] + mw_get_num_of_node_();
			nelem[i+1] = nelem[i] + mw_get_num_of_element_();
			n = mw_get_num_of_neibpe_( &part[i] );
			for( j = 0; j < n; j++ ) {
				iproc = mw_get_transrank_( &part[i], &j );
				neighbor[iproc] = 1;
			}
		}
	} else {
		npart = 1;
		part = malloc( sizeof(*part) );
		if( !part ) memory_error( "part" );
		nnode = malloc( sizeof(*nnode) * 2 );
		if( !nnode ) memory_error( "nnode" );
		nelem = malloc( sizeof(*nelem) * 2 );
		if( !nelem ) memory_error( "nelem" );
		nnode[0] = 0;
		nelem[0] = 0;
		part[0] = partID;
		mw_select_mesh_part_( &part[0] );
		nnode[1] = mw_get_num_of_node_();
		nelem[1] = mw_get_num_of_element_();
		n = mw_get_num_of_neibpe_( &part[0] );
		for( i = 0; i < n; i++ ) {
			iproc = mw_get_transrank_( &part[0], &i );
			neighbor[iproc] = 1;
		}
	}
	total_nnode = nnode[npart];
	total_nelem = nelem[npart];

	mesh = calloc( 1, sizeof(*mesh) );
	if( !mesh ) memory_error( "mesh" );
	mesh->n_node = total_nnode;
	i = 0;
	inode = mw_get_node_id_( &i );
	mesh->n_dof = mw_get_dof_( &inode );
	mesh->n_elem = total_nelem;
	mesh->ne_internal = total_nelem;
	mesh->HECMW_COMM = mw_mpi_comm_();
	mesh->my_rank = mynode;

	for( i = 0; i < nproc; i++ ) {
		if( neighbor[i] ) mesh->n_neighbor_pe++;
	}
	import = calloc( total_nnode, sizeof(*import) );
	if( !import ) memory_error( "import" );
	export = calloc( total_nnode, sizeof(*export) );
	if( !export ) memory_error( "export" );
	importid = malloc( sizeof(*importid) * total_nnode );
	if( !importid ) memory_error( "importid" );
	external_node = calloc( total_nnode, sizeof(*external_node) );
	if( !external_node ) memory_error( "external_node" );

	if( mesh->n_neighbor_pe ) {
		mesh->neighbor_pe = malloc( sizeof(*mesh->neighbor_pe) * mesh->n_neighbor_pe );
		if( !mesh->neighbor_pe ) memory_error( "mesh->neighbor_pe" );
		count = 0;
		for( i = 0; i < nproc; i++ ) {
			if( neighbor[i] ) mesh->neighbor_pe[count++] = i;
		}

		count = 0;
		for( i = 0; i < npart; i++ ) {
			mw_select_mesh_part_( &part[i] );
			n = mw_get_num_of_comm_mesh_();
			for( j = 0; j < n; j++ ) {
				iproc = mw_get_transrank_( &part[i], &j );
				m = mw_get_num_of_comm_node_( &j );
				if( iproc > mynode ) {
					sendbuff = malloc( sizeof(*sendbuff) * m );
					for( k = 0; k < m; k++ ) {
						inode = mw_get_node_id_comm_node_( &j, &k );
						inode = mw_get_node_index_( &inode );
						export[count+inode] = iproc + 1;
						sendbuff[k] = count + inode;
					}
					HECMW_Isend( sendbuff, m, HECMW_INT, iproc, 0, mesh->HECMW_COMM, &request );
					free( sendbuff );
				}
			}
			for( j = 0; j < n; j++ ) {
				iproc = mw_get_transrank_( &part[i], &j );
				m = mw_get_num_of_comm_node_( &j );
				if( iproc < mynode ) {
					recvbuff = malloc( sizeof(*recvbuff) * m );
					HECMW_Irecv( recvbuff, m, HECMW_INT, iproc, 0, mesh->HECMW_COMM, &request );
					for( k = 0; k < m; k++ ) {
						inode = mw_get_node_id_comm_node_( &j, &k );
						inode = mw_get_node_index_( &inode );
						import[count+inode] = iproc + 1;
						importid[count+inode] = recvbuff[k];
						external_node[count+inode] = iproc + 1;
					}
					free( recvbuff );
				}
			}
			count += mw_get_num_of_node_();
		}

		mesh->import_index = calloc( mesh->n_neighbor_pe + 1, sizeof(*mesh->import_index) );
		if( !mesh->import_index ) memory_error( "mesh->import_index" );
		mesh->export_index = calloc( mesh->n_neighbor_pe + 1, sizeof(*mesh->export_index) );
		if( !mesh->export_index ) memory_error( "mesh->export_index" );
		for( i = 0; i < total_nnode; i++ ) {
			if( import[i] ) {
				for( j = 0; j < mesh->n_neighbor_pe; j++ ) {
					if( import[i] == (mesh->neighbor_pe[j]+1) ) mesh->import_index[j]++;
				}
			}
			if( export[i] ) {
				for( j = 0; j < mesh->n_neighbor_pe; j++ ) {
					if( export[i] == (mesh->neighbor_pe[j]+1) ) mesh->export_index[j]++;
				}
			}
		}
		for( i = mesh->n_neighbor_pe; i >= 0; i-- ) {
			mesh->import_index[i] = 0;
			mesh->export_index[i] = 0;
			for( j = 0; j < i; j++ ) {
				mesh->import_index[i] += mesh->import_index[j];
				mesh->export_index[i] += mesh->export_index[j];
			}
		}

		if( mesh->import_index[mesh->n_neighbor_pe] ) {
			mesh->import_item = malloc( sizeof(*mesh->import_item) * mesh->import_index[mesh->n_neighbor_pe] );
			if( !mesh->import_item ) memory_error( "mesh->import_item" );
		}
		if( mesh->export_index[mesh->n_neighbor_pe] ) {
			mesh->export_item = malloc( sizeof(*mesh->export_item) * mesh->export_index[mesh->n_neighbor_pe] );
			if( !mesh->export_item ) memory_error( "mesh->export_item" );
		}
		ncount = 0;
		ecount = 0;
		for( i = 0; i < mesh->n_neighbor_pe; i++ ) {
			for( j = 0; j < total_nnode; j++ ) {
				if( import[j] == (mesh->neighbor_pe[i]+1) ) mesh->import_item[ncount++] = j + 1;
				if( export[j] == (mesh->neighbor_pe[i]+1) ) mesh->export_item[ecount++] = j + 1;
			}
		}

		mesh->nn_internal = total_nnode - mesh->import_index[mesh->n_neighbor_pe];
	} else {
		mesh->nn_internal = total_nnode;
	}

	mesh->node_internal_list = malloc( sizeof(*mesh->node_internal_list) * mesh->nn_internal );
	if( !mesh->node_internal_list ) memory_error( "mesh->node_internal_list" );
	mesh->node_ID = malloc( sizeof(*mesh->node_ID) * total_nnode * 2 );
	if( !mesh->node_ID ) memory_error( "mesh->node_ID" );
	mesh->global_node_ID = malloc( sizeof(*mesh->global_node_ID) * total_nnode );
	if( !mesh->global_node_ID ) memory_error( "mesh->global_node_ID" );
	mesh->node = malloc( sizeof(*mesh->node) * total_nnode * 3 );
	if( !mesh->node ) memory_error( "mesh->node" );
	reorderid = malloc( sizeof(*reorderid) * total_nnode );
	if( !reorderid ) memory_error( "reorderid" );
	ncount = 0;
	ecount = mesh->nn_internal;
	count = 0;
	for( i = 0; i < npart; i++ ) {
		mw_select_mesh_part_( &part[i] );
		n = mw_get_num_of_node_();
		for( j = 0; j < n; j++ ) {
			if( import[ncount] ) {
				mesh->node_ID[2*ecount] = importid[ncount] + 1;
				mesh->node_ID[2*ecount+1] = import[ncount] - 1;
				inode = mw_get_node_id_( &j );
				mesh->global_node_ID[ecount] = inode;
				mw_get_node_coord_( &inode, &mesh->node[3*ecount], &mesh->node[3*ecount+1], &mesh->node[3*ecount+2] );
				reorderid[ncount] = ecount;
				ecount++;
			} else {
				mesh->node_internal_list[count] = count + 1;
				mesh->node_ID[2*count] = count + 1;
				mesh->node_ID[2*count+1] = mesh->my_rank;
				inode = mw_get_node_id_( &j );
				mesh->global_node_ID[count] = inode;
				mw_get_node_coord_( &inode, &mesh->node[3*count], &mesh->node[3*count+1], &mesh->node[3*count+2] );
				reorderid[ncount] = count;
				count++;
			}
			ncount++;
		}
	}

	mesh->elem_internal_list = malloc( sizeof(*mesh->elem_internal_list) * mesh->ne_internal );
	if( !mesh->elem_internal_list ) memory_error( "mesh->elem_internal_list" );
	mesh->elem_ID = malloc( sizeof(*mesh->elem_ID) * total_nelem * 2 );
	if( !mesh->elem_ID ) memory_error( "mesh->elem_ID" );
	mesh->global_elem_ID = malloc( sizeof(*mesh->global_elem_ID) * total_nelem );
	if( !mesh->global_elem_ID ) memory_error( "mesh->global_elem_ID" );
	mesh->elem_type = malloc( sizeof(*mesh->elem_type) * total_nelem );
	if( !mesh->elem_type ) memory_error( "mesh->elem_type" );
	mesh->elem_node_index = malloc( sizeof(*mesh->elem_node_index) * (total_nelem+1) );
	if( !mesh->elem_node_index ) memory_error( "mesh->elem_node_index" );
	ecount = 0;
	mesh->elem_node_index[0] = 0;
	for( i = 0; i < npart; i++ ) {
		mw_select_mesh_part_( &part[i] );
		n = mw_get_num_of_element_();
		for( j = 0; j < n; j++ ) {
			mesh->elem_ID[2*ecount] = ecount + 1;
			mesh->elem_ID[2*ecount+1] = mesh->my_rank;
			mesh->elem_internal_list[ecount] = ecount + 1;
			mesh->global_elem_ID[ecount] = mw_get_element_id_( &j );
			mw_select_element_( &j );
			etype = mw_get_element_type_();
			mesh->elem_type[ecount] = mw_mw3_elemtype_to_fistr_elemtype_( &etype );
			mesh->elem_node_index[ecount+1] = mesh->elem_node_index[ecount] + mw_get_num_of_element_vert_();
			ecount++;
		}
	}

	mesh->elem_node_item = malloc( sizeof(*mesh->elem_node_item) * mesh->elem_node_index[total_nelem] );
	if( !mesh->elem_node_item ) memory_error( "mesh->elem_node_item" );
	count = 0;
	for( i = 0; i < npart; i++ ) {
		mw_select_mesh_part_( &part[i] );
		n = mw_get_num_of_element_();
		for( j = 0; j < n; j++ ) {
			mw_select_element_( &j );
			mw_get_element_vert_node_id_( &mesh->elem_node_item[mesh->elem_node_index[count]] );
			for( k = mesh->elem_node_index[count]; k < mesh->elem_node_index[count+1]; k++ ) {
				inode = mw_get_node_index_( &mesh->elem_node_item[k] ) + nnode[i];
				mesh->elem_node_item[k] = reorderid[inode] + 1;
			}
			count++;
		}
	}

/* begin for FEMAP neutral file */
	mesh->section_ID = malloc( sizeof(*mesh->section_ID) * total_nelem );
	if( !mesh->section_ID ) memory_error( "mesh->section_ID" );
	for( i = 0; i < total_nelem; i++ ) mesh->section_ID[i] = 0;

	mesh->section = calloc( 1, sizeof(struct hecmwST_section) );
	if( !mesh->section ) memory_error( "mesh->section" );
	mesh->section->n_sect = 1;
	mesh->section->sect_type = malloc( sizeof(*mesh->section->sect_type) );
	if( !mesh->section->sect_type ) memory_error( "mesh->section->sect_type" );
	mesh->section->sect_type[0] = 1;
	mesh->section->sect_opt = malloc( sizeof(*mesh->section->sect_opt) );
	if( !mesh->section->sect_opt ) memory_error( "mesh->section->sect_opt" );
	mesh->section->sect_opt[0] = 0;
	mesh->section->sect_mat_ID_index = malloc( sizeof(*mesh->section->sect_mat_ID_index) * 2 );
	if( !mesh->section->sect_mat_ID_index ) memory_error( "mesh->section->sect_mat_ID_index" );
	mesh->section->sect_mat_ID_index[0] = 0;
	mesh->section->sect_mat_ID_index[1] = 1;
	mesh->section->sect_mat_ID_item = malloc( sizeof(*mesh->section->sect_mat_ID_item) );
	if( !mesh->section->sect_mat_ID_item ) memory_error( "mesh->section->sect_mat_ID_item" );
	mesh->section->sect_mat_ID_item[0] = 0;

	mesh->material = calloc( 1, sizeof(struct hecmwST_material) );
	if( !mesh->material ) memory_error( "mesh->material" );
	mesh->material->n_mat = 1;
	mesh->material->mat_val = malloc( sizeof(*mesh->material->mat_val) * 4 );
	if( !mesh->material->mat_val ) memory_error( "mesh->material->mat_val" );
	for( i = 0; i < 4; i++ ) mesh->material->mat_val[i] = 0.0;
/* end for FEMAP neutral file */

	free( part );
	free( nnode );
	free( nelem );
	free( neighbor );
	free( import );
	free( importid );
	free( reorderid );
	free( export );

	return mesh;
}

void MW3_free_mesh( struct hecmwST_local_mesh* mesh )
{
	free( mesh->material->mat_val );
	free( mesh->material );
	free( mesh->section->sect_mat_ID_item );
	free( mesh->section->sect_mat_ID_index );
	free( mesh->section->sect_opt );
	free( mesh->section->sect_type );
	free( mesh->section );

	free( mesh->section_ID );
	if( mesh->n_neighbor_pe ) {
		free( mesh->export_item );
		free( mesh->export_index );
		free( mesh->import_item );
		free( mesh->import_index );
		free( mesh->neighbor_pe );
	}
	free( mesh->elem_node_item );
	free( mesh->elem_node_index );
	free( mesh->elem_type );
	free( mesh->global_elem_ID );
	free( mesh->elem_ID );
	free( mesh->elem_internal_list );
	free( mesh->node );
	free( mesh->global_node_ID );
	free( mesh->node_ID );
	free( mesh );
	free( external_node );

	return;
}

int MW3_visualize_init( void )
{
	HECMW_Comm hecmw_comm;
	HECMW_Group hecmw_group;
	int comm_size, comm_rank;
	int ierr;

	hecmw_comm = mw_mpi_comm_();
	comm_size = mw_get_num_of_process_();
	comm_rank = mw_get_rank_();
	hecmw_group = 0;
	hecmw_comm_init_if( &hecmw_comm, &comm_size, &comm_rank, &hecmw_group);

	ierr = HECMW_ctrl_init();
	if( ierr < 0 ) return ierr;

	HECMW_visualize_init();

	return 0;
}

int MW3_visualize_init_ex( char* ctrlfile )
{
	HECMW_Comm hecmw_comm;
	HECMW_Group hecmw_group;
	int comm_size, comm_rank;
	int ierr;

	hecmw_comm = mw_mpi_comm_();
	comm_size = mw_get_num_of_process_();
	comm_rank = mw_get_rank_();
	hecmw_group = 0;
	hecmw_comm_init_if( &hecmw_comm, &comm_size, &comm_rank, &hecmw_group);

	if( strlen(ctrlfile) ) {
		ierr = HECMW_ctrl_init_ex( ctrlfile );
	} else {
		ierr = HECMW_ctrl_init();
	}
	if( ierr < 0 ) return ierr;

	HECMW_visualize_init();

	return 0;
}

int MW3_visualize_finalize( void )
{
	HECMW_ctrl_finalize();
	HECMW_visualize_finalize();

	return 0;
}

int MW3_result_init( struct hecmwST_local_mesh* mesh, int tstep, char* header )
{
	HECMW_Comm hecmw_comm;
	HECMW_Group hecmw_group;
	int comm_size, comm_rank;
	int ierr;

	hecmw_comm = mw_mpi_comm_();
	comm_size = mw_get_num_of_process_();
	comm_rank = mw_get_rank_();
	hecmw_group = 0;
	hecmw_comm_init_if( &hecmw_comm, &comm_size, &comm_rank, &hecmw_group);

	ierr = HECMW_ctrl_init();
	if( ierr < 0 ) return ierr;

	HECMW_result_init( mesh, tstep, tstep, header);

	return 0;
}

int MW3_result_init_ex( char* ctrlfile, struct hecmwST_local_mesh* mesh, int tstep, char* header )
{
	HECMW_Comm hecmw_comm;
	HECMW_Group hecmw_group;
	int comm_size, comm_rank;
	int ierr;

	hecmw_comm = mw_mpi_comm_();
	comm_size = mw_get_num_of_process_();
	comm_rank = mw_get_rank_();
	hecmw_group = 0;
	hecmw_comm_init_if( &hecmw_comm, &comm_size, &comm_rank, &hecmw_group);

	if( strlen(ctrlfile) ) {
		ierr = HECMW_ctrl_init_ex( ctrlfile );
	} else {
		ierr = HECMW_ctrl_init();
	}
	if( ierr < 0 ) return ierr;

	HECMW_result_init( mesh, tstep, tstep, header);

	return 0;
}

int MW3_result_finalize( void )
{
	HECMW_ctrl_finalize();
	HECMW_result_finalize();

	return 0;
}

void mw3_visualize_init_( int* level, int* partID )
{
	if( MW3_visualize_init() < 0 ) HECMW_abort( HECMW_comm_get_comm() );
	static_mesh = MW3_get_mesh( *level, *partID );
	static_data = malloc( sizeof(*static_data) );
	HECMW_result_copy_f2c_init( static_data, static_mesh->n_node, static_mesh->n_elem );
}

void mw3_visualize_init_ex_( char* ctrlfile, int* level, int* partID, int len )
{
	char filename[HECMW_FILENAME_LEN+1];

	if( HECMW_strcpy_f2c_r( ctrlfile, len, filename, sizeof(filename) ) == NULL ) return;
	if( MW3_visualize_init_ex( filename ) < 0 ) HECMW_abort( HECMW_comm_get_comm() );
	static_mesh = MW3_get_mesh( *level, *partID );
	static_data = malloc( sizeof(*static_data) );
	HECMW_result_copy_f2c_init( static_data, static_mesh->n_node, static_mesh->n_elem );
}

void mw3_visualize_( int* timestep, int* max_timestep, int* is_force )
{
	double *data;
	int n, size, i, j, count, icount, ecount;

	n = 0;
	for( i = 0; i < static_data->nn_component; i++ ) n += static_data->nn_dof[i];
	size = sizeof(double) * static_mesh->n_node * n;
	data = malloc( size );
	if( !data ) memory_error( "node_conversion_data" );
	memcpy( data, static_data->node_val_item, size );
	count = 0;
	icount = 0;
	ecount = static_mesh->nn_internal * n;
	for( i = 0; i < static_mesh->n_node; i++ ) {
		if( external_node[i] ) {
			for( j = 0; j < n; j++ ) static_data->node_val_item[ecount++] = data[count++];
		} else {
			for( j = 0; j < n; j++ ) static_data->node_val_item[icount++] = data[count++];
		}
	}
	free( data );

	HECMW_visualize( static_mesh, static_data, *timestep, *max_timestep, *is_force );
}

void mw3_visualize_finalize_( void )
{
	if( MW3_visualize_finalize() < 0 ) HECMW_abort( HECMW_comm_get_comm() );
	HECMW_result_free( static_data );
	HECMW_result_copy_f2c_finalize();
	MW3_free_mesh( static_mesh );
}

void mw3_result_init_( int* level, int* partID, int* tstep, char* header, int len )
{
	char headername[HECMW_HEADER_LEN+1];

	if( HECMW_strcpy_f2c_r( header, len, headername, sizeof(headername) ) == NULL ) return;

	static_mesh = MW3_get_mesh( *level, *partID );

	if( MW3_result_init( static_mesh, *tstep, headername ) < 0 ) HECMW_abort( HECMW_comm_get_comm() );
}

void mw3_result_init_ex_( char* ctrlfile, int* level, int* partID, int* tstep, char* header, int lenc, int lenh )
{
	char filename[HECMW_FILENAME_LEN+1];
	char headername[HECMW_HEADER_LEN+1];

	if( HECMW_strcpy_f2c_r( ctrlfile, lenc, filename, sizeof(filename) ) == NULL ) return;
	if( HECMW_strcpy_f2c_r( header, lenh, headername, sizeof(headername) ) == NULL ) return;

	static_mesh = MW3_get_mesh( *level, *partID );

	if( MW3_result_init_ex( filename, static_mesh, *tstep, headername ) < 0 ) HECMW_abort( HECMW_comm_get_comm() );
}

void mw3_result_add_( int* node_or_elem, int* n_dof, char* label, double* ptr, int len)
{
	double *data;
	int size, i, j, count, icount, ecount;
	char labelname[HECMW_NAME_LEN+1];

	if( HECMW_strcpy_f2c_r( label, len, labelname, sizeof(labelname) ) == NULL ) return;

	size = sizeof(double) * static_mesh->n_node * (*n_dof);
	data = malloc( size );
	if( !data ) memory_error( "node_conversion_data" );
	if ( *node_or_elem == 1 ) {
		count = 0;
		icount = 0;
		ecount = static_mesh->nn_internal * (*n_dof);
		for( i = 0; i < static_mesh->n_node; i++ ) {
			if( external_node[i] ) {
				for( j = 0; j < *n_dof; j++ ) data[ecount++] = ptr[count++];
			} else {
				for( j = 0; j < *n_dof; j++ ) data[icount++] = ptr[count++];
			}
		}
		if( HECMW_result_add( *node_or_elem, *n_dof, labelname, data ) < 0 ) HECMW_abort( HECMW_comm_get_comm() );
	} else {
		memcpy( data, ptr, size );
		if( HECMW_result_add( *node_or_elem, *n_dof, labelname, data ) < 0 ) HECMW_abort( HECMW_comm_get_comm() );
	}
}

void mw3_result_write_by_name_( char* name_ID, int len )
{
	char name[HECMW_NAME_LEN+1];

	if( HECMW_strcpy_f2c_r( name_ID, len, name, sizeof(name) ) == NULL ) return;
	if( HECMW_result_write_by_name( name ) < 0 ) HECMW_abort( HECMW_comm_get_comm() );
}

void mw3_result_finalize_( void )
{
	if( MW3_result_finalize() < 0 ) HECMW_abort( HECMW_comm_get_comm() );
	MW3_free_mesh( static_mesh );
}
