/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Shin'ichi Ezure (RIST)                        *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/



#ifndef INC_HECMW_MESH_DEBUG_WRITE
#define INC_HECMW_MESH_DEBUG_WRITE

#include "hecmw_struct.h"

#ifdef DEBUG
#define HECMW_dbg_n_node( mesh ) HECMW_dbg_n_node_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_nn_internal( mesh ) HECMW_dbg_nn_internal_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_n_dof( mesh ) HECMW_dbg_n_dof_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_n_dof_grp( mesh ) HECMW_dbg_n_dof_grp_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_node( mesh ) HECMW_dbg_node_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_node_id( mesh ) HECMW_dbg_node_id_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_node_id_lid( mesh ) HECMW_dbg_node_id_lid_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_node_id_domain( mesh ) HECMW_dbg_node_id_domain_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_global_node_id( mesh ) HECMW_dbg_global_node_id_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_node_dof_item( mesh ) HECMW_dbg_node_dof_item_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_node_dof_index( mesh ) HECMW_dbg_node_dof_index_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_node_init_val_item( mesh ) HECMW_dbg_node_init_val_item_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_node_init_val_index( mesh ) HECMW_dbg_node_init_val_index_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_n_elem( mesh ) HECMW_dbg_n_elem_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_ne_internal( mesh ) HECMW_dbg_ne_internal_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_n_elem_type( mesh ) HECMW_dbg_n_elem_type_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_elem_type_item( mesh ) HECMW_dbg_elem_type_item_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_elem_type_index( mesh ) HECMW_dbg_elem_type_index_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_elem_type( mesh ) HECMW_dbg_elem_type_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_elem_node_item( mesh ) HECMW_dbg_elem_node_item_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_elem_node_index( mesh ) HECMW_dbg_elem_node_index_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_elem_id( mesh ) HECMW_dbg_elem_id_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_elem_id_lid( mesh ) HECMW_dbg_elem_id_lid_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_elem_id_domain( mesh ) HECMW_dbg_elem_id_domain_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_global_elem_id( mesh ) HECMW_dbg_global_elem_id_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_section_id( mesh ) HECMW_dbg_section_id_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_n_elem_mat_id( mesh ) HECMW_dbg_n_elem_mat_id_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_elem_mat_id_item( mesh ) HECMW_dbg_elem_mat_id_item_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_elem_mat_id_index( mesh ) HECMW_dbg_elem_mat_id_index_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_elem_internal_list( mesh ) HECMW_dbg_elem_internal_list_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_n_neighbor_pe( mesh ) HECMW_dbg_n_neighbor_pe_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_neighbor_pe( mesh ) HECMW_dbg_neighbor_pe_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_import_item( mesh ) HECMW_dbg_import_item_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_import_index( mesh ) HECMW_dbg_import_index_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_export_item( mesh ) HECMW_dbg_export_item_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_export_index( mesh ) HECMW_dbg_export_index_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_shared_item( mesh ) HECMW_dbg_shared_item_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_shared_index( mesh ) HECMW_dbg_shared_index_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_n_sect( mesh ) HECMW_dbg_n_sect_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_section( mesh ) HECMW_dbg_section_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_n_mat( mesh ) HECMW_dbg_n_mat_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_material( mesh ) HECMW_dbg_material_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_n_mpc( mesh ) HECMW_dbg_n_mpc_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_mpc( mesh ) HECMW_dbg_mpc_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_n_amp( mesh ) HECMW_dbg_n_amp_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_amp( mesh ) HECMW_dbg_amp_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_n_node_group( mesh ) HECMW_dbg_n_node_group_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_node_group( mesh ) HECMW_dbg_node_group_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_n_elem_group( mesh ) HECMW_dbg_n_elem_group_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_elem_group( mesh ) HECMW_dbg_elem_group_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_n_surf_group( mesh ) HECMW_dbg_n_surf_group_( mesh, __FILE__, __LINE__ )
#define HECMW_dbg_surf_group( mesh ) HECMW_dbg_surf_group_( mesh, __FILE__, __LINE__ )
#else
#define HECMW_dbg_n_node( mesh )
#define HECMW_dbg_nn_internal( mesh )
#define HECMW_dbg_n_dof( mesh )
#define HECMW_dbg_n_dof_grp( mesh )
#define HECMW_dbg_node( mesh )
#define HECMW_dbg_node_id( mesh )
#define HECMW_dbg_node_id_lid( mesh )
#define HECMW_dbg_node_id_domain( mesh )
#define HECMW_dbg_global_node_id( mesh )
#define HECMW_dbg_node_dof_item( mesh )
#define HECMW_dbg_node_dof_index( mesh )
#define HECMW_dbg_node_init_val_item( mesh )
#define HECMW_dbg_node_init_val_index( mesh )
#define HECMW_dbg_n_elem( mesh )
#define HECMW_dbg_ne_internal( mesh )
#define HECMW_dbg_n_elem_type( mesh )
#define HECMW_dbg_elem_type_item( mesh )
#define HECMW_dbg_elem_type_index( mesh )
#define HECMW_dbg_elem_type( mesh )
#define HECMW_dbg_elem_node_item( mesh )
#define HECMW_dbg_elem_node_index( mesh )
#define HECMW_dbg_elem_id( mesh )
#define HECMW_dbg_elem_id_lid( mesh )
#define HECMW_dbg_elem_id_domain( mesh )
#define HECMW_dbg_global_elem_id( mesh )
#define HECMW_dbg_section_id( mesh )
#define HECMW_dbg_n_elem_mat_id( mesh )
#define HECMW_dbg_elem_mat_id_item( mesh )
#define HECMW_dbg_elem_mat_id_index( mesh )
#define HECMW_dbg_elem_internal_list( mesh )
#define HECMW_dbg_n_neighbor_pe( mesh )
#define HECMW_dbg_neighbor_pe( mesh )
#define HECMW_dbg_import_item( mesh )
#define HECMW_dbg_import_index( mesh )
#define HECMW_dbg_export_item( mesh )
#define HECMW_dbg_export_index( mesh )
#define HECMW_dbg_shared_item( mesh )
#define HECMW_dbg_shared_index( mesh )
#define HECMW_dbg_n_sect( mesh )
#define HECMW_dbg_section( mesh )
#define HECMW_dbg_n_mat( mesh )
#define HECMW_dbg_material( mesh )
#define HECMW_dbg_n_mpc( mesh )
#define HECMW_dbg_mpc( mesh )
#define HECMW_dbg_n_amp( mesh )
#define HECMW_dbg_amp( mesh )
#define HECMW_dbg_n_node_group( mesh )
#define HECMW_dbg_node_group( mesh )
#define HECMW_dbg_n_elem_group( mesh )
#define HECMW_dbg_elem_group( mesh )
#define HECMW_dbg_n_surf_group( mesh )
#define HECMW_dbg_surf_group( mesh )
#endif

extern void
HECMW_dbg_n_node_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_nn_internal_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_n_dof_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_n_dof_grp_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_node_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_node_id_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_node_id_lid_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_node_id_domain_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_global_node_id_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_node_dof_item_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_node_dof_index_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_node_init_val_item_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_node_init_val_index_( struct hecmwST_local_mesh *local_mesh, char *file, int line );

extern void
HECMW_dbg_n_elem_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_ne_internal_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_n_elem_type_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_elem_type_item_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_elem_type_index_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_elem_type_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_elem_node_item_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_elem_node_index_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_elem_id_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_elem_id_lid_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_elem_id_domain_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_global_elem_id_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_elem_internal_list_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_section_id_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_n_elem_mat_id_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_elem_mat_id_item_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_elem_mat_id_index_( struct hecmwST_local_mesh *local_mesh, char *file, int line );

extern void
HECMW_dbg_n_neighbor_pe_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_neighbor_pe_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_import_item_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_import_index_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_export_item_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_export_index_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_shared_item_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_shared_index_( struct hecmwST_local_mesh *local_mesh, char *file, int line );

extern void
HECMW_dbg_n_sect_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_section_( struct hecmwST_local_mesh *local_mesh, char *file, int line );

extern void
HECMW_dbg_n_mat_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_material_( struct hecmwST_local_mesh *local_mesh, char *file, int line );

extern void
HECMW_dbg_n_mpc_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_mpc_( struct hecmwST_local_mesh *local_mesh, char *file, int line );

extern void
HECMW_dbg_n_amp_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_amp_( struct hecmwST_local_mesh *local_mesh, char *file, int line );

extern void
HECMW_dbg_n_node_group_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_node_group_( struct hecmwST_local_mesh *local_mesh, char *file, int line );

extern void
HECMW_dbg_n_elem_group_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_elem_group_( struct hecmwST_local_mesh *local_mesh, char *file, int line );

extern void
HECMW_dbg_n_surf_group_( struct hecmwST_local_mesh *local_mesh, char *file, int line );
extern void
HECMW_dbg_surf_group_( struct hecmwST_local_mesh *local_mesh, char *file, int line );

#endif
