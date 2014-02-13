/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuaki Sakane (RIST)                         *
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
#include <ctype.h>
#include "hecmw_struct.h"
#include "hecmw_util.h"
#include "hecmw_dist_alloc.h"

/*============================================================================*/
/*  initialize structures for local mesh                                      */
/*============================================================================*/
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  global information                                                        */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int
init_struct_global( struct hecmwST_local_mesh *mesh )
{
  HECMW_assert( mesh );

  memset( mesh->gridfile, 0, sizeof(mesh->gridfile) );
  mesh->hecmw_n_file = 0;
  mesh->files        = NULL;
  memset( mesh->header, 0, sizeof(mesh->header) );

  mesh->hecmw_flag_adapt     = 0;
  mesh->hecmw_flag_initcon   = 0;
  mesh->hecmw_flag_parttype  = 0;
  mesh->hecmw_flag_partdepth = 0;
  mesh->hecmw_flag_version   = 0;

  mesh->zero_temp = 0.0;

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  node information                                                          */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int
init_struct_node( struct hecmwST_local_mesh *mesh )
{
  HECMW_assert( mesh );

  mesh->n_node             = 0;
  mesh->n_node_gross       = 0;
  mesh->nn_internal        = 0;
  mesh->node_internal_list = NULL;

  mesh->node           = NULL;
  mesh->node_ID        = NULL;
  mesh->global_node_ID = NULL;

  mesh->n_dof          = 0;
  mesh->n_dof_grp      = 0;
  mesh->n_dof_tot      = 0;
  mesh->node_dof_index = NULL;
  mesh->node_dof_item  = NULL;

  mesh->node_val_index = NULL;
  mesh->node_val_item  = NULL;

  mesh->node_init_val_index = NULL;
  mesh->node_init_val_item  = NULL;

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  element information                                                       */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int
init_struct_elem( struct hecmwST_local_mesh *mesh )
{
  HECMW_assert( mesh );

  mesh->n_elem             = 0;
  mesh->n_elem_gross       = 0;
  mesh->ne_internal        = 0;
  mesh->elem_internal_list = NULL;

  mesh->elem_ID        = NULL;
  mesh->global_elem_ID = NULL;

  mesh->n_elem_type     = 0;
  mesh->elem_type       = NULL;
  mesh->elem_type_index = NULL;
  mesh->elem_type_item  = NULL;

  mesh->elem_node_index = NULL;
  mesh->elem_node_item  = NULL;

  mesh->section_ID = NULL;

  mesh->n_elem_mat_ID     = 0;
  mesh->elem_mat_ID_index = NULL;
  mesh->elem_mat_ID_item  = NULL;

  mesh->elem_mat_int_index = NULL;
  mesh->elem_mat_int_val   = NULL;

  mesh->elem_val_index = NULL;
  mesh->elem_val_item  = NULL;

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  domain & communication information                                        */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int
init_struct_comm( struct hecmwST_local_mesh *mesh )
{
  HECMW_assert( mesh );

  mesh->zero        = 0;
  mesh->PETOT       = 0;
  mesh->PEsmpTOT    = 0;
  mesh->my_rank     = 0;
  mesh->errnof      = 0;
  mesh->n_subdomain = 0;

  mesh->n_neighbor_pe = 0;
  mesh->neighbor_pe   = NULL;

  mesh->import_index = NULL;
  mesh->import_item  = NULL;
  mesh->export_index = NULL;
  mesh->export_item  = NULL;
  mesh->shared_index = NULL;
  mesh->shared_item  = NULL;

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  adaptation information                                                    */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int
init_struct_adapt( struct hecmwST_local_mesh *mesh )
{
  HECMW_assert( mesh );

  mesh->coarse_grid_level       = 0;
  mesh->n_adapt                 = 0;
  mesh->when_i_was_refined_node = NULL;
  mesh->when_i_was_refined_elem = NULL;
  mesh->adapt_parent_type       = NULL;
  mesh->adapt_type              = NULL;
  mesh->adapt_level             = NULL;
  mesh->adapt_parent            = NULL;
  mesh->adapt_children_index    = NULL;
  mesh->adapt_children_item     = NULL;

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  refinement information                                                    */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int
init_struct_refine( struct hecmwST_local_mesh *mesh )
{
  HECMW_assert( mesh );

  mesh->n_refine     = 0;
  mesh->node_old2new = NULL;
  mesh->node_new2old = NULL;
  mesh->elem_old2new = NULL;
  mesh->elem_new2old = NULL;
  mesh->refine_origin = NULL;

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  section information                                                       */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int
init_struct_sect( struct hecmwST_local_mesh *mesh )
{
  HECMW_assert( mesh && mesh->section );

  mesh->section->n_sect            = 0;
  mesh->section->sect_type         = NULL;
  mesh->section->sect_opt          = NULL;
  mesh->section->sect_mat_ID_index = NULL;
  mesh->section->sect_mat_ID_item  = NULL;
  mesh->section->sect_I_index      = NULL;
  mesh->section->sect_I_item       = NULL;
  mesh->section->sect_R_index      = NULL;
  mesh->section->sect_R_item       = NULL;

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  material information                                                      */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int
init_struct_mat( struct hecmwST_local_mesh *mesh )
{
  HECMW_assert( mesh && mesh->material );

  mesh->material->n_mat             = 0;
  mesh->material->n_mat_item        = 0;
  mesh->material->n_mat_subitem     = 0;
  mesh->material->n_mat_table       = 0;
  mesh->material->mat_name          = NULL;
  mesh->material->mat_item_index    = NULL;
  mesh->material->mat_subitem_index = NULL;
  mesh->material->mat_table_index   = NULL;
  mesh->material->mat_val           = NULL;
  mesh->material->mat_temp          = NULL;

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  MPC information                                                           */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int
init_struct_mpc( struct hecmwST_local_mesh *mesh )
{
  HECMW_assert( mesh && mesh->mpc );

  mesh->mpc->n_mpc     = 0;
  mesh->mpc->mpc_index = NULL;
  mesh->mpc->mpc_item  = NULL;
  mesh->mpc->mpc_dof   = NULL;
  mesh->mpc->mpc_val   = NULL;
  mesh->mpc->mpc_const = NULL;

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  amplitude information                                                     */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int
init_struct_amp( struct hecmwST_local_mesh *mesh )
{
  HECMW_assert( mesh && mesh->amp );

  mesh->amp->n_amp               = 0;
  mesh->amp->amp_name            = NULL;
  mesh->amp->amp_type_definition = NULL;
  mesh->amp->amp_type_time       = NULL;
  mesh->amp->amp_type_value      = NULL;
  mesh->amp->amp_index           = NULL;
  mesh->amp->amp_val             = NULL;
  mesh->amp->amp_table           = NULL;

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  node group information                                                    */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int
init_struct_node_grp( struct hecmwST_local_mesh *mesh )
{
  HECMW_assert( mesh && mesh->node_group );

  mesh->node_group->n_grp     = 0;
  mesh->node_group->grp_name  = NULL;
  mesh->node_group->grp_index = NULL;
  mesh->node_group->grp_item  = NULL;

  mesh->node_group->n_bc         = 0;
  mesh->node_group->bc_grp_ID    = 0;
  mesh->node_group->bc_grp_type  = 0;
  mesh->node_group->bc_grp_index = 0;
  mesh->node_group->bc_grp_dof   = 0;
  mesh->node_group->bc_grp_val   = 0;

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  element group information                                                 */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int
init_struct_elem_grp( struct hecmwST_local_mesh *mesh )
{
  HECMW_assert( mesh && mesh->elem_group );

  mesh->elem_group->n_grp     = 0;
  mesh->elem_group->grp_name  = NULL;
  mesh->elem_group->grp_index = NULL;
  mesh->elem_group->grp_item  = NULL;

  mesh->elem_group->n_bc         = 0;
  mesh->elem_group->bc_grp_ID    = NULL;
  mesh->elem_group->bc_grp_type  = NULL;
  mesh->elem_group->bc_grp_index = NULL;
  mesh->elem_group->bc_grp_val   = NULL;

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  surface group information                                                 */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int
init_struct_surf_grp( struct hecmwST_local_mesh *mesh )
{
  HECMW_assert( mesh && mesh->surf_group );

  mesh->surf_group->n_grp     = 0;
  mesh->surf_group->grp_name  = NULL;
  mesh->surf_group->grp_index = NULL;
  mesh->surf_group->grp_item  = NULL;

  mesh->surf_group->n_bc         = 0;
  mesh->surf_group->bc_grp_ID    = NULL;
  mesh->surf_group->bc_grp_type  = NULL;
  mesh->surf_group->bc_grp_index = NULL;
  mesh->surf_group->bc_grp_val   = NULL;

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  contact information                                                       */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int
init_struct_contact_pair( struct hecmwST_local_mesh *mesh )
{
  HECMW_assert( mesh && mesh->contact_pair );

  mesh->contact_pair->n_pair        = 0;
  mesh->contact_pair->name     = NULL;
  mesh->contact_pair->type     = NULL;
  mesh->contact_pair->slave_grp_id  = NULL;
  mesh->contact_pair->master_grp_id = NULL;

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  initialize structures for local mesh                                      */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
int
HECMW_dist_init( struct hecmwST_local_mesh *mesh )
{
  HECMW_assert(mesh);

  /* global information */
  if( init_struct_global( mesh ) ) {
    return -1;
  }

  /* node information */
  if( init_struct_node( mesh ) ) {
    return -1;
  }

  /* element information */
  if( init_struct_elem( mesh ) ) {
    return -1;
  }

  /* domain & communication information */
  if( init_struct_comm( mesh ) ) {
    return -1;
  }

  /* adaptation information */
  if( init_struct_adapt( mesh ) ) {
    return -1;
  }

  /* refinement information */
  if( init_struct_refine( mesh ) ) {
    return -1;
  }

  /* section information */
  if( init_struct_sect( mesh ) ) {
    return -1;
  }

  /* material information */
  if( init_struct_mat( mesh ) ) {
    return -1;
  }

  /* MPC information */
  if( init_struct_mpc( mesh ) ) {
    return -1;
  }

  /* amplitude information */
  if( init_struct_amp( mesh ) ) {
    return -1;
  }

  /* node group information */
  if( init_struct_node_grp( mesh ) ) {
    return -1;
  }

  /* element group information */
  if( init_struct_elem_grp( mesh ) ) {
    return -1;
  }

  /* surface group information */
  if( init_struct_surf_grp( mesh ) ) {
    return -1;
  }

  /* contact information */
  if( init_struct_contact_pair( mesh ) ) {
    return -1;
  }

  return 0;
}

/*============================================================================*/
/*  allocate structures for local mesh                                        */
/*============================================================================*/
struct hecmwST_local_mesh *
HECMW_dist_alloc( )
{
  struct hecmwST_local_mesh *mesh;

  /* local mesh < hecmwST_local_mesh > */
  if( ( mesh = (struct hecmwST_local_mesh *)HECMW_calloc( 1, sizeof( struct hecmwST_local_mesh ) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return NULL;
  }

  /* section information < hecmwST_section > */
  if( ( mesh->section = (struct hecmwST_section *)HECMW_calloc( 1, sizeof( struct hecmwST_section ) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return NULL;
  }

  /* material information < hecmwST_material > */
  if( ( mesh->material = (struct hecmwST_material *)HECMW_calloc( 1, sizeof( struct hecmwST_material ) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return NULL;
  }

  /* MPC information < hecmwST_mpc > */
  if( ( mesh->mpc = (struct hecmwST_mpc *)HECMW_calloc( 1, sizeof( struct hecmwST_mpc ) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return NULL;
  }

  /* amplitude information < hecmwST_amplitude > */
  if( ( mesh->amp = (struct hecmwST_amplitude *)HECMW_calloc( 1, sizeof( struct hecmwST_amplitude ) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return NULL;
  }

  /* node group information < hecmwST_node_grp > */
  if( ( mesh->node_group = (struct hecmwST_node_grp *)HECMW_calloc( 1, sizeof( struct hecmwST_node_grp ) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return NULL;
  }

  /* element group information < hecmwST_elem_grp > */
  if( ( mesh->elem_group = (struct hecmwST_elem_grp *)HECMW_calloc( 1, sizeof( struct hecmwST_elem_grp ) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return NULL;
  }

  /* surface group information < hecmwST_surf_grp > */
  if( ( mesh->surf_group = (struct hecmwST_surf_grp *)HECMW_calloc( 1, sizeof( struct hecmwST_surf_grp ) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return NULL;
  }

  /* contact information < hecmwST_contact_pair > */
  if( ( mesh->contact_pair = (struct hecmwST_contact_pair *)HECMW_calloc( 1, sizeof( struct hecmwST_contact_pair ) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return NULL;
  }

  /* initialization */
  if( HECMW_dist_init( mesh ) ) {
    return NULL;
  }

  return mesh;
}
