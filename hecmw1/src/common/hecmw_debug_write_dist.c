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



#include <stdio.h>
#include "hecmw_struct.h"
#include "hecmw_util.h"
#include "hecmw_debug_write_dist.h"

/*============================================================================*/
/*                                                                            */
/*  node information                                                          */
/*                                                                            */
/*============================================================================*/
/*----------------------------------------------------------------------------*/
/*  number of nodes < n_node >                                                */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_n_node_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_node = %d", file, line, local_mesh->n_node );
}

/*----------------------------------------------------------------------------*/
/*  number of internal nodes < nn_internal >                                  */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_nn_internal_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: nn_internal = %d", file, line, local_mesh->nn_internal );
}

/*----------------------------------------------------------------------------*/
/*  maximal number of DOF < n_dof >                                           */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_n_dof_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_dof = %d", file, line, local_mesh->n_dof );
}

/*----------------------------------------------------------------------------*/
/*  number of DOF groups < n_dof_grp >                                        */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_n_dof_grp_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_dof_grp = %d", file, line, local_mesh->n_dof_grp );
}

/*----------------------------------------------------------------------------*/
/*  nodal coordinates < node >                                                */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_node_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_node = %d", file, line, local_mesh->n_node );

  if( local_mesh->node ) {
    for( i=0; i<local_mesh->n_node; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, node[3*i] = %E, node[3*i+1] = %E, node[3*i+2] = %E",
                 file, line, i,
                 local_mesh->node[3*i],
                 local_mesh->node[3*i+1],
                 local_mesh->node[3*i+2] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: node = NULL" );
  }
}

/*----------------------------------------------------------------------------*/
/*  global node id < global_node_ID >                                         */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_global_node_id_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_node = %d", file, line, local_mesh->n_node );

  if( local_mesh->global_node_ID ) {
    for( i=0; i<local_mesh->n_node; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, global_node_ID[i] = %d",
                 file, line, i,
                 local_mesh->global_node_ID[i] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: global_node_ID = NULL", file, line );
  }
}

/*----------------------------------------------------------------------------*/
/*  local node id & belonging domain of node < node_ID >                      */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_node_id_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_node = %d", file, line, local_mesh->n_node );

  if( local_mesh->node_ID ) {
    for( i=0; i<local_mesh->n_node; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, node_ID[2*i] = %d, node_ID[2*i+1] = %d",
                 file, line, i,
                 local_mesh->node_ID[2*i],
                 local_mesh->node_ID[2*i+1] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: node_ID = NULL", file, line );
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  local node id < node_ID[2*i] >                                            */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
extern void
HECMW_dbg_node_id_lid_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_node = %d", file, line, local_mesh->n_node );

  if( local_mesh->node_ID ) {
    for( i=0; i<local_mesh->n_node; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, node_ID[2*i] = %d",
                 file, line, i,
                 local_mesh->node_ID[2*i] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: node_ID = NULL", file, line );
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  belonging domain of node < node_ID[2*i+1] >                               */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
extern void
HECMW_dbg_node_id_domain_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_node = %d", file, line, local_mesh->n_node );

  if( local_mesh->node_ID ) {
    for( i=0; i<local_mesh->n_node; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, node_ID[2*i+1] = %d",
                 file, line, i,
                 local_mesh->node_ID[2*i+1] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: node_ID = NULL", file, line );
  }
}

/*----------------------------------------------------------------------------*/
/*  node DOF group < node_dof_index, node_dof_item >                          */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_node_dof_item_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_dof_grp = %d", file, line, local_mesh->n_dof_grp );

  if( local_mesh->node_dof_index && local_mesh->node_dof_item ) {
    for( i=0; i<local_mesh->n_dof_grp; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, node_dof_index[i] = %d, node_dof_index[i+1] = %d, node_dof_item[i] = %d",
                 file, line, i,
                 local_mesh->node_dof_index[i],
                 local_mesh->node_dof_index[i+1],
                 local_mesh->node_dof_item[i] );
    }
  } else {
    if( local_mesh->node_dof_index == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: node_dof_index = NULL", file, line );
    }
    if( local_mesh->node_dof_item == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: node_dof_item = NULL", file, line );
    }
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  node DOF group < node_dof_index >                                         */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
extern void
HECMW_dbg_node_dof_index_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_dof_grp = %d", file, line, local_mesh->n_dof_grp );

  if( local_mesh->node_dof_index ) {
    for( i=0; i<local_mesh->n_dof_grp; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, node_dof_index[i] = %d, node_dof_index[i+1] = %d",
                 file, line, i,
                 local_mesh->node_dof_index[i],
                 local_mesh->node_dof_index[i+1] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: node_dof_index = NULL", file, line );
  }
}

/*----------------------------------------------------------------------------*/
/*  initial condition of node < node_init_val_index, node_init_val_item >     */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_node_init_val_item_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i, j;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_node = %d", file, line, local_mesh->n_node );

  if( local_mesh->node_init_val_index && local_mesh->node_init_val_item ) {
    for( i=0; i<local_mesh->n_node; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, node_init_val_index[i] = %d, node_init_val_index[i+1] = %d",
                 file, line, i,
                 local_mesh->node_init_val_index[i],
                 local_mesh->node_init_val_index[i+1] );

      for( j=local_mesh->node_init_val_index[i]; j<local_mesh->node_init_val_index[i+1]; j++ ) {
        HECMW_log( HECMW_LOG_DEBUG, "%s:%d: j = %d, node_init_val_item[j] = %E",
                   file, line, j, local_mesh->node_init_val_item[j] );
      }
    }
  } else {
    if( local_mesh->node_init_val_index == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: node_init_val_index = NULL", file, line );
    }
    if( local_mesh->node_init_val_item == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: node_init_val_item = NULL", file, line );
    }
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  initial condition of node < node_init_val_index >                         */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
extern void
HECMW_dbg_node_init_val_index_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_node = %d", file, line, local_mesh->n_node );

  if( local_mesh->node_init_val_index ) {
    for( i=0; i<local_mesh->n_node; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, node_init_val_index[i] = %d, node_init_val_index[i+1] = %d",
                 file, line, i,
                 local_mesh->node_init_val_index[i],
                 local_mesh->node_init_val_index[i+1] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: node_init_val_index = NULL", file, line );
  }
}

/*============================================================================*/
/*                                                                            */
/*  element information                                                       */
/*                                                                            */
/*============================================================================*/
/*----------------------------------------------------------------------------*/
/*  number of elements < n_elem >                                             */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_n_elem_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_elem = %d", file, line, local_mesh->n_elem );
}

/*----------------------------------------------------------------------------*/
/*  number of internal elements < ne_internal >                               */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_ne_internal_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: ne_internal = %d", file, line, local_mesh->ne_internal );
}

/*----------------------------------------------------------------------------*/
/*  number of finite element types < n_elem_type >                            */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_n_elem_type_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_elem_type = %d", file, line, local_mesh->n_elem_type );
}

/*----------------------------------------------------------------------------*/
/*  finite element type < elem_type_index, elem_type_item >                   */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_elem_type_item_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_elem_type = %d", file, line, local_mesh->n_elem_type );

  if( local_mesh->elem_type_index && local_mesh->elem_type_item ) {
    for( i=0; i<local_mesh->n_elem_type; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, elem_type_index[i] = %d, elem_type_index[i+1] = %d: elem_type_item[i] = %d",
                 file, line, i,
                 local_mesh->elem_type_index[i],
                 local_mesh->elem_type_index[i+1],
                 local_mesh->elem_type_item[i] );
    }
  } else {
    if( local_mesh->elem_type_index == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: elem_type_index = NULL", file, line );
    }
    if( local_mesh->elem_type_item == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: elem_type_item = NULL", file, line );
    }
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  finite element type < elem_type_index >                                   */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
extern void
HECMW_dbg_elem_type_index_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_elem_type = %d", file, line, local_mesh->n_elem_type );

  if( local_mesh->elem_type_index ) {
    for( i=0; i<local_mesh->n_elem_type; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, elem_type_index[i] = %d, elem_type_index[i+1] = %d",
                 file, line, i,
                 local_mesh->elem_type_index[i],
                 local_mesh->elem_type_index[i+1] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: elem_type_index = NULL", file, line );
  }
}

/*----------------------------------------------------------------------------*/
/*  finite element type < elem_type >                                         */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_elem_type_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_elem = %d", file, line, local_mesh->n_elem );

  if( local_mesh->elem_type ) {
    for( i=0; i<local_mesh->n_elem; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, elem_type[i] = %d",
                 file, line, i,
                 local_mesh->elem_type[i] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: elem_type = NULL", file, line );
  }
}

/*----------------------------------------------------------------------------*/
/*  component node of element < elem_node_index, elem_node_item >             */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_elem_node_item_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i, j;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_elem = %d", file, line, local_mesh->n_elem );

  if( local_mesh->elem_node_index && local_mesh->elem_node_item ) {
    for( i=0; i<local_mesh->n_elem; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, elem_node_index[i] = %d, elem_node_index[i+1] = %d",
                 file, line, i,
                 local_mesh->elem_node_index[i],
                 local_mesh->elem_node_index[i+1] );

      for( j=local_mesh->elem_node_index[i]; j<local_mesh->elem_node_index[i+1]; j++ ) {
        HECMW_log( HECMW_LOG_DEBUG, "%s:%d: j = %d, elem_node_item[j] = %d",
                   file, line, j,
                   local_mesh->elem_node_item[j] );
      }
    }
  } else {
    if( local_mesh->elem_node_index == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: elem_node_index = NULL", file, line );
    }
    if( local_mesh->elem_node_item == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: elem_node_item = NULL", file, line );
    }
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  component node of element < elem_node_index >                             */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
extern void
HECMW_dbg_elem_node_index_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_elem = %d", file, line, local_mesh->n_elem );

  if( local_mesh->elem_node_index ) {
    for( i=0; i<local_mesh->n_elem; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, elem_node_index[i] = %d, elem_node_index[i+1] = %d",
                 file, line, i,
                 local_mesh->elem_node_index[i],
                 local_mesh->elem_node_index[i+1] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: elem_node_index = NULL", file, line );
  }
}

/*----------------------------------------------------------------------------*/
/*  local element id & belonging domain of element < elem_ID >                */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_elem_id_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_elem = %d", file, line, local_mesh->n_elem );

  if( local_mesh->elem_ID ) {
    for( i=0; i<local_mesh->n_elem; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, elem_ID[2*i] = %d, elem_ID[2*i+1] = %d",
                 file, line, i,
                 local_mesh->elem_ID[2*i],
                 local_mesh->elem_ID[2*i+1] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: elem_ID = NULL", file, line );
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  local element id < elem_ID[2*i] >                                         */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
extern void
HECMW_dbg_elem_id_lid_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_elem = %d", file, line, local_mesh->n_elem );

  if( local_mesh->elem_ID ) {
    for( i=0; i<local_mesh->n_elem; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, elem_ID[2*i] = %d",
                 file, line, i,
                 local_mesh->elem_ID[2*i] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: elem_ID = NULL", file, line );
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  belonging domain of element < elem_ID[2*i+1] >                            */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
extern void
HECMW_dbg_elem_id_domain_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_elem = %d", file, line, local_mesh->n_elem );

  if( local_mesh->elem_ID ) {
    for( i=0; i<local_mesh->n_elem; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, elem_ID[2*i+1] = %d",
                 file, line, i,
                 local_mesh->elem_ID[2*i+1] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: elem_ID = NULL", file, line );
  }
}

/*----------------------------------------------------------------------------*/
/*  global element id < global_element_ID >                                   */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_global_elem_id_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_elem = %d", file, line, local_mesh->n_elem );

  if( local_mesh->global_elem_ID ) {
    for( i=0; i<local_mesh->n_elem; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, global_elem_ID[i] = %d",
                 file, line, i,
                 local_mesh->global_elem_ID[i] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: global_elem_ID = NULL", file, line );
  }
}

/*----------------------------------------------------------------------------*/
/*  list of internal element < elem_internal_list >                           */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_elem_internal_list_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: ne_internal = %d", file, line, local_mesh->ne_internal );

  if( local_mesh->elem_internal_list ) {
    for( i=0; i<local_mesh->ne_internal; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, elem_internal_list[i] = %d",
                 file, line, i,
                 local_mesh->elem_internal_list[i] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: elem_internal_list = NULL", file, line );
  }
}

/*----------------------------------------------------------------------------*/
/*  section id < section_ID >                                                 */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_section_id_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_elem = %d", file, line, local_mesh->n_elem );

  if( local_mesh->section_ID ) {
    for( i=0; i<local_mesh->n_elem; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, section_ID[i] = %d",
                 file, line, i,
                 local_mesh->section_ID[i] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: section_ID = NULL", file, line );
  }
}

/*----------------------------------------------------------------------------*/
/*  number of material id < n_elem_mat_ID >                                   */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_n_elem_mat_id_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_elem_mat_ID = %d", file, line, local_mesh->n_elem_mat_ID );
}

/*----------------------------------------------------------------------------*/
/*  material id < elem_mat_ID_index, elem_mat_ID_item >                       */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_elem_mat_id_item_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i, j;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_elem = %d", file, line, local_mesh->n_elem );

  if( local_mesh->elem_mat_ID_index && local_mesh->elem_mat_ID_item ) {
    for( i=0; i<local_mesh->n_elem; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, elem_mat_ID_index[i] = %d, elem_mat_ID_index[i+1] = %d",
                 file, line, i,
                 local_mesh->elem_mat_ID_index[i],
                 local_mesh->elem_mat_ID_index[i+1] );

      for( j=local_mesh->elem_mat_ID_index[i]; j<local_mesh->elem_mat_ID_index[i+1]; j++ ) {
        HECMW_log( HECMW_LOG_DEBUG, "%s:%d: j = %d, elem_mat_ID_item[j] = %d",
                   file, line, j,
                   local_mesh->elem_mat_ID_item[j] );
      }
    }
  } else {
    if( local_mesh->elem_mat_ID_index == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: elem_mat_ID_index = NULL", file, line );
    }
    if( local_mesh->elem_mat_ID_item == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: elem_mat_ID_item = NULL", file, line );
    }
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  material id < elem_mat_ID_index >                                         */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
extern void
HECMW_dbg_elem_mat_id_index_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_elem = %d", file, line, local_mesh->n_elem );

  if( local_mesh->elem_mat_ID_index ) {
    for( i=0; i<local_mesh->n_elem; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, elem_mat_ID_index[i] = %d, elem_mat_ID_index[i+1] = %d",
                 file, line, i,
                 local_mesh->elem_mat_ID_index[i],
                 local_mesh->elem_mat_ID_index[i+1] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: elem_mat_ID_index = NULL", file, line );
  }
}

/*============================================================================*/
/*                                                                            */
/*  parallel & communication table information                                */
/*                                                                            */
/*============================================================================*/
/*----------------------------------------------------------------------------*/
/*  number of neighboring domains < n_neighbor_pe >                           */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_n_neighbor_pe_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_neighbor_pe = %d", file, line, local_mesh->n_neighbor_pe );
}

/*----------------------------------------------------------------------------*/
/*  neighbor domain < neighbor_pe >                                           */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_neighbor_pe_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_neighbor_pe = %d", file, line, local_mesh->n_neighbor_pe );

  if( local_mesh->neighbor_pe ) {
    for( i=0; i<local_mesh->n_neighbor_pe; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, neighbor_pe[i] = %d",
                 file, line, i,
                 local_mesh->neighbor_pe[i] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: neighbor_pe = NULL", file, line );
  }
}

/*----------------------------------------------------------------------------*/
/*  import information < import_index, import_item >                          */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_import_item_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i, j;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_neighbor_pe = %d", file, line, local_mesh->n_neighbor_pe );

  if( local_mesh->import_index && local_mesh->import_item ) {
    for( i=0; i<local_mesh->n_neighbor_pe; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, neighbor_pe[i] = %d, import_index[i] = %d, import_index[i+1] = %d",
                 file, line, i,
                 local_mesh->neighbor_pe[i],
                 local_mesh->import_index[i],
                 local_mesh->import_index[i+1] );

      for( j=local_mesh->import_index[i]; j<local_mesh->import_index[i+1]; j++ ) {
        HECMW_log( HECMW_LOG_DEBUG, "%s:%d: j = %d, import_item[j] = %d",
                   file, line, j,
                   local_mesh->import_item[j] );
      }
    }
  } else {
    if( local_mesh->import_index == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: import_index = NULL", file, line );
    }
    if( local_mesh->import_item == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: import_item = NULL", file, line );
    }
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  import information < import_index >                                       */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
extern void
HECMW_dbg_import_index_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_neighbor_pe = %d", file, line, local_mesh->n_neighbor_pe );

  if( local_mesh->import_index ) {
    for( i=0; i<local_mesh->n_neighbor_pe; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, neighbor_pe[i] = %d, import_index[i] = %d, import_index[i+1] = %d",
                 file, line, i,
                 local_mesh->neighbor_pe[i],
                 local_mesh->import_index[i],
                 local_mesh->import_index[i+1] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: import_index = NULL", file, line );
  }
}

/*----------------------------------------------------------------------------*/
/*  export information < export_index, export_item >                          */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_export_item_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i, j;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_neighbor_pe = %d", file, line, local_mesh->n_neighbor_pe );

  if( local_mesh->export_index && local_mesh->export_item ) {
    for( i=0; i<local_mesh->n_neighbor_pe; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, neighbor_pe[i] = %d, export_index[i] = %d, export_index[i+1] = %d",
                 file, line, i,
                 local_mesh->neighbor_pe[i],
                 local_mesh->export_index[i],
                 local_mesh->export_index[i+1] );

      for( j=local_mesh->export_index[i]; j<local_mesh->export_index[i+1]; j++ ) {
        HECMW_log( HECMW_LOG_DEBUG, "%s:%d: j = %d, export_item[j] = %d",
                   file, line, j,
                   local_mesh->export_item[j] );
      }
    }
  } else {
    if( local_mesh->export_index == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: export_index = NULL", file, line );
    }
    if( local_mesh->export_item == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: export_item = NULL", file, line );
    }
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  export information < export_index >                                       */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
extern void
HECMW_dbg_export_index_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_neighbor_pe = %d", file, line, local_mesh->n_neighbor_pe );

  if( local_mesh->export_index ) {
    for( i=0; i<local_mesh->n_neighbor_pe; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, neighbor_pe[i] = %d, export_index[i] = %d, export_index[i+1] = %d",
                 file, line, i,
                 local_mesh->neighbor_pe[i],
                 local_mesh->export_index[i],
                 local_mesh->export_index[i+1] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: export_index = NULL", file, line );
  }
}

/*----------------------------------------------------------------------------*/
/*  shared information < shared_index, export_item >                          */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_shared_item_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i, j;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_neighbor_pe = %d", file, line, local_mesh->n_neighbor_pe );

  if( local_mesh->shared_index && local_mesh->shared_item ) {
    for( i=0; i<local_mesh->n_neighbor_pe; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, neighbor_pe[i] = %d, shared_index[i] = %d, shared_index[i+1] = %d",
                 file, line, i,
                 local_mesh->neighbor_pe[i],
                 local_mesh->shared_index[i],
                 local_mesh->shared_index[i+1] );

      for( j=local_mesh->shared_index[i]; j<local_mesh->shared_index[i+1]; j++ ) {
        HECMW_log( HECMW_LOG_DEBUG, "%s:%d: j = %d, shared_item[j] = %d",
                   file, line, j,
                   local_mesh->shared_item[j] );
      }
    }
  } else {
    if( local_mesh->shared_index == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: shared_index = NULL", file, line );
    }
    if( local_mesh->shared_item == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: shared_item = NULL", file, line );
    }
  }
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  shared information < shared_index >                                       */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
extern void
HECMW_dbg_shared_index_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: n_neighbor_pe = %d", file, line, local_mesh->n_neighbor_pe );

  if( local_mesh->shared_index ) {
    for( i=0; i<local_mesh->n_neighbor_pe; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, neighbor_pe[i] = %d, shared_index[i] = %d, shared_index[i+1] = %d",
                 file, line, i,
                 local_mesh->neighbor_pe[i],
                 local_mesh->shared_index[i],
                 local_mesh->shared_index[i+1] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: shared_index = NULL", file, line );
  }
}

/*============================================================================*/
/*                                                                            */
/*  section information                                                       */
/*                                                                            */
/*============================================================================*/
/*----------------------------------------------------------------------------*/
/*  number of sections < section->n_sect >                                    */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_n_sect_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: section->n_sect = %d", file, line, local_mesh->section->n_sect );
}

/*----------------------------------------------------------------------------*/
/*  section information                                                       */
/*  < section->n_sect, section->sect_type, section->sect_opt,                 */
/*    section->sect_mat_ID_index, section->sect_mat_ID_item,                  */
/*    section->sect_I_index, section->sect_I_item,                            */
/*    section->sect_R_index, section->sect_R_item >                           */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_section_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i, j;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: section->n_sect = %d", file, line, local_mesh->section->n_sect );

  if( local_mesh->section->sect_type && local_mesh->section->sect_opt ) {
    for( i=0; i<local_mesh->section->n_sect; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, section->sect_type[i] = %d: section->sect_opt[i] = %d",
                 file, line, i,
                 local_mesh->section->sect_type[i],
                 local_mesh->section->sect_opt[i] );
    }
  } else {
    if( local_mesh->section->sect_type == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: section->sect_type = NULL", file, line );
    }
    if( local_mesh->section->sect_opt == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: section->sect_opt = NULL", file, line );
    }
  }

  if( local_mesh->section->sect_mat_ID_index && local_mesh->section->sect_mat_ID_item ) {
    for( i=0; i<local_mesh->section->n_sect; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, section->sect_mat_ID_index[i] = %d, section->sect_mat_ID_index[i+1] = %d",
                 file, line, i,
                 local_mesh->section->sect_mat_ID_index[i],
                 local_mesh->section->sect_mat_ID_index[i+1] );

      for( j=local_mesh->section->sect_mat_ID_index[i]; j<local_mesh->section->sect_mat_ID_index[i+1]; j++ ) {
        HECMW_log( HECMW_LOG_DEBUG, "%s:%d: j = %d, section->sect_mat_ID_item[j] = %d",
                   file, line, j,
                   local_mesh->section->sect_mat_ID_item[j] );
      }
    }
  } else {
    if( local_mesh->section->sect_mat_ID_index == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: section->sect_mat_ID_index = NULL", file, line );
    }
    if( local_mesh->section->sect_mat_ID_item == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: section->sect_mat_ID_item = NULL", file, line );
    }
  }

  if( local_mesh->section->sect_I_index && local_mesh->section->sect_I_item ) {
    for( i=0; i<local_mesh->section->n_sect; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, section->sect_I_index[i] = %d, section->sect_I_index[i+1] = %d",
                 file, line, i,
                 local_mesh->section->sect_I_index[i],
                 local_mesh->section->sect_I_index[i+1] );

      for( j=local_mesh->section->sect_I_index[i]; j<local_mesh->section->sect_I_index[i+1]; j++ ){
        HECMW_log( HECMW_LOG_DEBUG, "%s:%d: j = %d, section->sect_I_item[j] = %d",
                   file, line, j,
                   local_mesh->section->sect_I_item[j] );
      }
    }
  } else {
    if( local_mesh->section->sect_I_index == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: section->sect_I_index = NULL", file, line );
    }
    if( local_mesh->section->sect_I_item == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: section->sect_I_item = NULL", file, line );
    }
  }

  if( local_mesh->section->sect_R_index && local_mesh->section->sect_R_item ) {
    for( i=0; i<local_mesh->section->n_sect; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, section->sect_R_index[i] = %d, section->sect_R_index[i+1] = %d",
                 file, line, i,
                 local_mesh->section->sect_R_index[i],
                 local_mesh->section->sect_R_index[i+1] );

      for( j=local_mesh->section->sect_R_index[i]; j<local_mesh->section->sect_R_index[i+1]; j++ ) {
        HECMW_log( HECMW_LOG_DEBUG, "%s:%d: j = %d, section->sect_R_item[j] = %E",
                   file, line, j,
                   local_mesh->section->sect_R_item[j] );
      }
    }
  } else {
    if( local_mesh->section->sect_R_index == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: section->sect_R_index = NULL", file, line );
    }
    if( local_mesh->section->sect_R_item == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: section->sect_R_item = NULL", file, line );
    }
  }
}

/*============================================================================*/
/*                                                                            */
/*  material information                                                      */
/*                                                                            */
/*============================================================================*/
/*----------------------------------------------------------------------------*/
/*  number of materials < material->n_mat >                                   */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_n_mat_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: material->n_mat = %d", file, line, local_mesh->material->n_mat );
}


/*----------------------------------------------------------------------------*/
/*  material information                                                      */
/*  < material->n_mat, material->n_mat_item, material->n_mat_subitem,         */
/*    material->n_mat_table, material->mat_name, material->mat_item_index,    */                                  
/*    material->mat_subitem_index, material->mat_table_index,                 */
/*    material->mat_val, material->mat_temp >                                 */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_material_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i, j;

  /* material->n_mat */
  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: material->n_mat = %d", file, line, local_mesh->material->n_mat );

  /* material->mat_item_index */
  if( local_mesh->material->mat_item_index ) {
    for( i=0; i<local_mesh->material->n_mat; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, material->mat_item_index[i] = %d, material->mat_item_index[i+1] = %d",
                 file, line, i,
                 local_mesh->material->mat_item_index[i],
                 local_mesh->material->mat_item_index[i+1] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: material->mat_item_index = NULL", file, line );
  }

  /* material->mat_name */
  if( local_mesh->material->mat_name ) {
    for( i=0; i<local_mesh->material->n_mat; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, material->mat_name[i] = \"%s\"",
                 file, line, i,
                 local_mesh->material->mat_name[i] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: material->mat_name = NULL", file, line );
  }

  /* material->n_mat_item */
  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: material->n_mat_item = %d", file, line, local_mesh->material->n_mat_item );

  /* material->mat_subitem_index */
  if( local_mesh->material->mat_subitem_index ) {
    for( i=0; i<local_mesh->material->n_mat_item; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, material->mat_subitem_index[i] = %d, material->mat_subitem_index[i+1] = %d",
                 file, line, i,
                 local_mesh->material->mat_subitem_index[i],
                 local_mesh->material->mat_subitem_index[i+1] );
    }
  } else {
    HECMW_log( HECMW_LOG_DEBUG, "%s:%d: material->mat_subitem_index = NULL", file, line );
  }

  /* material->n_mat_subitem */
  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: material->n_mat_subitem = %d", file, line, local_mesh->material->n_mat_subitem );

  /* material->mat_table_index, material->mat_val, material->mat_temp */
  if( local_mesh->material->mat_table_index && local_mesh->material->mat_val && local_mesh->material->mat_temp ) {
    for( i=0; i<local_mesh->material->n_mat_subitem; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, material->mat_table_index[i] = %d, material->mat_table_index[i+1] = %d",
                 file, line, i,
                 local_mesh->material->mat_table_index[i],
                 local_mesh->material->mat_table_index[i+1] );
  
      for( j=local_mesh->material->mat_table_index[i]; j<local_mesh->material->mat_table_index[i+1]; j++ ) {
        HECMW_log( HECMW_LOG_DEBUG, "%s:%d: j = %d, material->mat_val[j] = %E: material->mat_temp[j] = %E",
                   file, line, j,
                   local_mesh->material->mat_val[j],
                   local_mesh->material->mat_temp[j] );
      }
    }
  } else {
    if( local_mesh->material->mat_table_index == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: material->mat_table_index = NULL", file, line );
    }
    if( local_mesh->material->mat_val == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: material->mat_val = NULL", file, line );
    }
    if( local_mesh->material->mat_temp == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: material->mat_temp = NULL", file, line );
    }
  }

  /* material->n_mat_table */
  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: material->n_mat_table = %d", file, line, local_mesh->material->n_mat_table );
}

/*============================================================================*/
/*                                                                            */
/*  MPC group information                                                     */
/*                                                                            */
/*============================================================================*/
/*----------------------------------------------------------------------------*/
/*  number of MPC groups < mpc->n_mpc >                                       */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_n_mpc_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: mpc->n_mpc = %d", file, line, local_mesh->mpc->n_mpc );
}


/*----------------------------------------------------------------------------*/
/*  MPC group information                                                     */
/*  < mpc->n_mpc, mpc->mpc_index,                                             */
/*    mpc->mpc_item, mpc->mpc_dof, mpc->mpc_val >                             */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_mpc_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i, j;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: mpc->n_mpc = %d", file, line, local_mesh->mpc->n_mpc );

  if( local_mesh->mpc->mpc_index && local_mesh->mpc->mpc_item &&
      local_mesh->mpc->mpc_dof && local_mesh->mpc->mpc_val ) {
    for( i=0; i<local_mesh->mpc->n_mpc; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, mpc->mpc_index[i] = %d, mpc->mpc_index[i+1] = %d",
                 file, line, i,
                 local_mesh->mpc->mpc_index[i],
                 local_mesh->mpc->mpc_index[i+1] );

      for( j=local_mesh->mpc->mpc_index[i]; j<local_mesh->mpc->mpc_index[i+1]; j++ ) {
        HECMW_log( HECMW_LOG_DEBUG, "%s:%d: j = %d, mpc->mpc_item[j] = %d: mpc->mpc_dof[j] = %d: mpc->mpc_val[j] = %E",
                   file, line, j,
                   local_mesh->mpc->mpc_item[j],
                   local_mesh->mpc->mpc_dof[j],
                   local_mesh->mpc->mpc_val[j] );
      }
    }
  } else {
    if( local_mesh->mpc->mpc_index == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: mpc->mpc_index = NULL", file, line );
    }
    if( local_mesh->mpc->mpc_item == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: mpc->mpc_item = NULL", file, line );
    }
    if( local_mesh->mpc->mpc_dof == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: mpc->mpc_dof = NULL", file, line );
    }
    if( local_mesh->mpc->mpc_val == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: mpc->mpc_val = NULL", file, line );
    }
  }
}

/*============================================================================*/
/*                                                                            */
/*  amplitude information                                                     */
/*                                                                            */
/*============================================================================*/
/*----------------------------------------------------------------------------*/
/*  number of amplitude groups < amp->n_amp >                                 */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_n_amp_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: amp->n_amp = %d", file, line, local_mesh->amp->n_amp );
}

/*----------------------------------------------------------------------------*/
/*  amplitude group information                                               */
/*  < amp->n_amp, amp->amp_index,                                             */
/*    amp->amp_type_definition, amp->amp_type_time, amp->amp_type_value,      */
/*    amp->amp_val, amp->amp_table >                                          */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_amp_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i, j;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: amp->n_amp = %d", file, line, local_mesh->amp->n_amp );

  if( local_mesh->amp->amp_index     && local_mesh->amp->amp_type_definition &&
      local_mesh->amp->amp_type_time && local_mesh->amp->amp_type_value      &&
      local_mesh->amp->amp_val       && local_mesh->amp->amp_table ) {
    for( i=0; i<local_mesh->amp->n_amp; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, amp->amp_index[i] = %d, amp->amp_index[i+1] = %d",
                 file, line, i,
                 local_mesh->amp->amp_index[i],
                 local_mesh->amp->amp_index[i+1] );

      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, amp->amp_type_definition[i] = %d: amp->amp_type_time[i] = %d: amp->amp_type_value[i] = %d",
                 file, line, i,
                 local_mesh->amp->amp_type_definition[i],
                 local_mesh->amp->amp_type_time[i],
                 local_mesh->amp->amp_type_value[i] );

      for( j=local_mesh->amp->amp_index[i]; j<local_mesh->amp->amp_index[i+1]; j++ ) {
        HECMW_log( HECMW_LOG_DEBUG, "%s:%d: j = %d, amp->amp_val[j] = %E, amp_table[j] = %E",
                   file, line, j,
                   local_mesh->amp->amp_val[j],
                   local_mesh->amp->amp_table[j] );
      }
    }
  } else {
    if( local_mesh->amp->amp_index == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: amp->amp_index = NULL", file, line );
    }
    if( local_mesh->amp->amp_type_definition == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: amp->amp_type_definition = NULL", file, line );
    }
    if( local_mesh->amp->amp_type_time == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: amp->amp_type_time = NULL", file, line );
    }
    if( local_mesh->amp->amp_type_value == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: amp->amp_type_value = NULL", file, line );
    }
    if( local_mesh->amp->amp_val == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: amp->amp_val = NULL", file, line );
    }
    if( local_mesh->amp->amp_table == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: amp->amp_table = NULL", file, line );
    }
  }
}

/*============================================================================*/
/*                                                                            */
/*  node group information                                                    */
/*                                                                            */
/*============================================================================*/
/*----------------------------------------------------------------------------*/
/*  number of node groups < node_group->n_grp >                               */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_n_node_grp_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: node_group->n_grp = %d", file, line, local_mesh->node_group->n_grp );
}

/*----------------------------------------------------------------------------*/
/*  node group information                                                    */
/*  < node_group->n_grp, node_group->grp_index,                               */
/*    node_group->grp_name, node_group->grp_item >                            */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_node_group_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i, j;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: node_group->n_grp = %d", file, line, local_mesh->node_group->n_grp );

  if( local_mesh->node_group->grp_index &&
      local_mesh->node_group->grp_name  &&
      local_mesh->node_group->grp_item ) {
    for( i=0; i<local_mesh->node_group->n_grp; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, node_group->grp_index[i] = %d, node_group->grp_index[i+1] = %d",
                 file, line, i,
                 local_mesh->node_group->grp_index[i],
                 local_mesh->node_group->grp_index[i+1] );

      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, node_group->grp_name[i] = \"%s\"",
                 file, line, i,
                 local_mesh->node_group->grp_name[i] );

      for( j=local_mesh->node_group->grp_index[i]; j<local_mesh->node_group->grp_index[i+1]; j++ ) {
        HECMW_log( HECMW_LOG_DEBUG, "%s:%d: j = %d, node_group->grp_item[j] = %d",
                   file, line, j,
                   local_mesh->node_group->grp_item[j] );
      }
    }
  } else {
    if( local_mesh->node_group->grp_index == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: node_group->grp_index == NULL", file, line );
    }
    if( local_mesh->node_group->grp_name == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: node_group->grp_name == NULL", file, line );
    }
    if( local_mesh->node_group->grp_item == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: node_group->grp_item == NULL", file, line );
    }
  }
}

/*============================================================================*/
/*                                                                            */
/*  element group information                                                 */
/*                                                                            */
/*============================================================================*/
/*----------------------------------------------------------------------------*/
/*  number of element groups < elem_group->n_grp >                            */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_n_elem_grp_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: elem_group->n_grp = %d", file, line, local_mesh->elem_group->n_grp );
}

/*----------------------------------------------------------------------------*/
/*  element group information                                                 */
/*  < elem_group->n_grp, elem_group->grp_index,                               */
/*    elem_group->grp_name, elem_group->grp_item >                            */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_elem_group_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i, j;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: elem_group->n_grp = %d", file, line, local_mesh->elem_group->n_grp );

  if( local_mesh->elem_group->grp_index &&
      local_mesh->elem_group->grp_name  &&
      local_mesh->elem_group->grp_item ) {
    for( i=0; i<local_mesh->elem_group->n_grp; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, elem_group->grp_index[i] = %d, elem_group->grp_index[i+1] = %d",
                 file, line, i,
                 local_mesh->elem_group->grp_index[i],
                 local_mesh->elem_group->grp_index[i+1] );

      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, elem_group->grp_name[i] = \"%s\"",
                 file, line, i,
                 local_mesh->elem_group->grp_name[i] );

      for( j=local_mesh->elem_group->grp_index[i]; j<local_mesh->elem_group->grp_index[i+1]; j++ ) {
        HECMW_log( HECMW_LOG_DEBUG, "%s:%d: j = %d, elem_group->grp_item[j] = %d",
                   file, line, j,
                   local_mesh->elem_group->grp_item[j] );
      }
    }
  } else {
    if( local_mesh->elem_group->grp_index == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: elem_group->grp_index = NULL", file, line );
    }
    if( local_mesh->elem_group->grp_name == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: elem_group->grp_name = NULL", file, line );
    }
    if( local_mesh->elem_group->grp_item == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: elem_group->grp_item = NULL", file, line );
    }
  }
}

/*============================================================================*/
/*                                                                            */
/*  surface group information                                                 */
/*                                                                            */
/*============================================================================*/
/*----------------------------------------------------------------------------*/
/*  number of surface groups < surf_group->n_grp >                            */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_n_surf_grp_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: surf_group->n_grp = %d", file, line, local_mesh->surf_group->n_grp );
}

/*----------------------------------------------------------------------------*/
/*  surface group information                                                 */
/*  < surf_group->n_grp, surf_group->grp_index,                               */
/*    surf_group->grp_name, surf_group->grp_item >                            */
/*----------------------------------------------------------------------------*/
extern void
HECMW_dbg_surf_group_( struct hecmwST_local_mesh *local_mesh, char *file, int line )
{
  int i, j;

  HECMW_log( HECMW_LOG_DEBUG, "%s:%d: surf_group->n_grp = %d", file, line, local_mesh->surf_group->n_grp );

  if( local_mesh->surf_group->grp_index &&
      local_mesh->surf_group->grp_name  &&
      local_mesh->surf_group->grp_item ) {
    for( i=0; i<local_mesh->surf_group->n_grp; i++ ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, surf_group->grp_index[i] = %d, surf_group->grp_index[i+1] = %d",
                 file, line, i,
                 local_mesh->surf_group->grp_index[i],
                 local_mesh->surf_group->grp_index[i+1] );

      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: i = %d, surf_group->grp_name[i] = \"%s\"",
                 file, line, i,
                 local_mesh->surf_group->grp_name[i] );

      for( j=local_mesh->surf_group->grp_index[i]; j<local_mesh->surf_group->grp_index[i+1]; j++ ) {
        HECMW_log( HECMW_LOG_DEBUG, "%s:%d: j = %d, surf_group->grp_item[2*i] = %d, surf_group->grp_item[2*i+1] = %d",
                   file, line, j,
                   local_mesh->surf_group->grp_item[2*j],
                   local_mesh->surf_group->grp_item[2*j+1] );
      }
    }
  } else {
    if( local_mesh->surf_group->grp_index == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: surf_group->grp_index = NULL", file, line );
    }
    if( local_mesh->surf_group->grp_name == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: surf_group->grp_name = NULL", file, line );
    }
    if( local_mesh->surf_group->grp_item == NULL ) {
      HECMW_log( HECMW_LOG_DEBUG, "%s:%d: surf_group->grp_item = NULL", file, line );
    }
  }
}
