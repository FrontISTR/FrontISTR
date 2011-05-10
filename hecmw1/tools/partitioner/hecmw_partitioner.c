/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.1                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : HEC-MW Utility                                    *
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
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>

#include "hecmw_util.h"
#include "hecmw_io.h"
#include "hecmw_init_for_partition.h"
#include "hecmw_partition.h"


int
main( int argc, char **argv )
{
  struct hecmwST_local_mesh *global_mesh = NULL;
  struct hecmwST_local_mesh *local_mesh = NULL;
  int rtc;

  
  rtc = HECMW_init( &argc, &argv );
  if( rtc != 0 )  goto error;

  /* set log information */
/*  HECMW_setloglv( HECMW_LOG_DEBUG ); */
/*  HECMW_setloglv( HECMW_LOG_INFO ); */
/*  HECMW_openlog( "partition_debug.log", HECMW_LOG_ALL, 1 ); */

  
  if( HECMW_comm_get_rank( ) == 0 ) {

    
    rtc = HECMW_init_for_partition( argc, argv );
    if( rtc != 0 )  goto error;

    
    HECMW_log( HECMW_LOG_INFO, "Reading mesh file..." );
    global_mesh = HECMW_get_mesh( "part_in" );
    if( global_mesh == NULL )  goto error;

    
    local_mesh = HECMW_partition( global_mesh );
    if( local_mesh == NULL )  goto error;

    
    HECMW_dist_free( global_mesh );

  }

  
/*  HECMW_closelog( "partition_debug.log" ); */

  
  HECMW_finalize( );

  return 0;

error:
  HECMW_dist_free( global_mesh );
  HECMW_finalize( );
  HECMW_abort( HECMW_comm_get_comm( ) );

  return -1;
}
