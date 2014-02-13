/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
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

#include "hecmw_msgno.h"
#include "hecmw_malloc.h"
#include "hecmw_error.h"

#include "hecmw_part_define.h"
#include "hecmw_part_get_control.h"
#include "hecmw_init_for_partition.h"

#define DEFAULT_CONTROL_FILE_NAME "hecmw_part_ctrl.dat"


static void
print_usage( void )
{
    fprintf( stderr, "Usage: hecmw_part [-f file]\n" );
    fprintf( stderr, "\n" );
    fprintf( stderr, "  -f  specify control file name\n" );
    fprintf( stderr, "  -h  print usage\n" );
}


extern int
HECMW_init_for_partition( int argc, char **argv )
{
    int counter;
    char control_file_name[HECMW_FILENAME_LEN+1];

    
    if( argc == 1 ) {
        strcpy( control_file_name, DEFAULT_CONTROL_FILE_NAME );
        if( HECMW_part_set_ctrl_file_name( control_file_name ) )  goto error;

    
    } else {
        counter=1;
        while( counter < argc ) {

            
            if( !strcmp( argv[counter], "-f" ) ) {
                counter++;
                if( counter >= argc ) {
                    print_usage( );
                    goto error;
                }
                if( strlen( argv[counter] ) > HECMW_FILENAME_LEN ) {
                    HECMW_set_error( HECMW_PART_E_TOO_LONG_FNAME,
                                     "%s", "control file for partitioner" );
                    goto error;
                }
                strcpy( control_file_name, argv[counter] );
                if( HECMW_part_set_ctrl_file_name( control_file_name ) )  goto error;
                counter++;

            
            } else {
                print_usage( );
                goto error;
            }
        }
    }

    return 0;

error:
    return -1;
}
