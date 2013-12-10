/*=====================================================================*
!                                                                      !
! Software Name : FrontISTR Ver. 4.2                                   !
!                                                                      !
!      Module Name : Starter Main                                      !
!                                                                      !
!            Written by Keiji Suemitsu (Advancesoft)                   !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
*=====================================================================*/

#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define BUFFSIZE 256

int main( int argc , char **argv )
{
	char path[BUFFSIZE], filename[BUFFSIZE], ctrlfile[BUFFSIZE], buff[BUFFSIZE];
	char *str;
	int ver;
	FILE *fp;

	if( argc > 2 ) {
		fprintf( stderr, "Error : incorrect command line parameter\n" );
		return -1;
	}

	if( argc > 1 ) {
		strcpy( filename, argv[1] );
	} else {
		strcpy( filename, "hecmw_ctrl.dat" );
	}

	if( ( fp = fopen( filename, "r" ) ) == 0 ) {
		fprintf( stderr, "Error : fail at opening total control file\n" );
		return -1;
	}
	while( fgets( buff, BUFFSIZE, fp ) ) {
		if( strstr( buff, "!CONTROL" ) ) {
			fgets( buff, BUFFSIZE, fp );
			buff[ strlen(buff)-1 ] = 0;
			if( strstr( buff, "\r" ) ) buff[ strlen(buff)-1 ] = 0;
			str = strrchr( buff, ' ' );
			strcpy( ctrlfile, str+1 );
			break;
		}
	}
	if( fclose( fp ) ) {
		fprintf( stderr, "Error : fail at closing total control file\n" );
		return -1;
	}

	if( ( fp = fopen( ctrlfile, "r" ) ) == 0 ) {
		fprintf( stderr, "Error : fail at opening analysis control file\n" );
		return -1;
	}
	ver = 0;
	while( fgets( buff, BUFFSIZE, fp ) ) {
		if( strstr( buff, "!SOLUTION" ) ) {
			if( !strstr( buff, "STATIC" ) && !strstr( buff, "NLSTATIC" ) ) {
				ver = 3;
				break;
			}
		} else if( strstr( buff, "!VISCOELASTIC" ) || strstr( buff, "!CREEP" ) ) {
				ver = 3;
				break;
		} else if( strstr( buff, "!CONTACT" ) ) {
			ver = 3;
			break;
		}
	}
	if( ver == 0 ) ver = 4;
	if( fclose( fp ) ) {
		fprintf( stderr, "Error : fail at closing analysis control file\n" );
		return -1;
	}

	if( ver == 3 ) {
		sprintf( path, "%s1", argv[0] );
		if( argc > 1 ) {
			execl( path, path, argv[1], NULL );
		} else {
			execl( path, path, NULL );
		}
	} else {
		sprintf( path, "%s2", argv[0] );
		if( argc > 1 ) {
			execl( path, path, argv[1], NULL );
		} else {
			execl( path, path, NULL );
		}
	}

	fprintf( stderr, "Error : fail at executing load module\n" );
	return -1;
}
