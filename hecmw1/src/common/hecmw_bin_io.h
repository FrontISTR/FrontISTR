/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2007/12/03                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Noboru Imai (Univ. of Tokyo)                  *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/

#ifndef hecmw_bin_ioH
#define hecmw_bin_ioH

#include <stdio.h>

/*---------------------------------------------------------------------------*/
/* CAUTION) hecmw_set_endian_info must be executed before calling following functions. */

void hecmw_set_endian_info(void);

/*---------------------------------------------------------------------------*/

int hecmw_write_bin_value( unsigned char* x, int size, FILE* fp );
int hecmw_write_bin_int( int x, FILE* fp );
int hecmw_write_bin_int_arr( int* x, int n, FILE* fp );
int hecmw_write_bin_double( double x, FILE* fp );
int hecmw_write_bin_double_arr( double* x, int n, FILE* fp );

/*---------------------------------------------------------------------------*/

int hecmw_read_bin_value( unsigned char* x, int size, FILE* fp );
int hecmw_read_bin_int( int* x, FILE* fp );
int hecmw_read_bin_int_arr( int* x, int n, FILE* fp );
int hecmw_read_bin_double( double* x, FILE* fp );
int hecmw_read_bin_double_arr( double* x, int n, FILE* fp );

/*---------------------------------------------------------------------------*/

/*
 * fmt : type and array size of arguments
 *       Each cherecter in fmt specifies type of argument and
 *       number after the cheracter means array size.
 * meening of character in fmt
 * 'I' : int
 * 'F' : double
 * 'S' : string (char*)
 * exp)
 *     int n, i[10]; char s[20] = "123";
 *     hecmw_write_bin(fp, "II10S", n, i, s);
 */

int hecmw_write_bin( FILE* fp, const char* fmt, ... );
int hecmw_read_bin( FILE* fp, const char* fmt, ... );


#endif
