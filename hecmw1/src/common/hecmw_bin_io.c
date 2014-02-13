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

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>

#include "hecmw_bin_io.h"


/*
	format of fmt
	'I'    : int
	'F'    : double
	'I<n>' : int[n]
	'F<n>' : double[n]
	'S'    : char[]
*/


#ifndef SUCCESS
#define SUCCESS 0
#define FAIL -1
#endif

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

static
int fg_little_endian = TRUE;

#define RES_BIN_HEADER "HECMW_BINARY_RESULT"

/*---------------------------------------------------------------------------*/

union endian_check_u {
   	char  c[2];
	short s;
};

static
void check_endian(void)
{
	union endian_check_u u;
	u.c[0] = 0;
	u.c[1] = 0;
	u.s = 1;
	fg_little_endian = (u.c[0] == 1);
}

void hecmw_set_endian_info(void)
{
	check_endian();
}

/*---------------------------------------------------------------------------*/

int hecmw_write_bin_value( unsigned char* x, int size, FILE* fp )
{
	unsigned char *c;
	int i;

	if( fg_little_endian ){
		c = x;
		for(i=0; i<size; i++ ){
			if( putc( *c, fp )==EOF ) return FAIL;
			c++;
		}
	} else {
		c = x;
		c += size-1;
		for(i=0; i<size; i++ ){
			if( putc( *c, fp )==EOF ) return FAIL;
			c--;
		}
	}
	return SUCCESS;
}


int hecmw_write_bin_int( int x, FILE* fp )
{
	long xx = x;
	return hecmw_write_bin_value( (unsigned char*)&xx, sizeof(long), fp);
}


int hecmw_write_bin_int_arr( int* x, int n, FILE* fp )
{
	long xx;
	int i;
	for(i=0; i<n; i++, x++) {
		xx = *x;
		if(hecmw_write_bin_value( (unsigned char*)&xx, sizeof(long), fp)!=SUCCESS) return FAIL;
	}
	return SUCCESS;
}


int hecmw_write_bin_double( double x, FILE* fp )
{
	return hecmw_write_bin_value((unsigned char*)&x, sizeof(double), fp);
}



int hecmw_write_bin_double_arr( double* x, int n, FILE* fp )
{
	int i;

	for(i=0; i<n; i++, x++) {
		if(hecmw_write_bin_value((unsigned char*)&x, sizeof(double), fp)!=SUCCESS) return FAIL;
	}
	return SUCCESS;
}


/*---------------------------------------------------------------------------*/


int hecmw_read_bin_value( unsigned char* x, int size, FILE* fp )
{
	unsigned char *c;
	int data;
	int i;

	if( fg_little_endian ){
		c = x;
		for(i=0; i<size; i++ ){
			data = getc(fp);
			if( data == EOF ) return FAIL;
			*c = (unsigned char)data;
			c++;
		}
	} else {
		c = x;
		c += size-1;
		for(i=0; i<size; i++ ){
			data = getc(fp);
			if( data == EOF ) return FAIL;
			*c = (unsigned char)data;
			c--;
		}
	}
	return SUCCESS;
}


int hecmw_read_bin_int( int* x, FILE* fp )
{
	long xx;
	if(hecmw_read_bin_value( (unsigned char*)&xx, sizeof(long), fp)!=SUCCESS) {
		return FAIL;
	}
	*x = xx;
	return SUCCESS;
}



int hecmw_read_bin_int_arr( int* x, int n, FILE* fp )
{
	long xx;
	int i;
	for( i=0; i<n; i++, x++) {
		if(hecmw_read_bin_value( (unsigned char*)&xx, sizeof(long), fp)!=SUCCESS) {
			return FAIL;
		}
		*x = xx;
	}
	return SUCCESS;
}


int hecmw_read_bin_double( double* x, FILE* fp )
{
	if(hecmw_read_bin_value( (unsigned char*)x, sizeof(double), fp)!=SUCCESS) {
		return FAIL;
	}
	return SUCCESS;
}


int hecmw_read_bin_double_arr( double* x, int n, FILE* fp )
{
	int i;
	for( i=0; i<n; i++, x++) {
		if(hecmw_read_bin_value( (unsigned char*)x, sizeof(double), fp)!=SUCCESS) return FAIL;
	}
	return SUCCESS;
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/


static
int get_fmt_type_size( char** fmt_p, char* type, int* size )
{
	int i;
	char s[256];
	char* sp = s;
	char* p = *fmt_p;

	if( !*p ) return 0; /* FAIL */

	*type = *p;
	p++;

	i = 0;
	while( *p && '0' <= *p && *p <= '9' ) {
		*sp = *p;
		p++;
		sp++;
		i++;
	}
	*sp = 0;

	 p = *fmt_p;
	if( i>0 ){
		sscanf( s, "%d", size );
		p += i+1;
	} else {
		*size = 0;
		p++;
	}
	*fmt_p = p;
	return 1; /* SUCCESS */
}



int hecmw_write_bin( FILE* fp, const char* fmt, ... )
{
	int i, n;
	char type;
	int size;
	char* fmt_p;
	int i_data;
	double d_data;
	int* i_ptr;
	double* d_ptr;
	char* s_ptr;

	va_list va;
	va_start(va, fmt);
	n = strlen((char*)fmt);
	i = 0;
	fmt_p = (char*)fmt;

	while( get_fmt_type_size( &fmt_p, &type, &size )) {
		switch(type){
		case 'I':case 'i': /* int */
			if( size > 0 ){
				i_ptr = va_arg(va, int*);
				if(hecmw_write_bin_int_arr( i_ptr, size, fp )) goto ERROR_EXIT;
			} else {
				i_data = va_arg(va, int);
				if(hecmw_write_bin_int( i_data, fp )) goto ERROR_EXIT;
			}
			break;
		case 'F':case 'f': /* double */
			if( size > 0 ){
				d_ptr = va_arg(va, double*);
				if(hecmw_write_bin_double_arr( d_ptr, size, fp )) goto ERROR_EXIT;
			} else {
				d_data = va_arg(va, double);
				if(hecmw_write_bin_double( d_data, fp )) goto ERROR_EXIT;
			}
			break;
		case 'S':case 's': /* string */
			if( size > 0 ){
				s_ptr = va_arg(va, char*);
				if(fwrite( s_ptr, sizeof(char), size, fp )!=size) goto ERROR_EXIT;
			} else {
				s_ptr = va_arg(va, char*);
				while( *s_ptr ) {
					if(fwrite( s_ptr, sizeof(char),1, fp ) !=1 ) goto ERROR_EXIT;
					s_ptr++;
				}
				if(fwrite( s_ptr, sizeof(char),1, fp ) !=1 ) goto ERROR_EXIT;
			}
			break;
		default:
			fprintf(stderr, "Illeagal type : %c (0x%x)\n", type, type );
			assert(0);
		}
	}
	va_end(va);
	return SUCCESS;

	ERROR_EXIT:
	va_end(va);
	return FAIL;
}



int hecmw_read_bin( FILE* fp, const char* fmt, ... )
{
	int i, n;
	char type;
	int size;
	char* fmt_p;
	int* i_ptr;
	double* d_ptr;
	char* s_ptr;
	char c;

	va_list va;
	va_start(va, fmt);
	n = strlen((char*)fmt);
	i = 0;
	fmt_p = (char*)fmt;

	while( get_fmt_type_size( &fmt_p, &type, &size )) {
		switch(type){
		case 'I':case 'i': /* int */
			if( size == 0 ) size = 1;
			i_ptr = va_arg(va, int*);
			if(hecmw_read_bin_int_arr( i_ptr, size, fp )) goto ERROR_EXIT;
			break;
		case 'F':case 'f': /* double */
			if( size == 0 ) size = 1;
			d_ptr = va_arg(va, double*);
			if(hecmw_read_bin_double_arr( d_ptr, size, fp )) goto ERROR_EXIT;
			break;
		case 'S':case 's': /* string */
			if( size > 0 ) {
				s_ptr = va_arg(va, char*);
				if(fread( s_ptr, sizeof(char), size, fp ) != size) goto ERROR_EXIT;
			} else {
				s_ptr = va_arg(va, char*);
				for(;;){
					if( fread( &c, sizeof(char),1, fp ) != 1) break;
					*s_ptr = c;
					if( c == 0 ) break;
					s_ptr++;
				}
			}
			break;
		default:
			assert(0);
		}
	}
	va_end(va);
	return SUCCESS;
ERROR_EXIT:
	va_end(va);
	return FAIL;
}



