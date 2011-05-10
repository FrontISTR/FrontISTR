/*=====================================================================*
 *                                                                     *
 *   Software Name : neu2fstr                                          *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : FrontSTR Input File Converter                     *
 *                                                                     *
 *            Written by Noboru Imai (Univ. of Tokyo)                  *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using Hight End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/

/*
	hecd_util ver.1.0
*/

#include "hecd_util.h"


namespace hecd_util {

void cleanup_token( char* s )
{
	char buff[256];
	cleanup_token( s, buff );
	strcpy( s, buff );
}

void cleanup_token( char* src, char* dest )
{
	#define is_skip_char(x) ( x==' ' || x=='=' || x=='\t' || x=='\r' || x=='\n' )

	char* s = src;
	while( *s && is_skip_char(*s)) s++;
	char* d = dest;
	while( *s && !is_skip_char(*s)) {
		*d = *s;
		d++;
		s++;
	}
	*d = 0;
}


void toupper( char* s )
{
	while(*s) {
		*s = (char)::toupper(*s);
		s++;
	}
}


void toupper( const char* src, char* dest )
{
	char* s = (char*)src;
	while(*s) {
		*dest = (char)::toupper(*s);
		s++;
		dest++;
	}
	*dest = 0;
}



void tolower( char* s )
{
	while(*s) {
		*s = (char)::tolower(*s);
		s++;
	}
}


void tolower( const char* src, char* dest )
{
	char* s = (char*)src;
	while(*s) {
		*dest = (char)::tolower(*s);
		s++;
		dest++;
	}
	*dest = 0;
}


void remove_cr( char* s )
{
	while(*s) {
		if( *s == '\r' || *s == '\n' ) {
			*s = 0;
			return;
		}
		s++;
	}
}

// note)
//  I/O of HEC-MW does not support '1e+2' formated number.
//  Then ftos converts '1e+2' to '1.0e+2'

void ftos( double x, char* s )
{
	char buff[256];
	sprintf( buff, "%.10lg", x );

	char* p = buff;
	bool fg_dot = false;
	while(*p) {
		if( *p == '.' ) {
			fg_dot = true;
		} else if( *p == 'e' || *p == 'E' ){
			if( !fg_dot ){
				*s = '.'; s++;
				*s = '0'; s++;
			}
		}
		*s = *p;
		p++;
		s++;
	}
	*s = 0;
}


} // end of namespace hecd_util

