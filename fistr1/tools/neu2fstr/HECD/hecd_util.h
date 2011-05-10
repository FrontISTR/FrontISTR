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

#ifndef hecd_utilH
#define hecd_utilH


#include <stdio.h>
#include <string.h>
#include <ctype.h>


namespace hecd_util {

void cleanup_token( char* s );
void cleanup_token( char* src, char* dest );
void toupper( char* s );
void toupper( const char* src, char* dest );
void tolower( char* s );
void tolower( const char* src, char* dest );
void remove_cr( char* s );
void ftos( double x, char* s );

} // end of namespace hecd_util


#endif
