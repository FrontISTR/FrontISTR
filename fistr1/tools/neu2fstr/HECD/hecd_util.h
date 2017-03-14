/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
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
