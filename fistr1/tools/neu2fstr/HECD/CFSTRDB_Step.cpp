/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
	CFSTRDB_Step Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;
using namespace hecd_util;


CFSTRDB_Step::CFSTRDB_Step()
 : CFSTRDataBlock(FSTRDB_STEP),
   type(TYPE_STANDARD), incmax(100)
{
}


CFSTRDB_Step::~CFSTRDB_Step()
{
	Clear();
}


void CFSTRDB_Step::Clear()
{
	type = TYPE_STANDARD;
	incmax = 100;
}


void CFSTRDB_Step::Write( CHECData* hecd )
{
	if(type==TYPE_NLGEOM) {
		hecd->WriteHeader( "!STEP", "SI",
			"TYPE", "NLGEOM", "INCMAX", incmax );
	}
}


bool CFSTRDB_Step::Read( CHECData* hecd, char* header_line )
{
	char s[256] = "";
	char type_s[256];
	int rcode[10];

	if(!hecd->ParseHeader( header_line,rcode, "SI","TYPE", s, "INCMAX", &incmax )) return false;

	cleanup_token( s, type_s );
	toupper( type_s );

	if( strcmp(type_s, "STANDARD") == 0 ){
		type = TYPE_STANDARD;
	} else if( strcmp(type_s, "NLGEOM") == 0 ){
		type = TYPE_NLGEOM;
	} else
		return false;

	return true;
}







