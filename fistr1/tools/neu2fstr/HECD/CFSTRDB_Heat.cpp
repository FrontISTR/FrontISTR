/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
	CFSTRDB_Heat Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;
using namespace hecd_util;


CFSTRDB_Heat::CFSTRDB_Heat()
 : CFSTRDataBlock(FSTRDB_HEAT),
   restart(0), dt(0), etime(0), dtime(0), deltmx(0),
   itmax(20), eps(1e-6)
{
}


CFSTRDB_Heat::~CFSTRDB_Heat()
{
	Clear();
}


void CFSTRDB_Heat::Clear()
{
	restart = 0;
	dt = 0;
	etime = 0;
	dtime = 0;
	deltmx = 0;
	itmax = 20;
	eps = 1e-6;
}


void CFSTRDB_Heat::Write( CHECData* hecd )
{
	if( restart ){
		hecd->WriteHeader( "!HEAT", "S", "READ", "RESTART" );
	} else {
		hecd->WriteHeader( "!HEAT" );
	}
	hecd->WriteData( "FFFFIF", dt, etime, dtime, deltmx, itmax, eps );
}


bool CFSTRDB_Heat::Read( CHECData* hecd, char* header_line )
{
	int rcode[10];
	char s[256]="";
	char restart_s[256];

	if(! hecd->ParseHeader( header_line, rcode, "S", "READ", s )) return false;
	cleanup_token(s, restart_s);
	toupper( restart_s );
	restart = ( strcmp( restart_s, "RESTERT" )==0 );

	hecd->ReadData( rcode, "FFFFIF", &dt, &etime, &dtime, &deltmx, &itmax, &eps );

	return true;
}


