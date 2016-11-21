/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
	CFSTRDB_Static Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;


CFSTRDB_Static::CFSTRDB_Static()
 : CFSTRDataBlock(FSTRDB_STATIC),
   dtime(0), etime(0), itmax(20), eps(1e-6)
{
}

CFSTRDB_Static::~CFSTRDB_Static()
{
	Clear();
}

void CFSTRDB_Static::Clear()
{
	dtime = 0;
	etime = 0;
	itmax = 20;
	eps = 1e-6;
}

void CFSTRDB_Static::Write( CHECData* hecd )
{
	hecd->WriteHeader( "!STATIC" );
	hecd->WriteData( "FFIF", dtime, etime, itmax, eps );
}


bool CFSTRDB_Static::Read( CHECData* hecd, char* header_line )
{
	int rcode[10];
	hecd->ReadData( rcode, "FFIF", &dtime, &etime, &itmax, &eps );
	return true;
}

