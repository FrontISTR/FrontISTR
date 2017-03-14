/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
	CFSTRDB_Echo Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;


CFSTRDB_Echo::CFSTRDB_Echo()
 : CFSTRDataBlock(FSTRDB_ECHO),
   echo(1)
{
}

CFSTRDB_Echo::~CFSTRDB_Echo()
{
	Clear();
}

void CFSTRDB_Echo::Clear()
{
	echo = 1;
}

void CFSTRDB_Echo::Write( CHECData* hecd )
{
	if( echo ) hecd->WriteHeader( "!ECHO" );
}


bool CFSTRDB_Echo::Read( CHECData* hecd, char* header_line )
{
	echo = 1;
	return true;
}

