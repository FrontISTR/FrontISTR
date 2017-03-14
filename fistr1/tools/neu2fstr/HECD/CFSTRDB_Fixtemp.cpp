/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
	CFSTRDB_Fixtemp Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;


CFSTRDB_Fixtemp::CFSTRDB_Fixtemp()
 : CFSTRDataBlock(FSTRDB_FIXTEMP), ItemList()
{
	amp[0] = 0;
}


CFSTRDB_Fixtemp::~CFSTRDB_Fixtemp()
{
	Clear();
}


void CFSTRDB_Fixtemp::Clear()
{
	ItemList.clear();
	amp[0] = 0;
}


void CFSTRDB_Fixtemp::Write( CHECData* hecd )
{
	if( ItemList.size() == 0 ) return;

	if( amp[0] == 0 ) {
		hecd->WriteHeader( "!FIXTEMP" );
	} else {
		hecd->WriteHeader( "!FIXTEMP", "S", "AMP", amp );
	}

	vector<CItem>::iterator iter;
	for(iter = ItemList.begin(); iter != ItemList.end(); iter++){
		hecd->WriteData( "SF", iter->ngrp, iter->value );
	}
}


bool CFSTRDB_Fixtemp::Read( CHECData* hecd, char* header_line )
{
	int rcode[10];
	amp[0] = 0;
	if(!hecd->ParseHeader( header_line, rcode, "S", "AMP", amp )) return false;

	while(1) {
		CItem item;
		if( !hecd->ReadData( rcode, "SF", item.ngrp, &item.value )) break;
		ItemList.push_back(item);
	}
	return true;
}
