/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
	CFSTRDB_Temperature Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;


CFSTRDB_Temperature::CFSTRDB_Temperature()
 : CFSTRDataBlock(FSTRDB_TEMPERATURE), ItemList()
{
}


CFSTRDB_Temperature::~CFSTRDB_Temperature()
{
	Clear();
}


void CFSTRDB_Temperature::Clear()
{
	ItemList.clear();
}


void CFSTRDB_Temperature::Write( CHECData* hecd )
{
	if( ItemList.size() == 0 ) return;

	hecd->WriteHeader( "!TEMPERATURE" );

	vector<CItem>::iterator iter;
	for(iter = ItemList.begin(); iter != ItemList.end(); iter++){
		hecd->WriteData( "SF", iter->ngrp, iter->value );
	}
}


bool CFSTRDB_Temperature::Read( CHECData* hecd, char* header_line )
{
	int rcode[5];
	while(1) {
		CItem item;
		if(!hecd->ReadData( rcode, "SF", item.ngrp, &item.value )) break;
		ItemList.push_back( item );
	}
	return true;
}

