/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
	CHECDB_NGroup Ver.1.0
*/

#include "CHECDB.h"
#include "CHECData.h"

using namespace std;
using namespace hecd_util;


CHECDB_NGroup::CHECDB_NGroup()
 : CHECDataBlock(HECDB_NGROUP), NodeList()
{
	name[0] = 0;
}


CHECDB_NGroup::~CHECDB_NGroup()
{
}


void CHECDB_NGroup::Clear()
{
	NodeList.clear();
}


void CHECDB_NGroup::Write( CHECData* hecd )
{
	if( NodeList.size() == 0 ) return;

	hecd->WriteHeader( "!NGROUP", "S", "NGRP", name );

	set<int>::iterator iter;
	for(iter = NodeList.begin(); iter != NodeList.end(); iter++){
		hecd->WriteData( "I", *iter );
	}
}


bool CHECDB_NGroup::Read( CHECData* hecd, char* header_line )
{
	int rcode[5];
	if(! hecd->ParseHeader( header_line, rcode, "S", "NGRP", name )) return false;

	char line[256];
	const int max_id_n = 100;
	int id[max_id_n];
	int i,n;

	while(1) {
		if(!hecd->ReadLine(line)) break;
		if( line[0] == '!' ) {
			hecd->PushReadLine(line);
			break;
		}
		n = hecd->ParseIntDataArray( line, id );
		if(n<0) return false;
		for(i=0;i<n;i++) {
			NodeList.insert( id[i] );
		}
	}
	return true;

}
