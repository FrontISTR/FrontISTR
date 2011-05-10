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
	CFSTRDB_CLoad Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;


CFSTRDB_CLoad::CFSTRDB_CLoad()
 : CFSTRDataBlock(FSTRDB_CLOAD), ItemList()
{
}


CFSTRDB_CLoad::~CFSTRDB_CLoad()
{
	Clear();
}


void CFSTRDB_CLoad::Clear()
{
	ItemList.clear();
}


void CFSTRDB_CLoad::Write( CHECData* hecd )
{
	if( ItemList.size() == 0 ) return;
	hecd->WriteHeader( "!CLOAD" );

	vector<CItem>::iterator iter;
	for(iter = ItemList.begin(); iter != ItemList.end(); iter++){
		hecd->WriteData( "SIF", iter->ngrp, iter->dof_id, iter->value );
	}	
}


bool CFSTRDB_CLoad::Read( CHECData* hecd, char* header_line )
{
	int rcode[10];

	while(1) {
		CItem item;
		if(!hecd->ReadData( rcode, "SIF", &item.ngrp, &item.dof_id, &item.value ))
			break;
		ItemList.push_back( item );
	}
	return true;
}

