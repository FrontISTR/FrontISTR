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
	CFSTRDB_SFilm Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;
using namespace hecd_util;



CFSTRDB_SFilm::CFSTRDB_SFilm()
 : CFSTRDataBlock(FSTRDB_SFILM), ItemList()
{
	amp1[0] = 0;
	amp2[0] = 0;
}


CFSTRDB_SFilm::~CFSTRDB_SFilm()
{
	Clear();
}


void CFSTRDB_SFilm::Clear()
{
	ItemList.clear();
	amp1[0] = 0;
	amp2[0] = 0;
}


void CFSTRDB_SFilm::Write( CHECData* hecd )
{
	char buff[256];

	if( ItemList.size() == 0 ) return;

	strcpy( buff, "!SFILM");
	if( amp1[0] != 0 ) {
		strcat( buff, ",AMP1=");
		strcat( buff, amp1);
	}
	if( amp2[0] != 0 ) {
		strcat( buff, ",AMP2=");
		strcat( buff, amp2);
	}
	hecd->WriteHeader( buff );

	vector<CItem>::iterator iter;
	for(iter = ItemList.begin(); iter != ItemList.end(); iter++){
		hecd->WriteData( "SFF", iter->sgrp, iter->value, iter->sink);
	}	
}


bool CFSTRDB_SFilm::Read( CHECData* hecd, char* header_line )
{
	int rcode[10];

	amp1[0] = 0;
	amp2[0] = 0;
	if(!hecd->ParseHeader( header_line, rcode, "SS", "AMP1", amp1, "AMP2", amp2 )) return false;

	while(1) {
		CItem item;
		bool fg = hecd->ReadData( rcode, "SFF", item.sgrp, &item.value, &item.sink );
		if( !fg ) break;
		ItemList.push_back( item );
	}
	return true;
}

