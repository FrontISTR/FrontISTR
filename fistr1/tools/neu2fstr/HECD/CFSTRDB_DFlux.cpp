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
	CFSTRDB_DFlux Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;
using namespace hecd_util;


// static method


const char* CFSTRDB_DFlux::LoadTypeName( int type )
{
	const char* pn[] = {
		"S0", "S1", "S2", "S3", "S4", "S5", "S6",
		"BF", 
		"unknown",
	};

	if( type < 0 || type >= TypeNumber() ) return "";
	return pn[type];
}




CFSTRDB_DFlux::CFSTRDB_DFlux()
 : CFSTRDataBlock(FSTRDB_DFLUX), ItemList()
{
	amp[0] = 0;
}


CFSTRDB_DFlux::~CFSTRDB_DFlux()
{
	Clear();
}


void CFSTRDB_DFlux::Clear()
{
	ItemList.clear();
	amp[0] = 0;
}


void CFSTRDB_DFlux::Write( CHECData* hecd )
{
	if( ItemList.size() == 0 ) return;
	if( amp[0] == 0 ) {
		hecd->WriteHeader( "!DFLUX" );	
	} else {
		hecd->WriteHeader( "!DFLUX", "S", "AMP", amp );
	}

	vector<CItem>::iterator iter;
	for(iter = ItemList.begin(); iter != ItemList.end(); iter++){
		hecd->WriteData( "SSF", iter->egrp, LoadTypeName( iter->type ), iter->value);
	}	
}


bool CFSTRDB_DFlux::Read( CHECData* hecd, char* header_line )
{
	int rcode[10];
	char s[256];
	int type;

	amp[0] = 0;
	if(!hecd->ParseHeader( header_line, rcode, "S", "AMP", amp )) return false;

	while(1) {
		CItem item;
		bool fg = hecd->ReadData( rcode, "SSF", item.egrp, s, &item.value );
		if( !fg ) break;

		cleanup_token(s);
		toupper(s);
		for( type=0; type<TypeNumber(); type++) {
			if( strcmp( LoadTypeName(type), s )==0) break;
		}
		if( type == TypeNumber()) return false;
		item.type = type; 
		ItemList.push_back( item );
	}
	return true;
}

