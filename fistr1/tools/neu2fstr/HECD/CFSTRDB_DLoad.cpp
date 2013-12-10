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
	CFSTRDB_DLoad Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;
using namespace hecd_util;


// static method
int CFSTRDB_DLoad::ParamNumber( int type )
{
	const int pn[] = {
		1, // TYPE_P0
		1, // TYPE_P1
		1, // TYPE_P2
		1, // TYPE_P3
		1, // TYPE_P4
		1, // TYPE_P5
		1, // TYPE_P6
		1, // TYPE_BX
		1, // TYPE_BY
		1, // TYPE_BZ
		4, // TYPE_GRAV
		7, // TYPE_CENT
		0 // TYPE_UNKNOWN
	};

	if( type < 0 || type >= TypeNumber() ) return 0;
	return pn[type];
}


const char* CFSTRDB_DLoad::LoadTypeName( int type )
{
	const char* pn[] = {
		"P0", "P1", "P2", "P3", "P4", "P5", "P6",
		"BX", "BY", "BZ",
		"GRAV", "CENT",
		"unknown"
	};

	if( type < 0 || type >= TypeNumber() ) return "";
	return pn[type];
}




CFSTRDB_DLoad::CFSTRDB_DLoad()
 : CFSTRDataBlock(FSTRDB_DLOAD), ItemList()
{
}


CFSTRDB_DLoad::~CFSTRDB_DLoad()
{
	Clear();
}


void CFSTRDB_DLoad::Clear()
{
	ItemList.clear();
}


void CFSTRDB_DLoad::Write( CHECData* hecd )
{
	if( ItemList.size() == 0 ) return;

	hecd->WriteHeader( "!DLOAD" );

	vector<CItem>::iterator iter;
	for(iter = ItemList.begin(); iter != ItemList.end(); iter++){
		hecd->ClearDataLineBuffer();
		hecd->AddDataLineItems( "SS", iter->egrp, LoadTypeName( iter->type ));
		int n = ParamNumber( iter->type );
		for(int i=0; i<n; i++) {
			hecd->AddDataLineItems( "F", iter->param[i] );
		}
		hecd->WriteDataLine();
	}	
}


bool CFSTRDB_DLoad::Read( CHECData* hecd, char* header_line )
{
	int rcode[10];
	char s[256];
	int type;

	while(1) {
		CItem item;
		bool fg = hecd->ReadData( rcode, "SSFFFFFFF",
			item.egrp,
			s,
			&item.param[0],
			&item.param[1],
			&item.param[2],
			&item.param[3],
			&item.param[4],
			&item.param[5],
			&item.param[6]
		);
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

