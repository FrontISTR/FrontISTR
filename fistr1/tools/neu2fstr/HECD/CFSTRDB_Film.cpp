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
	CFSTRDB_Film Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;
using namespace hecd_util;


// static method


const char* CFSTRDB_Film::LoadTypeName( int type )
{
	const char* pn[] = {
		"F0", "F1", "F2", "F3", "F4", "F5", "F6",
		"unknown"
	};

	if( type < 0 || type >= TypeNumber() ) return "";
	return pn[type];
}


CFSTRDB_Film::CFSTRDB_Film()
 : CFSTRDataBlock(FSTRDB_FILM), ItemList()
{
	amp1[0] = 0;
	amp2[0] = 0;
}


CFSTRDB_Film::~CFSTRDB_Film()
{
	Clear();
}


void CFSTRDB_Film::Clear()
{
	ItemList.clear();
	amp1[0] = 0;
	amp2[0] = 0;
}


void CFSTRDB_Film::Write( CHECData* hecd )
{
	char buff[256];

	if( ItemList.size() == 0 ) return;

	strcpy( buff, "!FILM");
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
		hecd->WriteData( "SSFF",
			iter->egrp, LoadTypeName( iter->type ), iter->value, iter->sink );
	}	
}



bool CFSTRDB_Film::Read( CHECData* hecd, char* header_line )
{
	int rcode[10];
	char s[256];
	int type;

	amp1[0] = 0;
	amp2[0] = 0;
	if(!hecd->ParseHeader( header_line, rcode, "SS", "AMP1", amp1, "AMP2", amp2 )) return false;

	while(1) {
		CItem item;
		bool fg = hecd->ReadData( rcode, "SSFF", item.egrp, s, &item.value, &item.sink );
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

