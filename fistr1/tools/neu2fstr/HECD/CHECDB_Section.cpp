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
	CHECDB_Section Ver.1.0
*/

#include "CHECDB.h"
#include "CHECData.h"

using namespace std;
using namespace hecd_util;


CHECDB_Section::CHECDB_Section()
 : CHECDataBlock(HECDB_SECTION),
   type( TYPE_UNKNOWN ),
   n_comp(0),
   secopt(0),
   thickness(1.0),
   integpoints(0),
   gapcon(0),
   gaprad1(0),
   gaprad2(0)
{
	egrp[0] = 0;
	material[0] = 0;
}


CHECDB_Section::~CHECDB_Section()
{
}


void CHECDB_Section::Clear()
{
	type = TYPE_UNKNOWN;
	n_comp = 0;
	secopt = 0;
	thickness = 1.0;
	integpoints = 0;
	gapcon = 0;
	gaprad1 = 0;
	gaprad2 = 0;
	egrp[0] = 0;
	material[0] = 0;
}


void CHECDB_Section::Write( CHECData* hecd )
{
	switch( type ){
	case TYPE_SOLID:
		hecd->WriteHeader( "!SECTION", "SSS", "TYPE", "SOLID", "EGRP", egrp, "MATERIAL", material );
		hecd->WriteData("F", thickness );
		break;
	case TYPE_SHELL:
		hecd->WriteHeader( "!SECTION", "SSS", "TYPE", "SHELL", "EGRP", egrp, "MATERIAL", material );
		hecd->WriteData("FI", thickness, integpoints );
		break;
	case TYPE_INTERFACE:
		hecd->WriteHeader( "!SECTION", "SSS", "INTERFACE", "SOLID", "EGRP", egrp, "MATERIAL", material );
		hecd->WriteData("FFFF", thickness, gapcon, gaprad1, gaprad2 );
	default:
		assert(0);
	}
}


bool CHECDB_Section::Read( class CHECData* hecd, char* header_line )
{
	int rcode[10];
	char s[256], type_s[256];

	if(!hecd->ParseHeader( header_line, rcode, "SSS",
		"TYPE", s, "EGRP", egrp, "MATERIAL", material)) return false;

	cleanup_token( s, type_s );
	toupper( type_s );
	if( strcmp( type_s, "SOLID" )==0) {
		type = TYPE_SOLID;
		if(! hecd->ReadData( rcode, "F", &thickness )) return false;
	} else if( strcmp( type_s, "SHELL" )==0) {
		type = TYPE_SHELL;
		if(! hecd->ReadData( rcode, "FI", &thickness, &integpoints )) return false;
	} else if( strcmp( type_s, "INTERFACE" )==0) {
		type = TYPE_INTERFACE;
		if(! hecd->ReadData( rcode, "FFFF", &thickness, &gapcon, &gaprad1, &gaprad2 )) return false;
	} else {
		return false;
	}

	return true;
}


