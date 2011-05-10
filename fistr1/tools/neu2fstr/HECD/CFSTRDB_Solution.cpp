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
	CFSTRDB_Solution Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;
using namespace hecd_util;


CFSTRDB_Solution::CFSTRDB_Solution()
 : CFSTRDataBlock(FSTRDB_SOLUTION), type(TYPE_STATIC)
{
}

CFSTRDB_Solution::~CFSTRDB_Solution()
{
	Clear();
}

void CFSTRDB_Solution::Clear()
{
	type = TYPE_STATIC;
}

void CFSTRDB_Solution::Write( CHECData* hecd )
{
	if( type == TYPE_UNKNOWN ) return;

	switch(type){
	case TYPE_STATIC:
		hecd->WriteHeader( "!SOLUTION", "S", "TYPE", "STATIC" );
		break;
	case TYPE_HEAT:
		hecd->WriteHeader( "!SOLUTION", "S", "TYPE", "HEAT" );
		break;
	case TYPE_EIGEN:
		hecd->WriteHeader( "!SOLUTION", "S", "TYPE", "EIGEN" );
		break;
	default:
		assert(0);
	}
}


bool CFSTRDB_Solution::Read( CHECData* hecd, char* header_line )
{
	int rcode[5];
	char s[256];
	char type_s[256];
	
	if(! hecd->ParseHeader( header_line, rcode, "S", "TYPE", s)) return false;

	cleanup_token( s, type_s );
	toupper( type_s );

	if( strcmp(type_s, "STATIC") == 0 ){
		type = TYPE_STATIC;
	} else if( strcmp(type_s, "HEAT") == 0 ){
		type = TYPE_HEAT;
	} else if( strcmp(type_s, "EIGEN") == 0 ){
		type = TYPE_EIGEN;
	} else
		return false;

	return true;
}

