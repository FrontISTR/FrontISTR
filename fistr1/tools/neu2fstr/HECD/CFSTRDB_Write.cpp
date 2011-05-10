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
	CFSTRDB_Write Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;


CFSTRDB_Write::CFSTRDB_Write()
 : CFSTRDataBlock(FSTRDB_WRITE),
   result(0), visual(0)
{
}

CFSTRDB_Write::~CFSTRDB_Write()
{
	Clear();
}

void CFSTRDB_Write::Clear()
{
	result = visual = 0;
}

void CFSTRDB_Write::Write( CHECData* hecd )
{
	char header_s[256];

	strcpy( header_s, "!WRITE" );
	if( result ) strcat( header_s, ",RESULT");
	if( visual ) strcat( header_s, ",VISUAL");

	hecd->WriteHeader( header_s ) ;
}


bool CFSTRDB_Write::Read( CHECData* hecd, char* header_line )
{
	int rcode[10];

	if(! hecd->ParseHeader( header_line, rcode, "EE", "RESULT",&result,"VISUAL",&visual)) return false;
	return true;
}

