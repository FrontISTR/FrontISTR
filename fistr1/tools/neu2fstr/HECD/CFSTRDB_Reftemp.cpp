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
	CFSTRDB_Reftemp Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;


CFSTRDB_Reftemp::CFSTRDB_Reftemp()
 : CFSTRDataBlock(FSTRDB_REFTEMP),
   value(0)
{
}

CFSTRDB_Reftemp::~CFSTRDB_Reftemp()
{
	Clear();
}

void CFSTRDB_Reftemp::Clear()
{
	value = 0;
}

void CFSTRDB_Reftemp::Write( CHECData* hecd )
{
	hecd->WriteHeader( "!REFTEMP" );
	hecd->WriteData( "F", value );
}


bool CFSTRDB_Reftemp::Read( CHECData* hecd, char* header_line )
{
	int rcode[5];
	return hecd->ReadData( rcode, "F", &value );
}
