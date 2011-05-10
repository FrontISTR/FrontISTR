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
	CFSTRDB_Echo Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;


CFSTRDB_Echo::CFSTRDB_Echo()
 : CFSTRDataBlock(FSTRDB_ECHO),
   echo(1)
{
}

CFSTRDB_Echo::~CFSTRDB_Echo()
{
	Clear();
}

void CFSTRDB_Echo::Clear()
{
	echo = 1;
}

void CFSTRDB_Echo::Write( CHECData* hecd )
{
	if( echo ) hecd->WriteHeader( "!ECHO" );
}


bool CFSTRDB_Echo::Read( CHECData* hecd, char* header_line )
{
	echo = 1;
	return true;
}

