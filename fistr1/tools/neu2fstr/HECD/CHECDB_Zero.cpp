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
	CHECDB_Zero Ver.1.0
*/

#include "CHECDB.h"
#include "CHECData.h"


CHECDB_Zero::CHECDB_Zero()
 : CHECDataBlock(HECDB_ZERO), zero(0)
{
}



CHECDB_Zero::~CHECDB_Zero()
{
	Clear();
}


void CHECDB_Zero::Clear()
{
	zero = 0;
}


void CHECDB_Zero::Write( CHECData* hecd )
{
	hecd->WriteHeader("!ZERO");
	hecd->WriteData( "F", zero );
}


bool CHECDB_Zero::Read( CHECData* hecd, char* header_line )
{
	int rcode[10];
	return hecd->ReadData( rcode, "F", &zero );
}



