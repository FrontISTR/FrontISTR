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
	CHECDB_Header Ver.1.0
*/

#include "CHECDB.h"
#include "CHECData.h"


CHECDB_Header::CHECDB_Header()
 : CHECDataBlock( HECDB_HEADER )
{
	title[0] = 0;
}

CHECDB_Header::~CHECDB_Header()
{
	Clear();
}


void CHECDB_Header::Clear()
{
	title[0] = 0;
}

void CHECDB_Header::Write( CHECData* hecd )
{
	if( title[0] == 0 ) return;
	hecd->WriteHeader("!HEADER");
	hecd->WriteLine( title );
}


bool CHECDB_Header::Read( CHECData* hecd, char* header_line )
{
	char line[256];
	if(!hecd->ReadLine(line)) return false;
	if( line[0] == '!') {
		hecd->PushReadLine(line);
		return true;
	}
	strcpy( title, line );
	return true;
}

