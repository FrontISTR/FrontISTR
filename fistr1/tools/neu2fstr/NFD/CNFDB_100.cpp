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
	CNFDB_100 Ver.1.0
	-----------------------------
	100 Header
*/

#include "CNFData.h"
#include "CNFDB_100.h"

// 100 Header

CNFDB_100::CNFDB_100()
 : CNFDataBlock(100)
{
	title[0] = 0;
	version = -1;
}


void CNFDB_100::Read( CNFData* nfd )
{
	char buff[256];

	// #1
		nfd->ReadLineEx( buff );
		nfd->ReadStr( buff, title, sizeof(title));

	// #2
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "F", &version );
}


void CNFDB_100::WriteData( CNFData* nfd, FILE* fp )
{
	// #1
		nfd->WriteStr( fp, title );

	// #2
		fprintf(fp, "%.1lf,\n", version );
}




