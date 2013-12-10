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
	CNFDB_405 Ver.1.0
*/


// 405 Coordnate Systems

#include "CNFData.h"
#include "CNFDB_405.h"



CNFDB_405::CNFDB_405()
 : CNFDataBlock(405)
{}


void CNFDB_405::Read( CNFData* nfd )
{
	char buff[256];

	// #1
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIIII", 
			&ID,
			&define_sys,
			&type,
			&color,
			&layer );
	// #2
		nfd->ReadLineEx( buff );
		nfd->ReadStr(buff, title, sizeof(title));
	// #3
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFF", &origin[0], &origin[1], &origin[2] );
	// #4
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFF", &rot[0], &rot[1], &rot[2] );
}


void CNFDB_405::WriteData( class CNFData* nfd, FILE* fp )
{
	// #1
		nfd->WriteData( fp, "IIIIIn", 
			ID,
			define_sys,
			type,
			color,
			layer );
	// #2
		nfd->WriteStr( fp, title);
	// #3
		nfd->WriteData( fp, "FFFn", origin[0], origin[1], origin[2] );
	// #4
		nfd->WriteData( fp, "FFFn", rot[0], rot[1], rot[2] );
}

