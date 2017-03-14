/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
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




