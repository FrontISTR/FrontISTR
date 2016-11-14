/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
	CNFDB_403 Ver.1.0
*/


// 403 Node

#include "CNFData.h"
#include "CNFDB_403.h"


CNFDB_403::CNFDB_403()
 : CNFDataBlock(403)
{
	node_type = 0;
}


void CNFDB_403::Read( CNFData* nfd )
{
	char buff[256];

	// #1
		nfd->ReadLineEx( buff );
		if( nfd->version >= 4.4 ) {
			nfd->ReadRecord( buff, "IIIIIBBBBBBFFFI",
				&ID,
				&define_sys,
				&output_sys,
				&layer,
				&color,
				&permbc[0],
				&permbc[1],
				&permbc[2],
				&permbc[3],
				&permbc[4],
				&permbc[5],
				&x,
				&y,
				&z,
				&node_type
			);
		} else {
			nfd->ReadRecord( buff, "IIIIIBBBBBBFFF",
				&ID,
				&define_sys,
				&output_sys,
				&layer,
				&color,
				&permbc[0],
				&permbc[1],
				&permbc[2],
				&permbc[3],
				&permbc[4],
				&permbc[5],
				&x,
				&y,
				&z
			);
		}
}


void CNFDB_403::WriteData( class CNFData* nfd, FILE* fp )
{
	// #1
		if( nfd->version >= 4.4 ) {
			nfd->WriteData( fp, "IIIIIBBBBBBFFFIn",
				ID,
				define_sys,
				output_sys,
				layer,
				color,
				permbc[0],
				permbc[1],
				permbc[2],
				permbc[3],
				permbc[4],
				permbc[5],
				x,
				y,
				z,
				node_type
			);
		} else {
			nfd->WriteData( fp, "IIIIIBBBBBBFFFn",
				ID,
				define_sys,
				output_sys,
				layer,
				color,
				permbc[0],
				permbc[1],
				permbc[2],
				permbc[3],
				permbc[4],
				permbc[5],
				x,
				y,
				z
			);
		}
}


