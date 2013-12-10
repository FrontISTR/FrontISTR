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
	CNFDB_601 Ver.1.0
	-----------------------------
	601 Material
*/


#include <vector>
#include "CNFData.h"
#include "CNFDB_601.h"

using namespace std;

// 601 Material

CNFDB_601::CNFDB_601()
 : CNFDataBlock(601)
{}



void CNFDB_601::Read( CNFData* nfd )
{
	int i;
	char buff[256];

	// #1
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIIIIII", 
			&ID,
			&format,
			&color,
			&type,
			&subtype,
			&layer,
			&FunctionCount );
	// #2
		nfd->ReadLineEx( buff );
		nfd->ReadStr( buff, title, sizeof(title));
	// #3
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "I", &Bcount);
	// #4
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIIIIIIIII", 
			&bval[0], &bval[1], &bval[2], &bval[3], &bval[4],
			&bval[5], &bval[6], &bval[7], &bval[8], &bval[9] );
	// #5
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "I", &Icount);
	// #6
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIIIIIIIII", 
			&ival[0], &ival[1], &ival[2], &ival[3], &ival[4],
			&ival[5], &ival[6], &ival[7], &ival[8], &ival[9] );
	// #7
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIIIIIIIII", 
			&ival[10], &ival[11], &ival[12], &ival[13], &ival[14],
			&ival[15], &ival[16], &ival[17], &ival[18], &ival[19] );
	// #8
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIIII", 
			&ival[20], &ival[21], &ival[22], &ival[23], &ival[24] );
	// #9
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "I", &Mcount);
	// #10-29
		for(i=0; i<20; i++) {
			int j = i * 10;
			nfd->ReadLineEx( buff );
			nfd->ReadRecord( buff, "FFFFFFFFFF", 
				&mval[j  ], &mval[j+1], &mval[j+2], &mval[j+3], &mval[j+4],
				&mval[j+5], &mval[j+6], &mval[j+7], &mval[j+8], &mval[j+9] );
		}
	// #30
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "I", &Fcount);
	// #31-35
		for(i=0; i<5; i++) {
			int j = i * 10;
			nfd->ReadLineEx( buff );
			nfd->ReadRecord( buff, "IIIIIIIIII", 
				&fval[j  ], &fval[j+1], &fval[j+2], &fval[j+3], &fval[j+4],
				&fval[j+5], &fval[j+6], &fval[j+7], &fval[j+8], &fval[j+9] );
		}
	// #36
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "I", &Tcount);
	// #37-43
		int rec_n = Tcount / 10;
		for(i=0; i<rec_n; i++) {
			int j = i * 10;
			nfd->ReadLineEx( buff );
			nfd->ReadRecord( buff, "IIIIIIIIII", 
				&tval[j  ], &tval[j+1], &tval[j+2], &tval[j+3], &tval[j+4],
				&tval[j+5], &tval[j+6], &tval[j+7], &tval[j+8], &tval[j+9] );
		}
	// # ---------------------
		for(i=0; i<FunctionCount; i++){
			cfunc_rec rec;
			// ##1
			nfd->ReadLineEx( buff );
			nfd->ReadRecord( buff, "II", &rec.ID, &rec.type);
			// ##2
			nfd->ReadLineEx( buff );
			nfd->ReadStr( buff, rec.title, sizeof( rec.title ));
			// ## -------------------
			while(1){
				nf_int index;
				nf_float x, y;
				nfd->ReadLineEx( buff );
				nfd->ReadRecord( buff, "IFF", &index, &x, &y );
				if( index == -1 ) break;
				cfunc_rec::cdata_pair data(index,x,y);
				rec.data.push_back(data);
			}
			func_list.push_back( rec );
		}
}






void CNFDB_601::WriteData( CNFData* nfd, FILE* fp )
{
	int i;

	// #1
		nfd->WriteData( fp, "IIIIIIIn", 
			ID,
			format,
			color,
			type,
			subtype,
			layer,
			FunctionCount );
	// #2
		nfd->WriteStr( fp, title );
	// #3
		nfd->WriteData( fp, "In", Bcount);
	// #4
		nfd->WriteData( fp, "IIIIIIIIIIn", 
			bval[0], bval[1], bval[2], bval[3], bval[4],
			bval[5], bval[6], bval[7], bval[8], bval[9] );
	// #5
		nfd->WriteData( fp, "In", Icount);
	// #6
		nfd->WriteData( fp, "IIIIIIIIIIn", 
			ival[0], ival[1], ival[2], ival[3], ival[4],
			ival[5], ival[6], ival[7], ival[8], ival[9] );
	// #7
		nfd->WriteData( fp, "IIIIIIIIIIn", 
			ival[10], ival[11], ival[12], ival[13], ival[14],
			ival[15], ival[16], ival[17], ival[18], ival[19] );
	// #8
		nfd->WriteData( fp, "IIIIIn", 
			ival[20], ival[21], ival[22], ival[23], ival[24] );
	// #9
		nfd->WriteData( fp, "In", Mcount);
	// #10-29
		for(i=0; i<20; i++) {
			int j = i * 10;
			nfd->WriteData( fp, "FFFFFFFFFFn", 
				mval[j  ], mval[j+1], mval[j+2], mval[j+3], mval[j+4],
				mval[j+5], mval[j+6], mval[j+7], mval[j+8], mval[j+9] );
		}
	// #30
		nfd->WriteData( fp, "In", Fcount);
	// #31-35
		for(i=0; i<5; i++) {
			int j = i * 10;
			nfd->WriteData( fp, "IIIIIIIIIIn", 
				fval[j  ], fval[j+1], fval[j+2], fval[j+3], fval[j+4],
				fval[j+5], fval[j+6], fval[j+7], fval[j+8], fval[j+9] );
		}
	// #36
		nfd->WriteData( fp, "In", Tcount);
	// #37-43
		int rec_n = Tcount / 10;
		for(i=0; i<rec_n; i++) {
			int j = i * 10;
			nfd->WriteData( fp, "IIIIIIIIIIn", 
				tval[j  ], tval[j+1], tval[j+2], tval[j+3], tval[j+4],
				tval[j+5], tval[j+6], tval[j+7], tval[j+8], tval[j+9] );
		}
	// # ---------------------

		vector<cfunc_rec>::iterator iter;

		for(iter = func_list.begin(); iter != func_list.end(); iter++ ){
			// ##1
			nfd->WriteData( fp, "IIn", iter->ID, iter->type);
			// ##2
			nfd->WriteStr( fp, iter->title );
			// ## -------------------
			vector<cfunc_rec::cdata_pair>::iterator id;
			for(id = iter->data.begin(); id != iter->data.end(); id++ ){
				nfd->WriteData( fp, "IFFn", id->index, id->x, id->y );
			}
			nfd->WriteData( fp, "IFFn", -1, 0.0, 0.0);
		}
}



