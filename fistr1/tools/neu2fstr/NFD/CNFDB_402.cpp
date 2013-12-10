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
	CNFDB_402 Ver.1.0
	-----------------------------
	402 Properties ( of element )
*/



#include "CNFData.h"
#include "CNFDB_402.h"

// 402 Properties ( of element )


CNFDB_402::CNFDB_402()
 : CNFDataBlock(402)
{
	num_lam = 0;
	lam_MID = 0;
	num_val = 0;
	Value = 0;
	num_outline = 0;
	u = 0;
	v = 0;
	draw = 0;
	num_outline_2 = 0;
	u_2 = 0;
	v_2 = 0;
	draw_2 = 0;
}


CNFDB_402::~CNFDB_402()
{
	Clear();
}


void CNFDB_402::Clear()
{
	delete[] lam_MID;
	delete[] Value;
	delete[] u;
	delete[] v;
	delete[] draw;
	delete[] u_2;
	delete[] v_2;
	delete[] draw_2;

	num_lam = 0;
	lam_MID = 0;
	num_val = 0;
	Value = 0;
	num_outline = 0;
	u = 0;
	v = 0;
	draw = 0;
	num_outline_2 = 0;
	u_2 = 0;
	v_2 = 0;
	draw_2 = 0;
}



void CNFDB_402::Read( CNFData* nfd )
{
	char buff[256];
	int i;

	Clear();

	// #1
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIIIII", 
				&ID,
				&color,
				&matID,
				&type,
				&layer,
				&refCS
		);
	// #2
		nfd->ReadLineEx( buff );
		nfd->ReadStr( buff, title, sizeof(title));
	// #3
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIII", &floag[0], &floag[1], &floag[2] , &floag[3] );
	// #4
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "I", &num_lam );
	// # ----------------------
	//	8 values par record;
		lam_MID = new nf_int[num_lam];
		nfd->ReadMultRec( 'I', 8, num_lam, lam_MID );
	// # ----------------------
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "I", &num_val );
	// # ----------------------
	//	5 values par record;
		Value = new nf_float[num_val];
		nfd->ReadMultRec( 'F', 5, num_val, Value );

	// ======= Ver.6.0 ========================
	if( nfd->version < 6.0 ) return;

	// # ----------------------
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "I", &num_outline );
	// # ----------------------
		u = new nf_float[num_outline];
		v = new nf_float[num_outline];
		draw = new nf_int[num_outline];
		for(i=0; i<num_outline; i++ ){
			nfd->ReadLineEx( buff );
			nfd->ReadRecord( buff, "FFI", &u[i], &v[i], &draw[i] );
		}

	// ======= Ver.8.1 ========================
	if( nfd->version < 8.1 ) return;

	// # ----------------------
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "I", &num_outline_2 );
	// # ----------------------
		u_2 = new nf_float[num_outline_2];
		v_2 = new nf_float[num_outline_2];
		draw_2 = new nf_int[num_outline_2];
		for(i=0; i<num_outline_2; i++ ){
			nfd->ReadLineEx( buff );
			nfd->ReadRecord( buff, "FFI", &u_2[i], &v_2[i], &draw_2[i] );
		}
}


//*****************************************************************************


void CNFDB_402::WriteData( class CNFData* nfd, FILE* fp )
{
	int i;

	// #1
		nfd->WriteData( fp, "IIIIIIn", 
				ID,
				color,
				matID,
				type,
				layer,
				refCS
		);
	// #2
		nfd->WriteStr( fp , title );
	// #3
		nfd->WriteData( fp, "IIIIn", floag[0], floag[1], floag[2] , floag[3] );
	// #4
		nfd->WriteData( fp, "In", num_lam );
	// # ----------------------
	//	8 values par record;
		for( i=0; i<num_lam; i++) {
			nfd->WriteData(fp, "I", lam_MID[i] );
			if( i%8 == 7 ) nfd->WriteData( fp, "n" );
		}
		if( i%8 != 0 ) nfd->WriteData( fp, "n" );
	// # ----------------------
		nfd->WriteData( fp, "In", num_val );
	// # ----------------------
	//	5 values par record;
		for( i=0; i<num_val; i++) {
			nfd->WriteData(fp, "F", Value[i] );
			if( i%5 == 4 ) nfd->WriteData( fp, "n" );
		}
		if( i%5 != 0 ) nfd->WriteData( fp, "n" );

	// ======= Ver.6.0 ========================
	if( nfd->version < 6.0 ) return;

	// # ----------------------
		nfd->WriteData( fp, "In", num_outline );
	// # ----------------------
		for(i=0; i<num_outline; i++ ){
			nfd->WriteData( fp, "FFIn", u[i], v[i], draw[i] );
		}

	// ======= Ver.8.1 ========================
	if( nfd->version < 8.1 ) return;

	// # ----------------------
		nfd->WriteData( fp, "In", num_outline_2 );
	// # ----------------------
		for(i=0; i<num_outline_2; i++ ){
			nfd->WriteData( fp, "FFIn", u_2[i], v_2[i], draw_2[i] );
		}
}

