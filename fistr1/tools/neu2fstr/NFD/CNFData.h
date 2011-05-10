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
	CNFData Ver.1.0
*/


#ifndef CNFDataH
#define CNFDataH

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif


#include <set>
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "CNFMessage.h"
#include "CNFDataBlock.h"
#include "CNFDB_100.h"
#include "CNFDB_402.h"
#include "CNFDB_403.h"
#include "CNFDB_404.h"
#include "CNFDB_405.h"
#include "CNFDB_408.h"
#include "CNFDB_506.h"
#include "CNFDB_507.h"
#include "CNFDB_601.h"

const double DefaultCNFDataVersion = 8.0;

// results of ReadLine method
const int READLINE_SUCESS    =  1;
const int READLINE_SEPARATOR = -1;
const int READLINE_EOF       = -2;
const int READLINE_ERROR     =  0;


const int NFD_SupportedBlockListSize = 9;
const int NFD_SupportedBlockList[] = {
	100, 402, 403, 404, 405, 408, 506, 507, 601 };

class CNFData {
public:
	double version;		// set by instance of CNFDB_100
	char title[256];	// set by instance of CNFDB_100

	CNFData();
	CNFData( const char* fname );
	virtual ~CNFData();

	virtual void Clear();
	virtual void Load( const char* fname );
	virtual void Save( const char* fname );

public:
	//-----------------------------------------------------
	// Indivisual DataBlock Strage and Utilities
	// CAUTION) No strage for CNFDB_100
	//-----------------------------------------------------

	#define GENRATE_CODE( x ) \
		std::vector< CNFDB_##x *> DB_##x;\
		void Clear_##x ();

	GENRATE_CODE( 402 )
	GENRATE_CODE( 403 )
	GENRATE_CODE( 404 )
	GENRATE_CODE( 405 )
	GENRATE_CODE( 408 )
	GENRATE_CODE( 506 )
	GENRATE_CODE( 507 )
	GENRATE_CODE( 601 )
	#undef GENRATE_CODE

public:
	//-----------------------------------------------------
	// DataBlock Utilities
	//-----------------------------------------------------
	virtual CNFDataBlock* CreateDataBlock( int block_id );
	virtual void StoreDataBlock( CNFDataBlock* block );
	virtual void SkipDataBlock();
	//-----------------------------------------------------
	// Basic Input/Output Methods
	//-----------------------------------------------------
	char neu_file[256];
	FILE* log_fp;
	FILE* fp;
	int line;
	char line_buff[512];
	int fg_line_buff_empty;
	// fmt: I:integer, B:bool, F:float
	int EndOfFile() { return feof(fp); }
	void ReadRecStart( char* buff );
	int ReadRecNext( char type, void* value );
	void ReadRecord( char* buff, const char* fmt, ... );
	int ReadLine( char* buff, int size=255 );
	void ReadLineEx( char* buff, int size=255 );
	void ReadStr( char* buff, char* s, int size );
	void PushBackLine( char* buff );
	void ReadMultRec( char type, int n_in_rec, int val_n, void* val );

	void WriteStr( FILE* fp, const char* s );
	void WriteData( FILE* fp, const char* fmt, ... );
	void WriteBlockSeparator( FILE* fp);

	bool WriteDataBlock( FILE* fp, int id );

	//-----------------------------------------------------
	// Report Status
	//-----------------------------------------------------
	void WriteSummary( FILE* fp=0 );
protected:
	char rec_buff[256];
	int fg_rec_first;
	int rec_column;

	virtual void PrintMessage( const char* msg );

protected:
	std::set<int> non_supported_block_list;
};


#endif
