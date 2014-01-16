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

#include <stdarg.h>

#include "CNFData.h"

using namespace std;


CNFData::CNFData()
 : version(DefaultCNFDataVersion),
   fp(0), line(0)
{
	neu_file[0] = 0;
	title[0] = 0;
	log_fp = stdout;
}

CNFData::CNFData( const char* fname )
 : version(DefaultCNFDataVersion), fp(0), line(0)
{
	title[0] = 0;
	log_fp = stdout;
	Load(fname);
}


CNFData::~CNFData()
{
	Clear();
}


void CNFData::Clear()
{
	#define GENRATE_CODE( x ) Clear_##x ();

	GENRATE_CODE( 402 )
	GENRATE_CODE( 403 )
	GENRATE_CODE( 404 )
	GENRATE_CODE( 405 )
	GENRATE_CODE( 408 )
	GENRATE_CODE( 506 )
	GENRATE_CODE( 507 )
	GENRATE_CODE( 601 )

	#undef GENRATE_CODE

	non_supported_block_list.clear();
}


#define GENRATE_CODE( x ) \
	void CNFData::Clear_##x () {\
		vector< CNFDB_##x *>::iterator iter;\
		for( iter = DB_##x.begin(); iter != DB_##x.end(); iter++ ) \
			delete *iter;\
		DB_##x.clear(); \
	}

GENRATE_CODE( 402 )
GENRATE_CODE( 403 )
GENRATE_CODE( 404 )
GENRATE_CODE( 405 )
GENRATE_CODE( 408 )
GENRATE_CODE( 506 )
GENRATE_CODE( 507 )
GENRATE_CODE( 601 )

#undef GENRATE_CODE



//*****************************************************************************
// Load
//*****************************************************************************


void CNFData::Load( const char* fname )
{
	char buff[256];
	nf_int BlockID;
	CNFDataBlock *block;
	int fg_block_continue;

	strcpy( neu_file, fname );

	fp = fopen( fname, "r");
	if(!fp) {
		char s[256];
		sprintf(s, ": %s", fname );
		throw CNFError(NFE_OPEN_ERROR , s);
	}

	non_supported_block_list.clear();
	fg_line_buff_empty = TRUE;
	line = 0;
	try {
		while(!EndOfFile()){
			int rcode = ReadLine( buff );
			if( rcode == READLINE_EOF ) {
				break;
			} else if( rcode == READLINE_ERROR ) {
				throw CNFError(NFE_READDATA_ERROR, line);
			} else if( rcode != READLINE_SEPARATOR ) {
				throw CNFError(NFE_DATA_BLOCK_REQUIRED, line);
			}
			ReadLineEx( buff );
			ReadRecord( buff, "I", &BlockID );
			do {
				block = CreateDataBlock(BlockID);
				if( block ){
					fprintf(log_fp,"line:%d reading block %d\n", (int)line, (int)BlockID);
					if( ReadLine(buff) == READLINE_SUCESS ){
						PushBackLine(buff);
						block->Read( this );
						StoreDataBlock( block );
						if( ReadLine(buff) != READLINE_SEPARATOR ){
							PushBackLine(buff);
							fg_block_continue = TRUE;
						} else {
							fg_block_continue = FALSE;
						}
					} else {
						fprintf(log_fp,"   ---- blank block\n");
						delete block;
						fg_block_continue = FALSE;
					}
				} else {
					non_supported_block_list.insert( BlockID );
					char s[256];
					sprintf(s,": %d", (int)BlockID);
					CNFWarning w( NFW_NON_SUPPORTED_DATA_BLOCK, s, line );
					PrintMessage( w.Msg());
					SkipDataBlock();
					fg_block_continue = FALSE;
				}
			} while( fg_block_continue );
		}
	} catch( CNFError e ){
		fclose(fp);
		throw e;
	}

	fclose(fp);
	fprintf(log_fp, "Data reading completed!\n");
}


//*****************************************************************************
// Save
//*****************************************************************************



void CNFData::Save( const char* fname )
{
	strcpy( neu_file, fname );

	fp = fopen( fname, "w");
	if(!fp) {
		char s[256];
		sprintf(s, ": %s", fname );
		throw CNFError(NFE_OPEN_ERROR , s);
	}

	#define GENRATE_CODE( x )\
	if(!CNFData::WriteDataBlock( fp, x )) { fclose(fp); throw CNFError( NFE_WRITEDATA_ERROR ); }

	GENRATE_CODE( 100 )
	GENRATE_CODE( 402 )
	GENRATE_CODE( 403 )
	GENRATE_CODE( 404 )
	GENRATE_CODE( 405 )
	GENRATE_CODE( 408 )
	GENRATE_CODE( 506 )
	GENRATE_CODE( 507 )
	GENRATE_CODE( 601 )

	#undef GENRATE_CODE

	fclose(fp);
}



//*****************************************************************************
// DataBlock Utilities
//*****************************************************************************


CNFDataBlock* CNFData::CreateDataBlock( int BlockID )
{
	#define GENRATE_CODE( x ) case x: return new CNFDB_##x;

	switch( BlockID ){
	GENRATE_CODE( 100 )
	GENRATE_CODE( 402 )
	GENRATE_CODE( 403 )
	GENRATE_CODE( 404 )
	GENRATE_CODE( 405 )
	GENRATE_CODE( 408 )
	GENRATE_CODE( 506 )
	GENRATE_CODE( 507 )
	GENRATE_CODE( 601 )
	default:
		return NULL;
	}

	#undef GENRATE_CODE
}


void CNFData::StoreDataBlock( CNFDataBlock* block )
{
	#define GENRATE_CODE( x ) case x: DB_##x.push_back((CNFDB_##x *)block); break;

	switch( block->DataBlockID ){
	case 100:
		version = ((CNFDB_100*)block)->version;
		strcpy( title, ((CNFDB_100*)block)->title);
		delete block;
		break;
	GENRATE_CODE( 402 )
	GENRATE_CODE( 403 )
	GENRATE_CODE( 404 )
	GENRATE_CODE( 405 )
	GENRATE_CODE( 408 )
	GENRATE_CODE( 506 )
	GENRATE_CODE( 507 )
	GENRATE_CODE( 601 )
	default:
		assert(0);
	}

	#undef GENRATE_CODE
}



void CNFData::SkipDataBlock()
{
	char buff[256];
	while(!EndOfFile()){
		if( ReadLine( buff )==READLINE_SEPARATOR ) return;
	}
}


//*****************************************************************************
// Utilities
//*****************************************************************************


void CNFData::PrintMessage( const char* msg )
{
	fprintf(stdout, "%s\n", msg );
	if( stdout != log_fp ) fprintf(log_fp, "%s\n", msg );
	fflush(stdout);
}



//*****************************************************************************
// Basic Input Method
//*****************************************************************************


void CNFData::ReadRecStart( char* buff )
{
	strcpy( rec_buff, buff );
	fg_rec_first = TRUE;
	rec_column = 0;
}


int CNFData::ReadRecNext( char type, void* value )
{
	rec_column++;
	char* token;

	if( fg_rec_first ){
		token = strtok( rec_buff, ",\n\r");
		fg_rec_first = false;
	} else {
		token = strtok( NULL, ",\n\r");
	}
	if(!token) return FALSE;

	switch(type){
	case 'I':
		if(sscanf( token, "%d", (int*)value) != 1){
			throw CNFError(NFE_INVALID_TOKEN, line, rec_column);
		}
		break;
	case 'B':
		if(sscanf( token, "%hhd", (unsigned char*)value) != 1){
			throw CNFError(NFE_INVALID_TOKEN, line, rec_column);
		}
		break;
	case 'F':
		if(sscanf( token, "%lf", (double*)value) != 1){
			throw CNFError(NFE_INVALID_TOKEN, line, rec_column);
		}
		break;
	default:
		assert(0);
	}

	return TRUE;
}


void CNFData::ReadRecord( char* buff, const char* fmt, ... )
{
	va_list va;
	va_start(va, fmt);
	int n = strlen(fmt);
	void* val_ptr[200];
	assert(!(n>200));
	for(int i=0; i<n; i++){
		val_ptr[i] = (void*)va_arg(va, void*);
	}
	va_end(va);

	try {
		ReadRecStart( buff );
		char* c = (char*)fmt;
		for(int i=0; i<n; i++, c++){
			if(!ReadRecNext( *c, val_ptr[i] )){
				throw CNFError( NFE_ITEM_REQUIRED, line, i+1);
			}
		}
	} catch( CNFError err ) {
		throw err;
	}
}


int CNFData::ReadLine( char* buff, int size )
{
	if( !fg_line_buff_empty ){
		strcpy( buff, line_buff );
		fg_line_buff_empty = TRUE;
		return READLINE_SUCESS;
	}
	do {
		line ++;
		if( !fgets( buff, size, fp )) return  READLINE_EOF;
		char* p = buff;
		while(*p) {
			if(*p == '\r' || *p == '\n' ) {
				*p = 0;
				break;
			}
			p++;
		}
	} while( buff[0] == 0 );
	if( strcmp( buff, "   -1") == 0 ){
		return READLINE_SEPARATOR;
	}
	return READLINE_SUCESS;
}


void CNFData::ReadLineEx( char* buff, int size )
{
	if( ReadLine( buff, size ) != READLINE_SUCESS ) {
		throw CNFError(NFE_LINE_REQUIRED, line);
	}
}



void CNFData::ReadStr( char* buff, char* s, int size )
{
	int i;

	if(strcmp(buff, "<NULL>")==0){
		s[0] = 0;
		return;
	}

	int n = strlen(buff);
	assert(!(n >= size));
	for( i=0; i<n; i++) s[i] = buff[i];
	s[i] = 0;
}


void CNFData::PushBackLine( char* buff )
{
	strcpy( line_buff, buff );
	fg_line_buff_empty = FALSE;
}



void CNFData::ReadMultRec( char type, int n_in_rec, int val_n, void* val )
{
	char buff[256];
	int i,j,n,id, end_id;
	char* p = (char*)val;
	int size;

	switch( type ){
	case 'I': size = sizeof(nf_int); break;
	case 'F': size = sizeof(nf_float); break;
	case 'B': size = sizeof(nf_bool); break;
	default:
		assert(0);
	}

	n = val_n / n_in_rec;
	id = 0;
	for(i=0; i<n; i++) {
		ReadLineEx( buff );
		ReadRecStart( buff );
		j = 0;
		end_id = id + n_in_rec;
		for(; id<end_id; id++, j++){
			if(!ReadRecNext( type, p )){
				throw CNFError( NFE_ITEM_REQUIRED, line, j+1);
			}
			p += size;
		}
	}
	if( id < val_n ){
		ReadLineEx( buff );
		ReadRecStart( buff );
		j=0;
		for(; id<val_n; id++, j++){
			if(!ReadRecNext( type, p )){
				throw CNFError( NFE_ITEM_REQUIRED, line, j+1);
			}
			p += size;
		}
	}
}

//-----------------------------------------------------------------------------


void CNFData::WriteStr( FILE* fp, const char* s )
{
	if( s==0 || s[0] == 0) {
		fprintf(fp, "<NULL>\n");
	} else {
		fprintf(fp, "%s\n", s);
	}
}


static
void float_write(  FILE* fp, double x )
{
	fprintf(fp, "%lg,", x);

	/*******************************************
	if( x != 0.0 && -1e-4 <= x && x <= 1e-4 ){
		fprintf(fp, "%le,", x );
		return;
	}

	char buff[256];
	sprintf(buff, "%lf", x);

	char* p = &buff[strlen(buff)-1];
	while( buff < p && *p == '0') p--;
	char* zp = p;
	zp++;
	while( buff < p && *p != '.') p--;

	if( *p == '.' ) {
		*zp = 0;
	}

	fprintf(fp, "%s,", buff);
	********************************************/
}


void CNFData::WriteData( FILE* fp, const char* fmt, ... )
{
	va_list va;
	va_start(va, fmt);
	int n = strlen(fmt);
	int b;
	double x;

	try {
		char* c = (char*)fmt;
		for(int i=0; i<n; i++, c++){
			switch( *c ) {
			case 'n': case 'N':
				fprintf(fp,"\n");
				break;
			case 'i': case 'I':
				fprintf(fp,"%d,", (int)(va_arg(va, nf_int)));
				break;
			case 'b': case 'B':
				b = va_arg(va, nf_int);
				fprintf(fp,"%d,", b);
				break;
			case 'f': case 'F':
				x = va_arg(va, nf_float);
				float_write( fp, x );
				break;
			default:
				assert(0);
			}
		}
	} catch( CNFError err ) {
		va_end(va);
		throw err;
	}
	va_end(va);
}



void CNFData::WriteBlockSeparator( FILE* fp)
{
	fprintf(fp, "   -1\n" );
}


bool CNFData::WriteDataBlock( FILE* fp, int id )
{
	#define GENRATE_CODE( x ) \
		case (x) : {\
			if( DB_##x.size() > 0 ) {\
				WriteBlockSeparator( fp );\
				fprintf(fp, "   %d\n", (x));\
				vector< CNFDB_##x *>::iterator iter;\
				for(iter = DB_##x.begin(); iter != DB_##x.end(); iter++ )\
					(*iter)->WriteData( this, fp );\
				WriteBlockSeparator( fp );\
			}\
			break;\
	}

	switch( id ){
	case 100: {
			CNFDB_100 DB100;
			DB100.version = version;
			strcpy( DB100.title, title );
			WriteBlockSeparator( fp );
			fprintf(fp, "   100\n" );\
			DB100.WriteData( this, fp );
			WriteBlockSeparator( fp );
		}
		break;
	GENRATE_CODE( 402 )
	GENRATE_CODE( 403 )
	GENRATE_CODE( 404 )
	GENRATE_CODE( 405 )
	GENRATE_CODE( 408 )
	GENRATE_CODE( 506 )
	GENRATE_CODE( 507 )
	GENRATE_CODE( 601 )
	default:
		return false;
	}

	#undef GENRATE_CODE

	return true;
}


//*****************************************************************************
// Report Status
//*****************************************************************************

void CNFData::WriteSummary( FILE* fp )
{
	if( !fp ) fp = stdout;

	fprintf(fp, "%8s  %8s\n", "BlockID", "num" );
	fprintf(fp, "====================\n");
	#define GENRATE_CODE( x ) \
		if((DB_##x).size() > 0 ) fprintf(fp, "%8d  %8lu\n", (x), (DB_##x).size());

	GENRATE_CODE( 402 )
	GENRATE_CODE( 403 )
	GENRATE_CODE( 404 )
	GENRATE_CODE( 405 )
	GENRATE_CODE( 408 )
	GENRATE_CODE( 506 )
	GENRATE_CODE( 507 )
	GENRATE_CODE( 601 )
	#undef GENRATE_CODE
	fprintf(fp, "====================\n\n");

	fprintf(fp,"non supporting block (skipped)\n");
	int i=0;
	set<int>::iterator iter;
	for(iter = non_supported_block_list.begin(); iter != non_supported_block_list.end(); iter++) {
		if(i==10) {
			fprintf(fp,"\n");
			i = 0;
		}
		fprintf(fp, "  %4d,", *iter );
		i++;
	}
	fprintf(fp,"\n");

	fflush(fp);
}



