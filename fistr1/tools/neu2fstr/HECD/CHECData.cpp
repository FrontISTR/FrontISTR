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
	CHECData Ver.1.0
*/


#include <stdarg.h>
#include <assert.h>
#include "CHECData.h"
#include "hecd_util.h"

using namespace std;
using namespace hecd_util;


//=============================================================================
// common
//=============================================================================


CHECData::CHECData()
{
}


CHECData::~CHECData()
{
	Clear();
}


void CHECData::Clear()
{
	vector<CHECDataBlock*>::iterator iter;
	for(iter = DB.begin(); iter!=DB.end(); iter++){
		delete *iter;
	}
	DB.clear();
}


void CHECData::StoreDataBlock( CHECDataBlock* block )
{
	DB.push_back( block );
}


//=============================================================================
// Save
//=============================================================================


bool CHECData::Save( const char* file_name )
{
	strcpy( fname, file_name);
	fp = fopen( fname, "w");
	if(!fp) return false;

	vector<CHECDataBlock*>::iterator iter;
	for(iter = DB.begin(); iter!=DB.end(); iter++){
		(*iter)->Write( this );
	}

	fclose(fp);
	fp = 0;
	return true;	
}


void CHECData::WriteLine( const char* s )
{
	fprintf(fp, "%s\n", s );
}


// ----------------------------------------------------------------------------
// Header Line
// ----------------------------------------------------------------------------
// fmt : format of parameters.
//       Each character specify the parameter's type
//         'I': integer, 'F':float, 'S':string
// ... : piars of parameter name and value


void CHECData::WriteHeader( const char* name, const char* fmt, ... )
{
	va_list va;
	va_start(va, fmt);

	int n = strlen(fmt);
	char* c = (char*)fmt;

	char line[256];
	char param[256];
	char s[256];

	strcpy( line, name );

	for(int i=0; i<n; i++, c++){
		strcpy( param, va_arg(va, char*));
		strcat(line, ",");
		strcat(line, param);
		strcat(line, "=");
		switch(*c){
		case 'I':
			sprintf(s, "%d", va_arg(va, int));
			break;
		case 'F':
			// sprintf(s, "%lg", va_arg(va, double));
			ftos( va_arg(va, double), s );
			break;
		case 'S':
			sprintf(s, "%s", va_arg(va, char*));
			break;
		default:
			assert(0);
		}
		strcat(line, s);
	}
	fprintf( fp, "%s\n", line );
	va_end(va);
}


// ----------------------------------------------------------------------------
// Parameter Line
// ----------------------------------------------------------------------------

void CHECData::WriteParameter(const char* fmt, ... )
{
	va_list va;
	va_start(va, fmt);

	int n = strlen(fmt);
	char* c = (char*)fmt;

	char line[256];
	char param[256];
	char s[256];

	strcpy( line, "!" );

	for(int i=0; i<n; i++, c++){
		strcpy( param, va_arg(va, char*));
		if( i!=0 )
			strcat(line, ",");
		strcat(line, param);
		strcat(line, "=");
		switch(*c){
		case 'I':
			sprintf(s, "%d", va_arg(va, int));
			break;
		case 'F':
			//sprintf(s, "%lg", va_arg(va, double));
			ftos( va_arg(va, double), s );
			break;
		case 'S':
			sprintf(s, "%s", va_arg(va, char*));
			break;
		default:
			assert(0);
		}
		strcat(line, s);
	}
	fprintf( fp, "%s\n", line );
	va_end(va);
}


// ----------------------------------------------------------------------------
// Data Line
// ----------------------------------------------------------------------------
// fmt : formt of data like one of WriteHeader
// ... : values to write the file.


void CHECData::WriteData( const char* fmt, ...)
{
	va_list va;
	va_start(va, fmt);

	int n = strlen(fmt);
	char* c = (char*)fmt;

	char line[256];
	char s[256];

	line[0] = 0;

	for(int i=0; i<n; i++, c++){
		if( i!=0 ) strcat(line, ",");
		switch(*c){
		case 'I':
			sprintf(s, "%d", va_arg(va, int));
			break;
		case 'F':
			//sprintf(s, "%lg", va_arg(va, double));
			ftos( va_arg(va, double), s );
			break;
		case 'S':
			sprintf(s, "%s", va_arg(va, char*));
			break;
		default:
			assert(0);
		}
		strcat(line, s);
	}
	fprintf( fp, "%s\n", line );
	va_end(va);
}


void CHECData::ClearDataLineBuffer()
{
	data_line_buffer[0] = 0;
}


void CHECData::AddDataLineItems( const char* fmt, ...)
{
	va_list va;
	va_start(va, fmt);

	int n = strlen(fmt);
	char* c = (char*)fmt;

	char s[256];

	for(int i=0; i<n; i++, c++){
		switch(*c){
		case 'I':
			sprintf(s, "%d", va_arg(va, int));
			break;
		case 'F':
			//sprintf(s, "%lg", va_arg(va, double));
			ftos( va_arg(va, double), s );
			break;
		case 'S':
			sprintf(s, "%s", va_arg(va, char*));
			break;
		default:
			assert(0);
		}
		strcat(data_line_buffer, s);
		strcat(data_line_buffer, ",");
	}
	va_end(va);
}


void CHECData::WriteDataLine()
{
	int len = strlen( data_line_buffer );
	char *p = &data_line_buffer[len-1];

	// remove ',' of tail
	while( p>=data_line_buffer ){
		if( *p == ',' ){
			*p = 0;
			break;
		} else if( *p == ' ' ){
			*p = 0;
		} else {
			break;
		}
		p--;
	}

	fprintf(fp, "%s\n", data_line_buffer );
	ClearDataLineBuffer();
}


//=============================================================================
// Load
//=============================================================================


bool CHECData::Load( const char* file_name )
{
	Clear();
	return AddLoad( file_name );
}


bool CHECData::AddLoad( const char* file_name )
{
	line_count = 0;
	fg_header_pushed = false;

	strcpy( fname, file_name);
	fp = fopen( fname, "r");
	if(!fp) return false;

	char line[256];
	char header_name[256];
	CHECDataBlock* block;

	while( ReadLine(line) ) {
		if(!GetHeaderName( line, header_name )) {
			fclose(fp);
			return false;
		}
		if( strcmp( header_name, "!END" )==0 ) break;
		block = CreateDataBlock( header_name );
		if(!block) {
			fclose(fp);
			return false;
		}
		if(!block->Read( this, line )){
			fclose(fp);
			return false;
		}
		DB.push_back(block);
	}
	
	fclose(fp);
	fp = 0;
	return true;	
}


bool CHECData::ReadLine( char* s, int size )
{
	line_count++;

	if(fg_header_pushed){
		strcpy( s, header_line_buffer );
		fg_header_pushed = false;
		return true;
	}

	while( fgets( s, size, fp ) ){
		if( s[0] == 0 || s[0]=='#'|| s[0]=='\r' || s[0]=='\n' ) continue;
		if( s[0] == '!' && s[1] =='!' ) continue;
		remove_cr( s );
		return true;
	}

	return false;
}


void CHECData::PushReadLine( const char* s )
{
	strcpy( header_line_buffer, s );
	fg_header_pushed = true;
	line_count--;
}


CHECDataBlock* CHECData::CreateDataBlock( const char* header_name )
{
	return CreateHECDataBlock( header_name );
}


bool CHECData::GetHeaderName( const char* header_line, char* header_name )
{
	#define is_separator(x) (x==',' || x==' ' || x=='\t' || x=='\r' || x=='\n' )

	char* p = (char*)header_line;
	while(*p && is_separator(*p)) p++;
	if( *p != '!' ) return false;

	char* bp = header_name;
	while(*p && !is_separator(*p)) {
		*bp = (char)toupper(*p);
		p++;
		bp++;
	}

	*bp = 0;
	return true;

	#undef is_separator
}


//-----------------------------------------------------------------------------
// ParseHeader
//-----------------------------------------------------------------------------
// rcode[i] : 1 -- set, 0 -- not set, -1 -- error
// fmt : string composed by following characters
//	'I' : int
//	'F' : double
//	'S' : char*
//	'E' : int (exist or not)
// ... : pairs of a parameter name and a pointer of the parameter
// example)
//	char* line = "!SOLVER, TYPE=CG, PRECOND=1, ITERLOG=YES";
//	int rcode[3];
//	char type[30]; int precond; char iterlog[30];
//	ParseHeader(line,rcode,"SIS","TYPE",type,"PRECOND",&precond,"ITERLOG",iterlog);
//-----------------------------------------------------------------------------


static
bool get_param_token( char*& p, char* token )
{
	#define is_separator(x) (x==',' || x==' ' || x=='\t' || x=='\r' || x=='\n' )

	while(*p && is_separator(*p)) p++;
	if( !*p ) return false;

	char* t = token;
	if( *p == '=' ) {
		t[0] = '=';
		t[1] = 0;
		p++;
		return true;
	}

	while(*p && is_separator(*p)) p++;
	if( !*p ) return false;

	while( *p && !is_separator(*p) && *p!='=') {
		*t = *p;
		t++;
		p++;
	}
	*t = 0;
	return true;
}


bool CHECData::vParseParameter( char* line, int* rcode, const char* fmt, va_list va )
{
	const int max_param_n = 40;
	char param_name[max_param_n][20];
	void* param[max_param_n];
	int param_n = strlen(fmt);
	int i;
	int* ip;

	for(i=0; i<param_n; i++) {
		strcpy( param_name[i], va_arg(va, char*));
		param[i] = va_arg(va, void*);
		rcode[i] = 0;
		if(fmt[i] == 'E' ){
			ip = (int*)param[i];
			*ip = 0;
		}
	}

	char* p = line;
	char p_token[50];
	char eq_token[50];
	char v_token[50];
	char p_str[50];

	while( get_param_token( p, p_token )) {
		cleanup_token( p_token, p_str );
		toupper( p_str );
		for( i=0; i<param_n; i++) {
			if( strcmp( p_str, param_name[i] ) == 0 ) break;
		}
		if( i==param_n ) return false;
		// ---------------------------------
		if( fmt[i] == 'E' ) {
			rcode[i] = 1;
			ip = (int*)param[i];
			*ip = 1;
			continue;
		}
		// ---------------------------------
		if(! get_param_token( p, eq_token ) || eq_token[0] != '=') {
			rcode[i] = -1;
			return false;
		}
		if(! get_param_token( p, v_token )) {
			rcode[i] = -1;
			return false;
		}
		switch( fmt[i] ){
		case 'S':
			if( sscanf( v_token, "%s", (char*)(param[i]) )!= 1) {
				rcode[i] = -1;
				return false;
			}
			break;
		case 'I':
			if( sscanf( v_token, "%d", (int*)(param[i]) )!= 1) {
				rcode[i] = -1;
				return false;
			}
			break;
		case 'F':
			if( sscanf( v_token, "%lf", (double*)(param[i]) )!= 1) {
				rcode[i] = -1;
				return false;
			}
			break;
		default:
			assert(0);
		}
		rcode[i] = 1;
	}
	return true;
}


bool CHECData::ParseParameter( char* line, int* rcode, const char* fmt, ...)
{
	va_list va;
	va_start(va, fmt);
	bool fg = vParseParameter( line, rcode, fmt, va );
	va_end(va);
	return fg;
}


bool CHECData::ParseHeader( char* header_line, int* rcode, const char* fmt, ... )
{
	char* p = header_line;
	while(*p) {
		if(*p == ',' || *p == '\r' || *p == '\n' ) {
			break;
		}
		p++;
	}

	va_list va;
	va_start(va, fmt);
	bool fg = vParseParameter( p, rcode, fmt, va );
	va_end(va);
	return fg;
}


bool CHECData::ReadParameter( int* rcode, const char* fmt, ... )
{
	char line[256];
	if(!ReadLine(line)) return false;

	char* p = line;
	if( *p == '!' ) p++;
	va_list va;
	va_start(va, fmt);
	bool fg = vParseParameter( p, rcode, fmt, va );
	va_end(va);
	return fg;
}


//-----------------------------------------------------------------------------
// ReadData
//-----------------------------------------------------------------------------
// fmt : string composed by 'I', 'F' or 'S'
// ... : pointers of parameter
//-----------------------------------------------------------------------------

bool CHECData::ReadData( int* rcode, const char* fmt, ... )
{
	char line[256];

	if(!ReadLine(line)) return false;
	if(line[0] == '!') {
		PushReadLine(line);
		return false;
	}

	const int max_param_n = 100;
	int i;
	int param_n = strlen(fmt);
	void* param[max_param_n];
	va_list va;
	va_start(va, fmt);
	for(i=0; i<param_n; i++) {
		param[i] = va_arg( va, void* );
		rcode[i] = 0;
	}
	va_end(va);

	char* token = strtok( line, ",\r\n" );
	i = 0;
	while( token && i<param_n ){
		switch( fmt[i] ){
		case 'I':
			if( sscanf( token, "%d", (int*)(param[i]) ) != 1 ) return false;
			break;
		case 'F':
			if( sscanf( token, "%lf", (double*)(param[i]) ) != 1 ) return false;
			break;
		case 'S':
			if( sscanf( token, "%s", (char*)(param[i]) ) != 1 ) return false;
			break;
		default:
			assert(0);
		}
		rcode[i] = 1;
		token = strtok( 0, ",\r\n" );
		i++;
	}

	return true;
}


int CHECData::ParseDoubleDataArray( char* line, double* data )
{
	int n = 0;
	char* token = strtok(line,",\r\n");
	while(token){
		if( sscanf( token, "%lf", &data[n])!=1 ) return -(n+1);
		n++;
		token = strtok(0, ",\r\n");
	}
	return n;
}

int CHECData::ParseIntDataArray( char* line, int* data )
{
	int n = 0;
	char* token = strtok(line,",\r\n");
	while(token){
		if( sscanf( token, "%d", &data[n])!=1 ) return -(n+1);
		n++;
		token = strtok(0, ",\r\n");
	}
	return n;
}


//=============================================================================
// Utility
//=============================================================================


template <class T>
T hecdata_get( CHECData* hd, int data_type, const char* name )
{
	vector<CHECDataBlock*>::iterator iter;
	for(iter = hd->DB.begin(); iter != hd->DB.end(); iter++) {
		if( (*iter)->data_type == data_type ){
			T block = (T)(*iter);
			if(strcmp(block->name, name) == 0 ){
				return block;
			}
		}
	}
	return 0;
}



CHECDB_Material* CHECData::GetMaterial( const char* name )
{
	return hecdata_get<CHECDB_Material*>( this, HECDB_MATERIAL, name );
}


CHECDB_NGroup* CHECData::GetNGroup( const char* name )
{
	return hecdata_get<CHECDB_NGroup*>( this, HECDB_NGROUP, name );
}


CHECDB_EGroup* CHECData::GetEGroup( const char* name )
{
	return hecdata_get<CHECDB_EGroup*>( this, HECDB_EGROUP, name );
}


CHECDB_SGroup* CHECData::GetSGroup( const char* name )
{
	return hecdata_get<CHECDB_SGroup*>( this, HECDB_SGROUP, name );
}



CHECDB_Node::CNodeItem* CHECData::GetNodeItem( int id )
{
	vector<CHECDataBlock*>::iterator iter;
	for(iter = DB.begin(); iter != DB.end(); iter++) {
		if( (*iter)->data_type != HECDB_NODE) continue;
		CHECDB_Node::CNodeItem* nitem = ((CHECDB_Node*)(*iter))->GetNode(id);
		if(nitem) return nitem;
	}
	return 0;
}


CHECDB_Element::CElemItem* CHECData::GetElemItem( int id )
{
	vector<CHECDataBlock*>::iterator iter;
	for(iter = DB.begin(); iter != DB.end(); iter++) {
		if( (*iter)->data_type != HECDB_ELEMENT) continue;
		CHECDB_Element::CElemItem* eitem = ((CHECDB_Element*)(*iter))->GetElem(id);
		if(eitem) return eitem;
	}
	return 0;
}


int CHECData::GetElemType( int id )
{
	vector<CHECDataBlock*>::iterator iter;
	for(iter = DB.begin(); iter != DB.end(); iter++) {
		if( (*iter)->data_type != HECDB_ELEMENT) continue;
		CHECDB_Element::CElemItem* eitem = ((CHECDB_Element*)(*iter))->GetElem(id);
		if(eitem) return ((CHECDB_Element*)(*iter))->type;
	}
	return 0;
}


