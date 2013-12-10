/*=====================================================================*
 *                                                                     *
 *   Software Name : neu_reporter                                      *
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
	neu_reporter ver.1.0
	-----------------------
	reporting neutral file
*/


#include <iostream>
#include <stdio.h>
#include "CNFData.h"

using namespace std;


char run_name[256];
CNFData data;
FILE* wfp;


//=============================================================================
// etc.
//=============================================================================

static
void get_fname_from_path( char* fname, const char* path )
{
	int n = strlen(path);
	char* p = (char*)&path[n-1];
	while( path < p && (*p != '\\') && (*p != '/')) p--;
	if((*p == '\\') || (*p == '/')) p++;
	strcpy( fname, p );
}


//-----------------------------------------------------------------------------

// return 0 : quit
static
int input_block_id()
{
	char buff[256];
	int block_id;

	while(1) {
		printf("block# ? (h:show selectable block id, e:exit) > ");
		fflush(stdout);
		fgets(buff, sizeof(buff), stdin);
		char c = buff[0];
		if( c == 'h' ){
			for(int i=0; i<NFD_SupportedBlockListSize; i++) {
				fprintf(wfp, "%8d\n", NFD_SupportedBlockList[i] );
			}
		} else if( c == 'e' ){
			return 0;
		} else {
			if( sscanf(buff, "%d", &block_id ) == 1 ){
				bool fg = false;
				for(int i=0; i<NFD_SupportedBlockListSize; i++) {
					if( NFD_SupportedBlockList[i] == block_id ) {
						fg = true;
						break;
					}
				}
				if( fg ) break;
				printf("Not supported block id\n");
			} else {
				printf("Cannot converting number\n");
			}
		}
	}

	return block_id;
}


//=============================================================================
// command record & list
//=============================================================================


bool fg_quit = false;

class ccmd_rec {
public:
	char name[256];
	char help[256];
	void (*cmd)();
	ccmd_rec();
	ccmd_rec(const char* n ); // for key
	ccmd_rec(const char* n, const char* h, void (*f)() );
	ccmd_rec( const ccmd_rec& rec);
};

set<ccmd_rec> cmd_list;

ccmd_rec::ccmd_rec()
 : cmd(0)
{
	name[0] = 0;
	help[0] = 0;
}


ccmd_rec::ccmd_rec(const char* n )
 : cmd(0)
{
	strcpy(name, n);
	help[0] = 0;
}

ccmd_rec::ccmd_rec(const char* n, const char* h, void (*f)() )
{
	strcpy(name, n);
	strcpy(help, h);
	cmd = f;
}

ccmd_rec::ccmd_rec( const ccmd_rec& rec)
{
	strcpy(name, rec.name);
	strcpy(help, rec.help);
	cmd = rec.cmd;
}


inline
bool operator==(const ccmd_rec& a, const ccmd_rec& b)
{
	return( strcmp( a.name, b.name ) == 0);
}
inline
bool operator<(const ccmd_rec& a, const ccmd_rec& b)
{
	int c = strcmp( a.name, b.name );
	return c<0;
}
inline
bool operator>(const ccmd_rec& a, const ccmd_rec& b)
{
	int c = strcmp( a.name, b.name );
	return c>0;
}


//=============================================================================
// commands
//=============================================================================


#define str_cmd_help	"h"
#define help_cmd_help	"show command list"
void cmd_help()
{
	printf("\n");
	set<ccmd_rec>::iterator iter;
	for( iter = cmd_list.begin(); iter != cmd_list.end(); iter++) {
		printf( "  %8s  : %s\n", iter->name, iter->help );
	}
	printf("\n");
	fflush(stdout);
}

//-----------------------------------------------------------------------------


#define str_cmd_quit	"quit"
#define help_cmd_quit	"quit"
void cmd_quit()
{
	fg_quit = true;
}


//-----------------------------------------------------------------------------


#define str_cmd_save	"save"
#define help_cmd_save	"save"
void cmd_save()
{
	char fname[256];
	printf("save file name >> ");
	fflush(stdout);
	fgets( fname, sizeof(fname), stdin);
	fname[ strlen(fname)-1 ] = 0;

	try {
		data.Save(fname);
	} catch(...){
		printf("error in saving..\n");
		fflush(stdout);
	}
}


//-----------------------------------------------------------------------------


#define str_cmd_open_outfile	"open"
#define help_cmd_open_outfile	"open output file"
void cmd_open_outfile()
{
	char fname[256];
	printf("output file? >");
	fflush(stdout);
	fgets(fname, sizeof(fname), stdin );
	fname[ strlen(fname)-1 ] = 0; // remove CR/LF

	FILE* fp = fopen( fname, "w" );
	if(!fp) {
		cout << "Cannot open file" << endl;
		return;
	}

	if( wfp != stdout && wfp != stderr ){
		fclose(wfp);
		printf("previous output file closed.\n");
	}
	wfp = fp;

	printf("%s is opened.\n", fname);
}


//-----------------------------------------------------------------------------

#define str_cmd_close_outfile	"close"
#define help_cmd_close_outfile	"close output file"
void cmd_close_outfile()
{
	if( wfp != stdout && wfp != stderr ){
		fclose(wfp);
		printf("output file closed.\n");
		wfp = stdout;
	}
}



//-----------------------------------------------------------------------------


#define str_cmd_write_summary	"s"
#define help_cmd_write_summary	"write summary of loaded neutral data"
void cmd_write_summary()
{
	char fname[256];
	get_fname_from_path( fname, data.neu_file );

	fprintf( wfp, "summary of %s\n", fname );
	data.WriteSummary(wfp);
}


//-----------------------------------------------------------------------------

#define str_cmd_write_block	"b"
#define help_cmd_write_block	"write data block"
void cmd_write_block()
{
	int block_id = input_block_id();
	if( block_id == 0 ) return;
	data.WriteDataBlock( wfp, block_id );
}



//=============================================================================
// commands registration & execution
//=============================================================================


void regist_commands()
{

	#define GENERATE_CODE( x ) {\
		ccmd_rec rec( str_##x, help_##x, x );\
		cmd_list.insert(rec);\
	}
			
	GENERATE_CODE( cmd_help )
	GENERATE_CODE( cmd_quit )
	GENERATE_CODE( cmd_save )
	GENERATE_CODE( cmd_open_outfile )
	GENERATE_CODE( cmd_close_outfile )
	GENERATE_CODE( cmd_write_summary )
	GENERATE_CODE( cmd_write_block )

	#undef 	GENERATE_CODE
}


//-----------------------------------------------------------------------------


bool execute_command( const char* name )
{
	ccmd_rec key(name);

	set<ccmd_rec>::iterator iter;
	iter = find( cmd_list.begin(), cmd_list.end(), key);
	if( iter == cmd_list.end()){
		printf( "not such command (h:help)\n");
		return false;
	}
	(*iter->cmd)();
	return true;
}


//-----------------------------------------------------------------------------


void command_line()
{
	char buff[256];

	regist_commands();
	cmd_help();
	do {
		printf( "neu_reporter>");
		fflush(stdout);
		fgets( buff, sizeof(buff), stdin);
		buff[ strlen(buff)-1] = 0; // remove CR/LF
		execute_command( buff );
		fflush(wfp);
	} while( !fg_quit );
	cmd_close_outfile();
	printf("end of neu reporter\n");
}



//=============================================================================
// main & etc.
//=============================================================================

void set_run_name( char* argv0 )
{
	get_fname_from_path( run_name, argv0 );
}


void help()
{
	cout << "Neutral File Reporter Ver.1.0"<< endl;
	cout << "Copyright (C) 2005 Noboru Imai"<< endl;
	cout << "[usage] " << run_name << " [NUE file]" << endl;
}


int main( int argc, char** argv )
{
	set_run_name(argv[0]);
	wfp = stdout;

	if( argc <= 1 ) {
		help();
		return -1;
	}

	try {
		data.Load( argv[1] );
	} catch( CNFError e ){
		cout << e.Msg() << endl;
		return -1;
	}

	cout << "command line start..." << endl;

	command_line();

	return 0;
}
