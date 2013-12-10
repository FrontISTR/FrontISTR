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


#include <iostream>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "CConvMessage.h"
#include "CNFData.h"
#include "CFSTRData.h"
#include "conv_neu2hec.h"
#include "conv_neu2fstr_static.h"
#include "conv_neu2fstr_dynamic.h"
#include "conv_neu2fstr_heat.h"


#define VER	"1.004"
#define AUTHOR	"The University of Tokyo, RSS21 Project"

static int solution = sol_static;
static bool fg_hecmw_ctrl = false;
static bool fg_res_name = false;
static bool fg_vis_name = false;
static char date_time[256];
static char neu_name[256];
static char mesh_name[256];
static char ctrl_name[256];
static char hecmw_ctrl_name[256] = "hecmw_ctrl.dat";
static char res_name[256]; // default: mesh_name + ".res"
static char vis_name[256]; // default: mesh_name + ".vis"
static bool fg_direct = false;

static CNFData nfdata;
static CFSTRData fstrdata;


using namespace std;


static
void remove_cr( char* s )
{
	int n = strlen(s);
	if( n==0 ) return;
	n--;
	if( s[n] == '\r' || s[n] == '\n' ) s[n] = 0;
}


//-----------------------------------------------------------------------------

void title()
{
	printf(	"NEU to FrontSTR Data Converter, Ver. %s \n", VER );
	printf(	"Copyright (C) 2006 %s, All right reserved.\n", AUTHOR );
}

//-----------------------------------------------------------------------------

void help()
{
	printf( "[usage]\n"
		"  neu2fstr ([options]) [sSeEhH] [neu] [mesh] [ctrl]\n"
		"  [options]\n"
		"     -h        : show help (this message)\n"
		"     -c        : create HEC-MW control file (hecmw_ctrl.dat)\n"
		"     -r [name] : write result file name in hecmw_ctrl.dat (default:[mesh].res)\n"
		"     -v [name] : write visual file name in hecmw_ctrl.dat (default:[mesh].vis)\n"
		"     -d        : set direct solver in FrontSTR control file\n"
		"  [sSeEhH]  : command character s,S,e,E,h or H (specify only one character)\n"
		"     s,S    : generate for static analysis\n"
		"     e,E    : generate for eigen analysis\n"
		"              CAUTION) Edit !EIGEN section in generated [ctrl]\n"
		"     h,H    : generate for heat transfer analysis\n"
		"  [neu]     : neutral file name\n"
		"  [mesh]    : mesh file name\n"
		"  [ctrl]    : FrontSTR control file name\n"
	);
	fflush(stdout);
}

//-----------------------------------------------------------------------------

bool set_params( int argc, char** argv )
{
	int count = 0;

	#define CMP(x) (strcmp( argv[i], (x)) == 0)

	for( int i=1; i<argc; i++) {
		if( CMP( "-h" )){
			return false;
		} else if( CMP( "-c" )){
			fg_hecmw_ctrl = true;
		} else if( CMP( "-d" )){
			fg_direct = true;
		} else if( CMP( "-v" )) {
			fg_vis_name = true;
			i++;
			if( i>=argc ) {
				fprintf(stderr,"##Error: Visual output file name is required\n");
				return false;
			}
			strcpy( vis_name, argv[i] );
		} else if( CMP( "-r" )) {
			fg_res_name = true;
			i++;
			if( i>=argc ) {
				fprintf(stderr,"##Error: Result file name is required\n");
				return false;
			}
			strcpy( res_name, argv[i] );
		} else {
			switch( count ){
			case 0:
				if( CMP("s") || CMP("S")){
					solution = sol_static;
				} else if(CMP("e") || CMP("E")){
					solution = sol_eigen;
				} else if(CMP("h") || CMP("H")){
					solution = sol_heat;
				} else {
					fprintf( stderr, "##Error: No such solution type %s\n", argv[i] );
					return false;
				}
			case 1:
				strcpy( neu_name, argv[i] );
				break;
			case 2:
				strcpy( mesh_name, argv[i] );
				break;
			case 3:
				strcpy( ctrl_name, argv[i] );
				break;
			default:
				fprintf( stderr, "##Error: Too many arguments\n");
				return false;
			}
			count++;
		}
	}

	#undef CMP

	switch( count ){
	case 0:
		fprintf(stderr, "##Error: one character, neu, mesh and ctrl file must be specified\n");
		return false;
	case 1:
		fprintf(stderr, "##Error: neu, mesh and ctrl file must be specified\n");
		return false;
	case 2:
		fprintf(stderr, "##Error: mesh and ctrl file must be specified\n");
		return false;
	case 3:
		fprintf(stderr, "##Error: ctrl file must be specified\n");
		return false;
	}

	if( !fg_res_name ) {
		strcpy( res_name, mesh_name );
		strcat( res_name, ".res" );
	}
	if( !fg_vis_name ) {
		strcpy( vis_name, mesh_name );
		strcat( vis_name, ".vis" );
	}
	return true;
}

//-----------------------------------------------------------------------------

void create_comment( char* s, bool fg_mesh )
{
	char ss[256];
	s[0] = 0;

	strcat(s, "###########################################################\n");
	if(fg_mesh){
		sprintf(ss, " HEC-MW Mesh File Generated by neu2fstr ver.%s\n", VER );
	} else {
		sprintf(ss, " FrontSTR Control File Generated by neu2fstr ver.%s\n", VER );
	}
	strcat(s,ss);
	sprintf(ss, " Date&Time: %s\n", date_time );
	strcat(s,ss);
	sprintf(ss, " Original : %s\n", neu_name );
	strcat(s,ss);
	if(fg_mesh){
		sprintf(ss, " Control  : %s\n", ctrl_name );
	} else {
		sprintf(ss, " Mesh     : %s\n", mesh_name );
	}
	strcat(s,ss);
	strcat(s, "###########################################################\n");
}

//-----------------------------------------------------------------------------

bool create_hecmw_ctrl()
{
	FILE* fp = fopen(hecmw_ctrl_name, "w");
	if(!fp) {
		fprintf( stderr, "##Error: Cannot create hecmw_ctrl.dat file\n");
		return false;
	}
	fprintf( fp, "###########################################################\n");
	fprintf( fp, "# HEC-MW Control File Generated by neu2fstr ver.%s\n", VER );
	fprintf( fp, "# Date&Time: %s\n", date_time );
	fprintf( fp, "# Original : %s\n", neu_name );
	fprintf( fp, "###########################################################\n");
	fprintf( fp, "!MESH, NAME=fstrMSH,TYPE=HECMW-ENTIRE\n");
	fprintf( fp, "  %s\n", mesh_name );
	fprintf( fp, "!CONTROL,NAME=fstrCNT\n");
	fprintf( fp, "  %s\n", ctrl_name );
	fprintf( fp, "!RESULT,NAME=fstrRES,IO=OUT\n");
	fprintf( fp, "  %s\n", res_name );
	fprintf( fp, "!RESULT,NAME=vis_out,IO=OUT\n");
	fprintf( fp, "  %s\n", vis_name );
	fclose(fp);

	return true;
}

//-----------------------------------------------------------------------------

void set_solution( CFSTRData& fstrdata, int solution )
{
	CFSTRDB_Solution* sol = new CFSTRDB_Solution();

	switch( solution ){
	case sol_static:
		sol->type = CFSTRDB_Solution::TYPE_STATIC;
		break;
	case sol_eigen:
		sol->type = CFSTRDB_Solution::TYPE_EIGEN;
		break;
	case sol_heat:
		sol->type = CFSTRDB_Solution::TYPE_HEAT;
		break;
	default:
		assert(0);
	}
	fstrdata.DB.push_back(sol);
}


//-----------------------------------------------------------------------------


void set_solver( CFSTRData& fstrdata, bool fg_direct = false )
{
	CFSTRDB_Solver* solver = new CFSTRDB_Solver();
	if( fg_direct ){
		strcpy( solver->method, "DIRECT");
	} else {
		strcpy( solver->method, "CG");
	}
	fstrdata.DB.push_back( solver );
}


//-----------------------------------------------------------------------------
// main
//-----------------------------------------------------------------------------

int main( int argc, char** argv )
{
	time_t t;
	time(&t);
	strcpy( date_time, ctime(&t));
	remove_cr( date_time );

	title();
	if(!set_params( argc, argv )) {
		fflush(stderr);
		help();
		return -1;
	}

	cout << "loading neu file ... " << endl;

	try {
		nfdata.Load( neu_name );
	} catch( CNFError e ){
		fprintf( stderr, "NEU loading error : %s\n", e.Msg());
		return -1;
	}

	try {
		cout << "converting to HEC-MW mesh ... " << endl;
		conv_neu2hec( nfdata, fstrdata, solution );
		cout << "converting to FrontSTR control data " << endl;
		conv_neu2fstr_static( nfdata, fstrdata );
		conv_neu2fstr_dynamic( nfdata, fstrdata );
		conv_neu2fstr_heat( nfdata, fstrdata );
		set_solution( fstrdata, solution );
		set_solver( fstrdata, fg_direct );
	} catch( CConvMessage& e ){
		fprintf( stderr, "Converting error : %s\n", e.Msg());
		return -1;
	}

	char comment[256];

	cout << "saving HEC-MW mesh file..." << endl;
	create_comment( comment, true );
	if( !fstrdata.SaveMesh( mesh_name, comment )) {
		fprintf( stderr, "Cannot save mesh data\n");
		return -1;
	}	

	cout << "saving FrontSTR control file..." << endl;
	create_comment( comment, false );
	if( !fstrdata.SaveCtrl( ctrl_name, comment )) {
		fprintf( stderr, "Cannot save control data\n");
		return -1;
	}	

	if( fg_hecmw_ctrl ){
		cout << "creating hecmw_ctrl.dat file..." << endl;
		if(!create_hecmw_ctrl()) return -1;
	}
	cout << "neu2fstr completed." << endl;
	return 0;
}

