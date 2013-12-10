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
	CFSTRDB_Solver Ver.1.0
*/

#include "CFSTRDB.h"
#include "CHECData.h"

using namespace std;
using namespace hecd_util;


CFSTRDB_Solver::CFSTRDB_Solver()
 : CFSTRDataBlock(FSTRDB_SOLVER)
{
	Clear();
}


CFSTRDB_Solver::~CFSTRDB_Solver()
{
	Clear();
}

void CFSTRDB_Solver::Clear()
{
	strcpy( method, "CG" );
	precond = 1;
	nset = 0;
	iterlog = 1;
	timelog = 1;
	// 2nd line
	nier = 10000;
	iterPREmax = 2;
	nrest = 10;
	// 3rd line
	resid = 1.0e-8;
	fsigma_diag = 1.0;
	sigma = 0.0;
	// 4th line
	thresh = 0.1;
	filter = 0.1;
}

void CFSTRDB_Solver::Write( CHECData* hecd )
{
	const char YorN[][20] = { "NO", "YES" };

	hecd->WriteHeader( "!SOLVER", "SISS",
			"METHOD",  method,
			"PRECOND", precond,
			"ITERLOG", YorN[iterlog],
			"TIMELOG", YorN[timelog]
	);

	// 2nd line ------------------------

	hecd->WriteData( "III", nier, iterPREmax, nrest );

	// 3rd line ------------------------

	hecd->WriteData( "FFF", resid, fsigma_diag, sigma );

	// 4th line ------------------------

	hecd->WriteData( "FF", thresh, filter );
}


bool CFSTRDB_Solver::Read( CHECData* hecd, char* header_line )
{
	int rcode[10];
	char is[256]="";
	char ts[256]="";
	char iterlog_s[256]="";
	char timelog_s[256]="";

	if(! hecd->ParseHeader( header_line, rcode, "SISS",
		"METHOD",  method, "PRECOND", &precond, "ITERLOG", is, "TIMELOG", ts )) return false;

	cleanup_token( is, iterlog_s ); toupper( iterlog_s );
	cleanup_token( ts, timelog_s ); toupper( timelog_s );

	if( strcmp( iterlog_s, "YES" )==0)     iterlog = 1;
	else if( strcmp( iterlog_s, "NO" )==0) iterlog = 0;
	else if( iterlog_s[0] != 0 ) return false;

	if( strcmp( timelog_s, "YES" )==0)     timelog = 1;
	else if( strcmp( timelog_s, "NO" )==0) timelog = 0;
	else if( timelog_s[0] != 0 ) return false;

	// 2nd line ------------------------

	if(! hecd->ReadData( rcode, "III", &nier, &iterPREmax, &nrest )) return true;

	// 3rd line ------------------------

	if(! hecd->ReadData( rcode, "FFF", &resid, &fsigma_diag, &sigma )) return true;

	// 4th line ------------------------

	if(! hecd->ReadData( rcode, "FF", &thresh, &filter )) return true;

	return true;
}


