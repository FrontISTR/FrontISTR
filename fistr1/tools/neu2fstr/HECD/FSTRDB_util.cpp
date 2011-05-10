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
	Utility for CFSTRDB Ver.1.0
*/


#include <string.h>
#include <ctype.h>
#include "CFSTRDB.h"
#include "hecd_util.h"

using namespace hecd_util;


CHECDataBlock* CreateFSTRDataBlock( const char* header_name )
{
	char name[80];
	char* np = name;
	char* p = (char*)header_name;

	if( *p == '!' ) p++;
	*np = (char)toupper(*p);
	p++;
	np++;
	while(*p) {
		*np = (char)tolower(*p);
		p++;
		np++;
	}
	*np = 0;

	#define GENERATE_CODE( x ) \
		else if( strcmp( #x , name )==0 ){\
			return new CFSTRDB_##x(); \
		}

	if(false); // dummy
	GENERATE_CODE( Solution )
	GENERATE_CODE( Solver )
	GENERATE_CODE( Write )
	GENERATE_CODE( Echo )
	GENERATE_CODE( Step )
	GENERATE_CODE( Static )
	GENERATE_CODE( Boundary )
	GENERATE_CODE( CLoad )
	GENERATE_CODE( DLoad )
	GENERATE_CODE( Temperature )
	GENERATE_CODE( Reftemp )
	GENERATE_CODE( Eigen )
	GENERATE_CODE( Heat )
	GENERATE_CODE( Fixtemp )
	GENERATE_CODE( CFlux )
	GENERATE_CODE( DFlux )
	GENERATE_CODE( SFlux )
	GENERATE_CODE( Film )
	GENERATE_CODE( SFilm )
	GENERATE_CODE( Radiate )
	GENERATE_CODE( SRadiate )

	#undef 	GENERATE_CODE

	return 0;
}


bool IsFSTRDataBlockName( const char* name )
{
	char s[256];
	if(name[0]=='!')
		toupper( &name[1], s );
	else
		toupper( name, s );
	
	#define GENERATE_CODE( x ) \
		else if( strcmp( #x , name )==0 ){\
			return true;\
		}

	if(false); // dummy
	GENERATE_CODE( SOLUTION )
	GENERATE_CODE( SOLVER )
	GENERATE_CODE( WRITE )
	GENERATE_CODE( ECHO )
	GENERATE_CODE( STEP )
	GENERATE_CODE( STATIC )
	GENERATE_CODE( BOUNDARY )
	GENERATE_CODE( CLOAD )
	GENERATE_CODE( DLOAD )
	GENERATE_CODE( TEMPERATURE )
	GENERATE_CODE( REFTEMP )
	GENERATE_CODE( EIGEN )
	GENERATE_CODE( HEAT )
	GENERATE_CODE( FIXTEMP )
	GENERATE_CODE( CFLUX )
	GENERATE_CODE( DFLUX )
	GENERATE_CODE( SFLUX )
	GENERATE_CODE( FILM )
	GENERATE_CODE( SFILM )
	GENERATE_CODE( RADIATE )
	GENERATE_CODE( SRADIATE )
	#undef 	GENERATE_CODE

	return false;
}



