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
	Utility for CHECDB Ver.1.0
*/


#include <string.h>
#include <ctype.h>
#include "CHECDB.h"
#include "hecd_util.h"

using namespace hecd_util;


CHECDataBlock* CreateHECDataBlock( const char* header_name )
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
			return new CHECDB_##x (); \
		}

	if(false); // dummy
	GENERATE_CODE( Header )
	GENERATE_CODE( Node )
	GENERATE_CODE( Element )
	GENERATE_CODE( Material )
	GENERATE_CODE( Section )
	GENERATE_CODE( NGroup )
	GENERATE_CODE( EGroup )
	GENERATE_CODE( SGroup )
	GENERATE_CODE( Amplitude )
	GENERATE_CODE( Zero )
	GENERATE_CODE( Visual )

	#undef 	GENERATE_CODE

	return 0;
}


bool IsHECDataBlockName( const char* name )
{
	char s[256];
	if(name[0]=='!')
		toupper( &name[1], s );
	else
		toupper( name, s );
	
	#define GENERATE_CODE( x ) \
		else if( strcmp( #x , name )==0 ){\
			return true; \
		}

	if(false); // dummy
	GENERATE_CODE( HEADER )
	GENERATE_CODE( NODE )
	GENERATE_CODE( ELEMENT )
	GENERATE_CODE( MATERIAL )
	GENERATE_CODE( SECTION )
	GENERATE_CODE( NGROUP )
	GENERATE_CODE( EGROUP )
	GENERATE_CODE( SGROUP )
	GENERATE_CODE( AMPLITUDE )
	GENERATE_CODE( ZERO )
	GENERATE_CODE( VISUAL )

	#undef 	GENERATE_CODE

	return false;
}

