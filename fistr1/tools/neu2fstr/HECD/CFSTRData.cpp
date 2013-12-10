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
	CFSTRData Ver.1.0
*/

#include <vector>
#include "CFSTRData.h"
#include "CFSTRDB.h"

using namespace std;

CFSTRData::CFSTRData()
 : CHECData()
{
}


//=============================================================================
// Save
//=============================================================================


bool CFSTRData::SaveMesh( const char* file_name, const char* comment )
{
	strcpy( fname, file_name);
	fp = fopen( fname, "w");
	if(!fp) return false;

	WriteComment(fp, comment);

	vector<CHECDataBlock*>::iterator iter;
	for(iter = DB.begin(); iter!=DB.end(); iter++){
		if(!(*iter)->IsMesh()) continue;
		if( (*iter)->data_type == HECDB_VISUAL ) continue;
		(*iter)->Write( this );
	}
	WriteHeader( "!END" );
	fclose(fp);
	fp = 0;
	return true;	
}



bool CFSTRData::SaveCtrl( const char* file_name, const char* comment )
{
	strcpy( fname, file_name);
	fp = fopen( fname, "w");
	if(!fp) return false;

	WriteComment(fp, comment);

	CHECDataBlock* vis = 0;

	vector<CHECDataBlock*>::iterator iter;
	for(iter = DB.begin(); iter!=DB.end(); iter++){
		if((*iter)->data_type == HECDB_VISUAL ) {
			vis = *iter;
			continue;
		}
		if((*iter)->IsMesh()) continue;
		(*iter)->Write( this );
	}
	if( vis ){
		vis->Write( this );
	}
	WriteHeader( "!END" );
	fclose(fp);
	fp = 0;
	return true;	
}


void CFSTRData::WriteComment( FILE* fp, const char* comment)
{
	int n = strlen(comment);
	if( n==0 ) return;

	char* buff = new char[n+1];

	strcpy( buff, comment );
	char* p = strtok( buff, "\r\n");
	while(p){
		if( *p == '#' ){
			fprintf(fp, "%s\n", p );
		} else {
			fprintf(fp, "#%s\n", p );
		}
		p = strtok( 0, "\r\n");
	}

	delete[] buff;
}



//=============================================================================
// Load
//=============================================================================


CHECDataBlock* CFSTRData::CreateDataBlock( const char* header_name )
{
	CHECDataBlock* block = CreateHECDataBlock( header_name );
	if( block ) return block;
	return CreateFSTRDataBlock( header_name );
}
