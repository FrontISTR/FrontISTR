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
	CHECDB_Element Ver.1.0
*/

#include <set>
#include "CHECDB.h"
#include "CHECData.h"

using namespace std;


int CHECDB_Element::NodeNumber(int type)
{
	switch(type) {
	case 111: return 2;
	case 112: return 3;
	case 231: return 3;
	case 232: return 6;
	case 241: return 4;
	case 242: return 8;
	case 341: return 4;
	case 342: return 10;
	case 351: return 6;
	case 352: return 15;
	case 361: return 8;
	case 362: return 20;
	case 541: return 8;
	case 542: return 16;
	case 611: return 2;
	case 612: return 3;
	case 731: return 3;
	case 732: return 6;
	case 741: return 4;
	case 742: return 8;
	default:
		return 0;
	}
}


int CHECDB_Element::FaceNumber(int type)
{
	switch(type) {
	case 111: return 0;
	case 112: return 0;
	case 231: return 3;
	case 232: return 3;
	case 241: return 4;
	case 242: return 4;
	case 341: return 4;
	case 342: return 4;
	case 351: return 5;
	case 352: return 5;
	case 361: return 6;
	case 362: return 6;
	case 541: return 2;
	case 542: return 2;
	case 611: return 0;
	case 612: return 0;
	case 731: return 1;
	case 732: return 1;
	case 741: return 1;
	case 742: return 1;
	default:
		return 0;
	}
}


int CHECDB_Element::DOFNumber(int type)
{
	switch(type) {
	case 111: return 1;
	case 112: return 1;
	case 231: return 2;
	case 232: return 2;
	case 241: return 2;
	case 242: return 2;
	case 341: return 3;
	case 342: return 3;
	case 351: return 3;
	case 352: return 3;
	case 361: return 3;
	case 362: return 3;
	case 541: return 3;
	case 542: return 3;
	case 611: return 3;
	case 612: return 3;
	case 731: return 6;
	case 732: return 6;
	case 741: return 6;
	case 742: return 6;
	default:
		return 0;
	}
}


int* Connectivity(int type, int face_id, int& n )
{
	static int tc111[][20] = {{0,1}};
	static int tn111[]     = {2,0};
	static int tc112[][20] = {{0,1,2}};
	static int tn112[]     = {2,0};
	static int tc231[][20] = {{1,2},{2,0},{0,1}};
	static int tn231[]     = {2,2,2,0};
	static int tc232[][20] = {{1,3,2},{2,4,0},{0,5,1}};
	static int tn232[]     = {3,3,3,0};
	static int tc241[][20] = {{3,0},{1,2},{0,1},{2,3}};
	static int tn241[]     = {2,2,2,2,0};
	static int tc242[][20] = {{3,7,0},{1,5,2},{0,4,1},{2,6,3}};
	static int tn242[]     = {3,3,3,3,0};
	static int tc341[][20] = {{1,2,3},{0,3,2},{0,1,3},{0,2,1}};
	static int tn341[]     = {3,3,3,3,0};
	static int tc342[][20] = {{1,4,2,9,3,8},{0,7,3,9,2,5},{0,6,1,8,3,7},{0,5,2,4,1,6}};
	static int tn342[]     = {6,6,6,6,0};
	static int tc351[][20] = {{1,2,5,4},{2,0,3,5},{0,1,4,3},{2,1,0},{3,4,5}};
	static int tn351[]     = {4,4,4,3,3,0};
	static int tc352[][20] = {{1,6,2,14,5,9,4,13},{2,7,0,12,3,10,5,14},{0,8,1,13,4,11,3,12},{2,6,1,8,0,7},{3,11,4,9,5,10}};
	static int tn352[]     = {8,8,8,6,6,0};
	static int tc361[][20] = {{3,0,4,7},{1,2,6,5},{0,1,5,4},{2,3,7,6},{3,2,1,0},{4,5,6,7}};
	static int tn361[]     = {4,4,4,4,4,4,0};
	static int tc362[][20] = {{3,11,0,16,4,15,7,19},{1,9,2,18,6,13,5,17},{0,8,1,17,5,12,4,16},{2,10,3,19,7,14,6,18},{3,10,2,9,1,8,0,11},{4,12,5,13,6,14,7,15}};
	static int tn362[]     = {8,8,8,8,8,8,0};
	static int tc541[][20] = {{3,2,1,0},{4,5,6,7}};
	static int tn541[]     = {4,4,0};
	static int tc542[][20] = {{3,10,2,9,1,8,0,11},{4,12,5,13,6,14,7,15}};
	static int tn542[]     = {8,8,0};
	static int tc611[][20] = {{0,1}};
	static int tn611[]     = {2,0};
	static int tc612[][20] = {{0,1,2}};
	static int tn612[]     = {3,0};
	static int tc731[][20] = {{0,1,2},{2,1,0}};
	static int tn731[]     = {3,3,0};
	static int tc732[][20] = {{0,3,1,4,2,5},{2,4,1,3,0,5}};
	static int tn732[]     = {6,6,0};
	static int tc741[][20] = {{0,1,2,3},{3,2,1,0}};
	static int tn741[]     = {4,4,0};
	static int tc742[][20] = {{0,4,1,5,2,6,3,7},{3,6,2,5,1,4,0,7}};
	static int tn742[]     = {4,4,0};

	#define CODE_GENERATE( x ) \
		case x: n = tn##x [face_id]; return tc##x[face_id];

	switch(type) {
	CODE_GENERATE( 111 )
	CODE_GENERATE( 112 )
	CODE_GENERATE( 231 )
	CODE_GENERATE( 232 )
	CODE_GENERATE( 241 )
	CODE_GENERATE( 242 )
	CODE_GENERATE( 341 )
	CODE_GENERATE( 342 )
	CODE_GENERATE( 351 )
	CODE_GENERATE( 352 )
	CODE_GENERATE( 361 )
	CODE_GENERATE( 362 )
	CODE_GENERATE( 541 )
	CODE_GENERATE( 542 )
	CODE_GENERATE( 611 )
	CODE_GENERATE( 612 )
	CODE_GENERATE( 731 )
	CODE_GENERATE( 732 )
	CODE_GENERATE( 741 )
	CODE_GENERATE( 742 )
	default:
		assert(0);
	}
	return 0;
}



bool CHECDB_Element::CheckType( int type )
{
	const int type_list[] = {
		111,112,231,232,241,242,351,352,361,362,541,542,611,612,731,732,741,742
	};
	int type_n = sizeof( type_list );
	for(int i=0; i<type_n; i++) {
		if( type == type_list[i] ) return true;
	}
	return false;
}


CHECDB_Element::CHECDB_Element()
 : CHECDataBlock(HECDB_ELEMENT), type(0), sec_id(0), option(0), ElemList()
{
}


CHECDB_Element::~CHECDB_Element()
{
	Clear();
}


void CHECDB_Element::Clear()
{
	ElemList.clear();
}


void CHECDB_Element::Write( CHECData* hecd )
{
	if( ElemList.size() == 0 ) return;

	hecd->WriteHeader( "!ELEMENT", "I", "TYPE", type );

	int n = NodeNumber(type);

	set<CElemItem>::iterator iter;
	for(iter = ElemList.begin(); iter != ElemList.end(); iter++){
		hecd->ClearDataLineBuffer();
		hecd->AddDataLineItems( "I", iter->ID );
		for(int i=0; i<n; i++) {
			hecd->AddDataLineItems( "I", iter->node[i] );
		}
		hecd->WriteDataLine();
	}	
}


CHECDB_Element::CElemItem* CHECDB_Element::GetElem( int id )
{
	set<CElemItem>::iterator iter;

	for(iter = ElemList.begin(); iter != ElemList.end(); iter++){
		if( iter->ID == id ) return (CElemItem*)&(*iter);
	}
	return 0;

/*
	iter = find( ElemList.begin(), ElemList.end(), checdb_id_class(id) );
	if( iter == ElemList.end()) return 0;
	return (CElemItem*)&(*iter);
*/
}


bool CHECDB_Element::Read( CHECData* hecd, char* header_line )
{
	int rcode[10];

	hecd->ParseHeader( header_line, rcode, "I", "TYPE", &type );
	int n = NodeNumber(type);
	if(n<=0) return false;

	CElemItem item(type);

	char line[256];
	while( hecd->ReadLine(line) ){
		if( line[0] == '!' ) {
			hecd->PushReadLine(line);
			break;
		}
		char* token = strtok( line, ",\r\n");
		if(!token) return false;
		if( sscanf( token, "%d", &item.ID ) != 1 ) return false;
		for(int i=0; i<n; i++ ){
			token = strtok( 0, ",\r\n");
			if(!token) return false;
			int nid;
			if( sscanf( token, "%d", &nid ) != 1 ) return false;
			item.node[i] = nid;
		}
		ElemList.insert( item );
	}
	return true;
}


