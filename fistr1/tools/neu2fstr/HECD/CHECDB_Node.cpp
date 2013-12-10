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
	CHECDB_Node Ver.1.0
*/

#include "CHECDB.h"
#include "CHECData.h"

using namespace std;


CHECDB_Node::CHECDB_Node()
 : CHECDataBlock(HECDB_NODE), NodeList()
{
}

CHECDB_Node::~CHECDB_Node()
{
	Clear();
}

void CHECDB_Node::Clear()
{
	NodeList.clear();
}


void CHECDB_Node::Write( CHECData* hecd )
{
	if( NodeList.size() == 0 ) return;

	hecd->WriteHeader( "!NODE" );

	set<CNodeItem>::iterator iter;
	for(iter = NodeList.begin(); iter != NodeList.end(); iter++){
		hecd->WriteData( "IFFF", iter->ID, iter->x, iter->y, iter->z );
	}	
}


CHECDB_Node::CNodeItem* CHECDB_Node::GetNode( int id )
{
	std::set<CNodeItem>::iterator iter;
	for(iter = NodeList.begin(); iter != NodeList.end(); iter++){
		if( iter->ID == id ){
			return (CNodeItem*)&(*iter);
		}
	}
	return 0;
/*
	iter = find( NodeList.begin(), NodeList.end(), CNodeItem(id));
	if( iter == NodeList.end()) return 0;
	return &(*iter);
*/
}


bool CHECDB_Node::Read( CHECData* hecd, char* header_line )
{
	CNodeItem item;

	char line[256];
	while( hecd->ReadLine(line) ){
		if( line[0] == '!' ) {
			hecd->PushReadLine(line);
			break;
		}
		char* token = strtok( line, ",\r\n");
		if(!token) return false;
		if( sscanf( token, "%d", &item.ID ) != 1 ) return false;
		// x
		token = strtok( 0, ",\r\n");
		if(!token) return false;
		if( sscanf( token, "%lf", &item.x ) != 1 ) return false;
		// y
		token = strtok( 0, ",\r\n");
		if(!token) return false;
		if( sscanf( token, "%lf", &item.y ) != 1 ) return false;
		// z
		item.z = 0;
		token = strtok( 0, ",\r\n");
		if( token ) {
			if( sscanf( token, "%lf", &item.z ) != 1 ) return false;
		}

		NodeList.insert( item );
	}
	return true;
}

