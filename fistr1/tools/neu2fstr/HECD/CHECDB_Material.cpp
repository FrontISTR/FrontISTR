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
	CHECDB_Material Ver.1.0
*/

#include "CHECDB.h"
#include "CHECData.h"

using namespace std;

CHECDB_Material::CHECDB_Material()
 : CHECDataBlock( HECDB_MATERIAL ), ItemList()
{
	name[0] = 0;
}


CHECDB_Material::~CHECDB_Material()
{
	Clear();
}


void CHECDB_Material::Clear()
{
	ItemList.clear();
	name[0] = 0;
}



void CHECDB_Material::CItem::Write( CHECData* hecd )
{
	hecd->WriteParameter("II", "ITEM", ID, "SUBITEM", SubItemNumber() );

	vector<CItemRec>::iterator ri;
	for(ri = RecList.begin(); ri != RecList.end(); ri++) {
		hecd->ClearDataLineBuffer();
		vector<double>::iterator pi;
		for( pi = ri->params.begin(); pi != ri->params.end(); pi++){
			hecd->AddDataLineItems( "F", *pi );
		}
		hecd->WriteDataLine();
	}
}


void CHECDB_Material::Write( CHECData* hecd )
{
	hecd->WriteHeader( "!MATERIAL", "SI", "NAME", name, "ITEM", ItemList.size());

	vector<CItem>::iterator iter;
	for( iter = ItemList.begin(); iter != ItemList.end(); iter++){
		iter->Write( hecd );
	}
}

//-----------------------------------------------------------------------------


bool CHECDB_Material::CItem::Read( CHECData* hecd )
{
	int subitem_n = 0;
	int rcode[10];
	int i, n;
	char line[256];
	const int max_data_n = 100;
	double data[max_data_n];

	if(!hecd->ReadParameter( rcode, "II", "ITEM", &ID, "SUBITEM", &subitem_n)) return false;
	while(1) {
		if( !hecd->ReadLine( line )) break;
		if( line[0] == '!' ) {
			hecd->PushReadLine(line);
			break;
		}
		n = hecd->ParseDoubleDataArray( line, data );
		if(n<0) return false;
		CItemRec rec;
		if( n == subitem_n ) {
			rec.last_is_temp = false;
		} else if( n == subitem_n+1 ) {
			rec.last_is_temp = true;
		} else {
			return false;
		}
		for(i=0; i<n; i++) {
			rec.params.push_back( data[i] );
		}
		RecList.push_back( rec );
	}
	return true;
}


bool CHECDB_Material::Read( CHECData* hecd, char* header_line )
{
	int item_n = 0;
	int rcode[10];
	if(!hecd->ParseHeader( header_line, rcode, "SI", "NAME", name, "ITEM", &item_n )) return false;

	for(int i=0; i<item_n; i++) {
		CItem item;
		if( !item.Read( hecd )) return false;
		ItemList.push_back( item );
	}
	return true;
}

