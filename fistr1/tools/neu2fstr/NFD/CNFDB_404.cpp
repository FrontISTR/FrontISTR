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
	CNFDB_404 Ver.1.0
*/


// 404 Element

#include "CNFData.h"
#include "CNFDB_404.h"

using namespace std;


CNFDB_404::CNFDB_404()
 : CNFDataBlock(404)
{}


CNFDB_404::~CNFDB_404()
{
	clear_ref_node_set();
}


void CNFDB_404::Read( CNFData* nfd )
{
	char buff[256];
	int i;

	clear_ref_node_set();

	// #1
		nfd->ReadLineEx( buff );
		if( nfd->version >= 8.0 ) {
			nfd->ReadRecord( buff, "IIIIIIIBIIIII", 
				&ID,
				&color,
				&propID,
				&type,
				&topology,
				&layer,
				&orientID,
				&matl_orflag,
				&geomID,
				&formulation,
				&contactsegment[0],
				&contactsegment[1],
				&formulation2
			);
		} else if( nfd->version >= 6.0 ) {
			nfd->ReadRecord( buff, "IIIIIIIBIIII", 
				&ID,
				&color,
				&propID,
				&type,
				&topology,
				&layer,
				&orientID,
				&matl_orflag,
				&geomID,
				&formulation,
				&contactsegment[0],
				&contactsegment[1]
			);
		} else if( nfd->version >= 5.0 ) {
			nfd->ReadRecord( buff, "IIIIIIIBI", 
				&ID,
				&color,
				&propID,
				&type,
				&topology,
				&layer,
				&orientID,
				&matl_orflag,
				&geomID
			);
		} else {
			nfd->ReadRecord( buff, "IIIIIIIB", 
				&ID,
				&color,
				&propID,
				&type,
				&topology,
				&layer,
				&orientID,
				&matl_orflag
			);
		}
	// #2
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIIIIIIIII", 
			&node[0], &node[1], &node[2], &node[3], &node[4],
			&node[5], &node[6], &node[7], &node[8], &node[9] );
	// #3
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIIIIIIIII", 
			&node[10], &node[11], &node[12], &node[13], &node[14],
			&node[15], &node[16], &node[17], &node[18], &node[19] );
	// #4
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFF", &orient[0], &orient[1], &orient[2]);
	// #5
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFF", &offset1[0], &offset1[1], &offset1[2]);
	// #6
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFF", &offset2[0], &offset2[1], &offset2[2]);
	// #7
		nfd->ReadLineEx( buff );
		if( nfd->version < 4.4 ) {
			nfd->ReadRecord( buff, "BBBBBBBBBBBB", 
				&release1[0], &release1[1], &release1[2],
				&release1[3], &release1[4], &release1[5],
				&release2[0], &release2[1], &release2[2],
				&release2[3], &release2[4], &release2[5] );
		} else {
			nfd->ReadRecord( buff, "BBBBBBBBBBBBIIII", 
				&release1[0], &release1[1], &release1[2],
				&release1[3], &release1[4], &release1[5],
				&release2[0], &release2[1], &release2[2],
				&release2[3], &release2[4], &release2[5],
				&list[0], &list[1], &list[2], &list[3] );
		}

	// Ver.4.5 ========================
	if( nfd->version < 4.5 ) return;

	int list_n = 0;
	for( i=0; i<4; i++) {
		if( list[i] != 0 ) list_n++;
	}
	if( list_n == 0 ) return;

	for(i=0; i<list_n; i++) {
		cref_node_list* node_list = make_ref_node_list( nfd );
		ref_node_set.push_back(node_list);
	}
}



void CNFDB_404::clear_ref_node_set()
{
	vector<cref_node_list*>::iterator iter;
	for(iter = ref_node_set.begin(); iter != ref_node_set.end(); iter++){
		delete *iter;
	}
	ref_node_set.clear();
}


CNFDB_404::cref_node_list* CNFDB_404::make_ref_node_list( CNFData* nfd )
{
	char buff[256];
	cref_node ref;
	cref_node_list* list = new cref_node_list();

	while(1) {
		nfd->ReadLineEx( buff );

		// check of end
		char temp[256];
		strcpy( temp, buff );
		char* token = strtok( temp, ",\n\r");
		int id;
		if(sscanf(token, "%d", &id) != 1) {
			delete list;
			throw CNFError( NFE_INVALID_TOKEN, nfd->line, 1 );
		}
		if( id == -1 ) break;

		nfd->ReadRecord( buff, "IIFIIIIII", 
			&ref.NodeID,
			&ref.faceID,
			&ref.weight,
			&ref.dof[1], &ref.dof[2], &ref.dof[3],
			&ref.dof[4], &ref.dof[5], &ref.dof[6]
		);
		list->ref_node.push_back( ref );
	}
	return list;
}



//*****************************************************************************

void CNFDB_404::WriteData( class CNFData* nfd, FILE* fp )
{

	// #1
		if( nfd->version >= 8.0 ) {
			nfd->WriteData( fp, "IIIIIIBIIIIIn", 
				ID,
				color,
				propID,
				topology,
				layer,
				orientID,
				matl_orflag,
				geomID,
				formulation,
				contactsegment[0],
				contactsegment[1],
				formulation2
			);
		} else if( nfd->version >= 6.0 ) {
			nfd->WriteData( fp, "IIIIIIBIIIIn", 
				ID,
				color,
				propID,
				topology,
				layer,
				orientID,
				matl_orflag,
				geomID,
				formulation,
				contactsegment[0],
				contactsegment[1]
			);
		} else if( nfd->version >= 5.0 ) {
			nfd->WriteData( fp, "IIIIIIBIn", 
				ID,
				color,
				propID,
				topology,
				layer,
				orientID,
				matl_orflag,
				geomID
			);
		} else {
			nfd->WriteData( fp, "IIIIIIBn", 
				ID,
				color,
				propID,
				topology,
				layer,
				orientID,
				matl_orflag
			);
		}
	// #2
		nfd->WriteData( fp, "IIIIIIIIIIn", 
			node[0], node[1], node[2], node[3], node[4],
			node[5], node[6], node[7], node[8], node[9] );
	// #3
		nfd->WriteData( fp, "IIIIIIIIIIn", 
			node[10], node[11], node[12], node[13], node[14],
			node[15], node[16], node[17], node[18], node[19] );
	// #4
		nfd->WriteData( fp, "FFFn", orient[0], orient[1], orient[2]);
	// #5
		nfd->WriteData( fp, "FFFn", offset1[0], offset1[1], offset1[2]);
	// #6
		nfd->WriteData( fp, "FFFn", offset2[0], offset2[1], offset2[2]);
	// #7
		if( nfd->version < 4.4 ) {
			nfd->WriteData( fp, "BBBBBBBBBBBBn", 
				release1[0], release1[1], release1[2],
				release1[3], release1[4], release1[5],
				release2[0], release2[1], release2[2],
				release2[3], release2[4], release2[5] );
		} else {
			nfd->WriteData( fp, "BBBBBBBBBBBBIIIIn", 
				release1[0], release1[1], release1[2],
				release1[3], release1[4], release1[5],
				release2[0], release2[1], release2[2],
				release2[3], release2[4], release2[5],
				list[0], list[1], list[2], list[3] );
		}

	// Ver.4.5 ========================
	if( nfd->version < 4.5 ) return;

	std::vector<cref_node_list*>::iterator iter;
	for(iter = ref_node_set.begin(); iter != ref_node_set.end(); iter++) {
		write_ref_node_list(  nfd, fp, *iter );
	}
}


void CNFDB_404::write_ref_node_list( CNFData* nfd, FILE* fp, CNFDB_404::cref_node_list* list )
{
	std::vector<cref_node>::iterator iter;
	for(iter = list->ref_node.begin(); iter != list->ref_node.end(); iter ++) {
		nfd->WriteData( fp, "IIFIIIIIIn", 
			iter->NodeID,
			iter->faceID,
			iter->weight,
			iter->dof[1], iter->dof[2], iter->dof[3],
			iter->dof[4], iter->dof[5], iter->dof[6]
		);
	}
	nfd->WriteData( fp, "In", -1 ); 
}



