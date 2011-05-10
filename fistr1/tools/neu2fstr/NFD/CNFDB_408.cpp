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
	CNFDB_408 Ver.1.0
*/


// 408 Group


#include <vector>
#include "CNFData.h"
#include "CNFDB_408.h"

using namespace std;


CNFDB_408::CNFDB_408()
 : CNFDataBlock(408)
{}

void CNFDB_408::Read( CNFData* nfd )
{
	char buff[256];
	int i;

	// #1
		nfd->ReadLineEx( buff );
		if( nfd->version >= 5.0 ) {
			nfd->ReadRecord( buff, "IIB", 
				&ID,
				&need_eval,
				&prev_enum );
		} else {
			nfd->ReadRecord( buff, "II", 
				&ID,
				&need_eval );
		}
	// #2
		nfd->ReadLineEx( buff );
		nfd->ReadStr( buff, title, sizeof(title));
	// #3
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "III", 
			&layer[0],
			&layer[1],
			&layer_method );
	// #4
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIIIFF", 
			&coclip_on,
			&coclip_dof,
			&coclip_meth,
			&coclip_csys,
			&coclip_min,
			&coclip_max );
	// #5
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "II", 
			&plclip_meth,
			&plclip_in );
	// #--------------------------------

		for( i=0; i<6; i++){
			// ##1
			nfd->ReadLineEx( buff );
			nfd->ReadRecord( buff, "II", 
				&plclip_set[i].plclip_on,
				&plclip_set[i].plclip_neg );
			// ##2
			nfd->ReadLineEx( buff );
			nfd->ReadRecord( buff, "FFF", 
				&plclip_set[i].plclip_base[0],
				&plclip_set[i].plclip_base[1],
				&plclip_set[i].plclip_base[2] );
			// ##3
			nfd->ReadLineEx( buff );
			nfd->ReadRecord( buff, "FFF", 
				&plclip_set[i].plclip_norm[0],
				&plclip_set[i].plclip_norm[1],
				&plclip_set[i].plclip_norm[2] );
		}

	// #--------------------------------
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "I", &max_rules );
	// #--------------------------------
		read_rule( nfd, rule_set );
	// #--------------------------------
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "I", &max_lists );
	// #--------------------------------
		read_list( nfd, list_set );
}


void CNFDB_408::read_rule( CNFData* nfd, std::vector<crule_rec>& list )
{
	char buff[256];
	crule_rec rec;
	list.clear();
	while(1) {
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "I", &rec.rule_type );
		if( rec.rule_type == -1 ) break;

		rec.entry_set.clear();
		crule_rec::centry_rec erec;
		while(1) {
			nfd->ReadLineEx( buff );
			nfd->ReadRecord( buff, "IIII",
				&erec.startID,
				&erec.stopID,
				&erec.incID,
				&erec.include );
			if( erec.startID == -1 ) break;
			rec.entry_set.push_back( erec );
		}
		list.push_back( rec );
	}
}


void CNFDB_408::read_list( CNFData* nfd, std::vector<clist_rec>& list )
{
	char buff[256];
	clist_rec rec;
	list.clear();
	while(1) {
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "I", &rec.list_type );
		if( rec.list_type == -1 ) break;

		rec.entityID.clear();
		nf_int eid;
		while(1) {
			nfd->ReadLineEx( buff );
			nfd->ReadRecord( buff, "I", &eid );
			if( eid == -1 ) break;
			rec.entityID.push_back( eid );
		}
		list.push_back( rec );
	}
}


//*****************************************************************************



void CNFDB_408::WriteData( class CNFData* nfd, FILE* fp )
{
	int i;

	// #1
		if( nfd->version >= 5.0 ) {
			nfd->WriteData( fp, "IIBn", 
				ID,
				need_eval,
				prev_enum );
		} else {
			nfd->WriteData( fp, "IIn", 
				ID,
				need_eval );
		}
	// #2
		nfd->WriteStr( fp, title );
	// #3
		nfd->WriteData( fp, "IIIn", 
			layer[0],
			layer[1],
			layer_method );
	// #4
		nfd->WriteData( fp, "IIIIFFn", 
			coclip_on,
			coclip_dof,
			coclip_meth,
			coclip_csys,
			coclip_min,
			coclip_max );
	// #5
		nfd->WriteData( fp, "IIn", 
			plclip_meth,
			plclip_in );
	// #--------------------------------

		for( i=0; i<6; i++){
			// ##1
			nfd->WriteData( fp, "IIn", 
				plclip_set[i].plclip_on,
				plclip_set[i].plclip_neg );
			// ##2
			nfd->WriteData( fp, "FFFn", 
				plclip_set[i].plclip_base[0],
				plclip_set[i].plclip_base[1],
				plclip_set[i].plclip_base[2] );
			// ##3
			nfd->WriteData( fp, "FFFn", 
				plclip_set[i].plclip_norm[0],
				plclip_set[i].plclip_norm[1],
				plclip_set[i].plclip_norm[2] );
		}

	// #--------------------------------
		nfd->WriteData( fp, "In", max_rules );
	// #--------------------------------
		write_rule( nfd, fp, rule_set );
	// #--------------------------------
		nfd->WriteData( fp, "In", max_lists );
	// #--------------------------------
		write_list( nfd, fp, list_set );
}



void CNFDB_408::write_rule( CNFData* nfd, FILE* fp, std::vector<crule_rec>& list )
{
	vector<crule_rec>::iterator iter;
	for(iter = list.begin(); iter != list.end(); iter++) {
		nfd->WriteData( fp, "In", iter->rule_type );
		vector<crule_rec::centry_rec>::iterator irec;
		for(irec = iter->entry_set.begin(); irec != iter->entry_set.end(); irec++){
			nfd->WriteData( fp, "IIIIn",
				irec->startID,
				irec->stopID,
				irec->incID,
				irec->include );
		}
		nfd->WriteData( fp, "IIIIn", -1,-1,-1,-1 );
	}
	nfd->WriteData( fp, "In", -1 );
}



void CNFDB_408::write_list( CNFData* nfd, FILE* fp, std::vector<clist_rec>& list )
{
	vector<clist_rec>::iterator iter;
	for(iter = list.begin(); iter != list.end(); iter++) {
		nfd->WriteData( fp, "In", iter->list_type );
		vector<nf_int>::iterator irec;
		for(irec = iter->entityID.begin(); irec != iter->entityID.end(); irec++){
			nfd->WriteData( fp, "In", *irec );
		}
		nfd->WriteData( fp, "In", -1 );
	}
	nfd->WriteData( fp, "In", -1 );
}
