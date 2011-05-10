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
	CNFDB_506 Ver.1.0
*/

// 506 Constraints


#include "CNFData.h"
#include "CNFDB_506.h"

using namespace std;



CNFDB_506::CNFDB_506()
 : CNFDataBlock(506)
{}

void CNFDB_506::Read( CNFData* nfd )
{
	char buff[256];

	// #1
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "I", &setID );
	// #2
		nfd->ReadLineEx( buff );
		nfd->ReadStr( buff, title, sizeof(title));
	// for nodes
		read_const_item( nfd, const_nodes );
	// for points
		read_const_item( nfd, const_points );
	// for curves
		read_const_item( nfd, const_curves );
	// for surfaces
		read_const_item( nfd, const_surfaces );
	// for equations
		read_const_eq( nfd, const_equations );

	int eq_n = const_equations.size();
	if(eq_n==0) return;

	// for number of coefficiency in equation
		read_num_co( nfd, eq_n, num_co );
	// for num_co records
		read_num_co_list( nfd, eq_n, num_co_list );
}


void CNFDB_506::read_const_item( CNFData* nfd, std::vector<cconst_item>& list )
{
	char buff[256];
	cconst_item item;
	list.clear();
	while(1) {
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIIBBBBBBB", 
			&item.ID,
			&item.color,
			&item.layer,
			&item.DOF[0], &item.DOF[1], &item.DOF[2],
			&item.DOF[3], &item.DOF[4], &item.DOF[5],
			&item.ex_geom
		);
		if( item.ID == -1 ) break;
		list.push_back( item );
	}
}



void CNFDB_506::read_const_eq( CNFData* nfd, std::vector<cconst_eq>& list )
{
	char buff[256];
	cconst_eq item;
	list.clear();
	while(1) {
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "III",
			&item.ID,
			&item.color,
			&item.layer
		);
		if( item.ID == -1 ) break;
		list.push_back( item );
	}
}


void CNFDB_506::read_num_co( CNFData* nfd, int eq_n, std::vector<nf_int>& list )
{
	char buff[256];
	nf_int x;
	list.clear();
	for(int i=0; i<eq_n; i++) {
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "I", &x );
		list.push_back( x );
	}
}


void CNFDB_506::read_num_co_list( CNFData* nfd, int eq_n, std::vector<cnum_co_rec>& list )
{
	char buff[256];
	cnum_co_rec item;
	list.clear();
	for(int i=0; i<eq_n; i++) {
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIF",
			&item.eqn_nodeID,
			&item.eqn_dof,
			&item.coeff
		);
		list.push_back( item );
	}
}


//*****************************************************************************


void CNFDB_506::WriteData( class CNFData* nfd, FILE* fp )
{
	// #1
		nfd->WriteData( fp, "In", setID );
	// #2
		nfd->WriteStr( fp, title );
	// for nodes
		write_const_item( nfd, fp, const_nodes );
	// for points
		write_const_item( nfd, fp, const_points );
	// for curves
		write_const_item( nfd, fp, const_curves );
	// for surfaces
		write_const_item( nfd, fp, const_surfaces );
	// for equations
		write_const_eq( nfd, fp, const_equations );

	int eq_n = const_equations.size();
	if(eq_n==0) return;

	// for number of coefficiency in equation
		write_num_co( nfd, fp, num_co );
	// for num_co records
		write_num_co_list( nfd, fp, num_co_list );
}


void CNFDB_506::write_const_item( CNFData* nfd, FILE* fp, std::vector<CNFDB_506::cconst_item>& list )
{
	vector<cconst_item>::iterator iter;
	for(iter = list.begin(); iter != list.end(); iter++ ){
		nfd->WriteData( fp, "IIIBBBBBBBn", 
			iter->ID,
			iter->color,
			iter->layer,
			iter->DOF[0], iter->DOF[1], iter->DOF[2],
			iter->DOF[3], iter->DOF[4], iter->DOF[5],
			iter->ex_geom
		);
	}
	nfd->WriteData( fp, "IIIBBBBBBBn", -1,-1,-1,0,0,0,0,0,0,0 );
}



void CNFDB_506::write_const_eq( CNFData* nfd, FILE* fp, std::vector<CNFDB_506::cconst_eq>& list )
{
	vector<cconst_eq>::iterator iter;
	for(iter = list.begin(); iter != list.end(); iter++ ){
		nfd->WriteData( fp, "IIIn",
			iter->ID,
			iter->color,
			iter->layer
		);
	}
	nfd->WriteData( fp, "IIIn", -1,-1,-1 );
}


void CNFDB_506::write_num_co( CNFData* nfd, FILE* fp, std::vector<nf_int>& list )
{
	vector<nf_int>::iterator iter;
	for(iter = list.begin(); iter != list.end(); iter++ ){
		nfd->WriteData( fp, "In", *iter );
	}
}


void CNFDB_506::write_num_co_list( CNFData* nfd, FILE* fp, std::vector<CNFDB_506::cnum_co_rec>& list )
{
	vector<cnum_co_rec>::iterator iter;
	for(iter = list.begin(); iter != list.end(); iter++ ){
		nfd->WriteData( fp, "IIFn",
			iter->eqn_nodeID,
			iter->eqn_dof,
			iter->coeff
		);
	}
}


