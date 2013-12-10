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

#ifndef CNFDB_506H
#define CNFDB_506H


#include <vector>
#include "CNFDataBlock.h"


// 506 Constraints

class CNFDB_506 : public CNFDataBlock {
public:
	CNFDB_506();
	virtual ~CNFDB_506() {}

	virtual void Read( class CNFData* nfd );
	virtual void WriteData( class CNFData* nfd, FILE* fp );

public:
	class cconst_item {
	public:
		nf_int ID;
		nf_int color;
		nf_int layer;
		nf_bool DOF[6];
		nf_bool ex_geom;
		cconst_item() {}
		cconst_item( const cconst_item& i)
		 : ID(i.ID), color(i.color), layer(i.layer), ex_geom(i.ex_geom) {
			DOF[0] = i.DOF[0];
			DOF[1] = i.DOF[1];
			DOF[2] = i.DOF[2];
			DOF[3] = i.DOF[3];
			DOF[4] = i.DOF[4];
			DOF[5] = i.DOF[5];
		}
		cconst_item& operator=( const cconst_item& i){
			ID = i.ID;
			color = i.color;
			layer = i.layer;
			ex_geom = i.ex_geom;
			DOF[0] = i.DOF[0];
			DOF[1] = i.DOF[1];
			DOF[2] = i.DOF[2];
			DOF[3] = i.DOF[3];
			DOF[4] = i.DOF[4];
			DOF[5] = i.DOF[5];
			return *this;
		}
	};
	class cconst_eq {
	public:
		nf_int ID;
		nf_int color;
		nf_int layer;
	};
	class cnum_co_rec {
	public:
		nf_int eqn_nodeID;
		nf_int eqn_dof;
		nf_float coeff;
	};

	// #1
		nf_int setID;
	// #2
		nf_char title[26];
	// for nodes
		std::vector<cconst_item> const_nodes;
	// for points
		std::vector<cconst_item> const_points;
	// for curves
		std::vector<cconst_item> const_curves;
	// for surfaces
		std::vector<cconst_item> const_surfaces;
	// for number of coefficiency in equation
		std::vector<nf_int> num_co;
	// for equations
		std::vector<cconst_eq> const_equations;
	// for num_co records
		std::vector<cnum_co_rec> num_co_list;
protected:
	void read_const_item( class CNFData* nfd, std::vector<cconst_item>& list );
	void read_const_eq( class CNFData* nfd, std::vector<cconst_eq>& list );
	void read_num_co( class CNFData* nfd, int eq_n, std::vector<nf_int>& list );
	void read_num_co_list( class CNFData* nfd, int eq_n, std::vector<cnum_co_rec>& list );
	void write_const_item( class CNFData* nfd, FILE* fp, std::vector<cconst_item>& list );
	void write_const_eq( class CNFData* nfd, FILE* fp, std::vector<cconst_eq>& list );
	void write_num_co( class CNFData* nfd, FILE* fp, std::vector<nf_int>& list );
	void write_num_co_list( class CNFData* nfd, FILE* fp, std::vector<cnum_co_rec>& list );
};


#endif


