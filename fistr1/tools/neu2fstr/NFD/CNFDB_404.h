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

#ifndef CNFDB_404H
#define CNFDB_404H

#include <vector>
#include "CNFDataBlock.h"


// 404 Element
class CNFDB_404 : public CNFDataBlock {
public:
	CNFDB_404();
	virtual ~CNFDB_404();

	virtual void Read( class CNFData* nfd );
	virtual void WriteData( class CNFData* nfd, FILE* fp );

public:
	// ID of topology
	enum {
		top_Line2 = 0,
		top_Line3,
		top_Tri3,
		top_Tri6,
		top_Quad4,
		top_Quad8,
		top_Tetra4,
		top_Wedge6,
		top_Brick8,
		top_Point,
		top_Tetra10,
		top_Wedge15,
		top_Brick20,
		top_RigidList,
		top_dummy2,
		top_MultiList,
		top_Contact
	};

	class cref_node {
	public:
		nf_int NodeID;
		nf_int faceID;
		nf_float weight;
		nf_int dof[7]; // dof[0]:dummy
		cref_node() {}
		cref_node( const cref_node& r )
		 : NodeID(r.NodeID), faceID(r.faceID),weight(r.weight){
			dof[0] = r.dof[0];
			dof[1] = r.dof[1];
			dof[2] = r.dof[2];
			dof[3] = r.dof[3];
			dof[4] = r.dof[4];
			dof[5] = r.dof[5];
			dof[6] = r.dof[6];
		}
		cref_node& operator=( const cref_node& r ){
			NodeID = r.NodeID;
			faceID = r.faceID;
			weight = r.weight;
			dof[0] = r.dof[0];
			dof[1] = r.dof[1];
			dof[2] = r.dof[2];
			dof[3] = r.dof[3];
			dof[4] = r.dof[4];
			dof[5] = r.dof[5];
			dof[6] = r.dof[6];
			return *this;
		}
	};

	class cref_node_list {
	public:
		std::vector<cref_node> ref_node;
	};

	// #1
		nf_int ID;
		nf_int color;
		nf_int propID;
		nf_int type;
		nf_int topology;
		nf_int layer;
		nf_int orientID;
		nf_bool matl_orflag;
		// Ver.5.0 ========================
		nf_int geomID;
		// Ver.6.0 ========================
		nf_int formulation;
		nf_int contactsegment[2];
		// Ver.8.0 ========================
		nf_int formulation2;
	// #2,3
		nf_int node[20];
	// #4
		nf_float orient[3];
	// #5
		nf_float offset1[3];
	// #6
		nf_float offset2[3];
	// #7
		nf_bool release1[6];
	// #8
		nf_bool release2[6];
		// Ver.4.4 ========================
		nf_int list[4];
	// Ver.4.5 ========================
	// every node record, if topology == top_MultiList or top_RigidList(Rigid?)
		std::vector<cref_node_list*> ref_node_set;

protected:
	cref_node_list* make_ref_node_list( CNFData* nfd );
	void clear_ref_node_set();
	void write_ref_node_list( CNFData* nfd, FILE* fp, cref_node_list* list );
};


#endif
