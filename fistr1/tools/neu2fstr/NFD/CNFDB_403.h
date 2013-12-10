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
	CNFDB_403 Ver.1.0
*/

#ifndef CNFDB_403H
#define CNFDB_403H


#include "CNFDataBlock.h"


// 403 Node
class CNFDB_403 : public CNFDataBlock {
public:
	CNFDB_403();
	virtual ~CNFDB_403() {}
	virtual void Read( class CNFData* nfd );
	virtual void WriteData( class CNFData* nfd, FILE* fp );

public:
	// #1
		nf_int ID;
		nf_int define_sys;
		nf_int output_sys;
		nf_int layer;
		nf_int color;
		nf_bool permbc[6];
		nf_float x;
		nf_float y;
		nf_float z;
	// ======= Ver.4.4 ========================
		nf_int node_type;
};


#endif
