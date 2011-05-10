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
	CNFDB_405 Ver.1.0
*/

#ifndef CNFDB_405H
#define CNFDB_405H


#include "CNFDataBlock.h"


// 405 Coordnate Systems
class CNFDB_405 : public CNFDataBlock {
public:
	CNFDB_405();
	virtual ~CNFDB_405() {}

	virtual void Read( class CNFData* nfd );
	virtual void WriteData( class CNFData* nfd, FILE* fp );

public:
	// #1
		nf_int ID;
		nf_int define_sys;
		nf_int type;
		nf_int color;
		nf_int layer;
	// #2
		nf_char title[26];
	// #3
		nf_float origin[3];
	// #4
		nf_float rot[3];
};


#endif
