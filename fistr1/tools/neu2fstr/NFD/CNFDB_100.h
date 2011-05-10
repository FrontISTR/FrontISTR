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
	CNFDB_100 Ver.1.0
	-----------------------------
	100 Header
*/


#ifndef CNFDB_100H
#define CNFDB_100H


#include "CNFDataBlock.h"


// 100 Header
class CNFDB_100 : public CNFDataBlock {
public:
	CNFDB_100();
	virtual ~CNFDB_100() {}

	virtual void Read( class CNFData* nfd );
	virtual void WriteData( class CNFData* nfd, FILE* fp );

public:
	// #1
		nf_char title[256];
	// #2
		nf_float version;
};


#endif
