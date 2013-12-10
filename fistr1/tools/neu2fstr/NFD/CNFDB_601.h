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
	CNFDB_601 Ver.1.0
	-----------------------------
	601 Material
*/

#ifndef CNFDB_601H
#define CNFDB_601H


#include "CNFDataBlock.h"


// 601 Material
class CNFDB_601 : public CNFDataBlock {
public:
	CNFDB_601();
	virtual ~CNFDB_601() {}
	
	virtual void Read( class CNFData* nfd );
	virtual void WriteData( class CNFData* nfd, FILE* fp );

public:
	class cfunc_rec {
	public:
		class cdata_pair {
		public:
			nf_int index;
			nf_float x;
			nf_float y;
			cdata_pair(nf_int I=0, nf_float X=0, nf_float Y=0)
			 : index(I), x(X), y(Y) {}
		};
		// ##1
		nf_int ID;
		nf_int type;
		// ##2
		nf_char title[80];
		// ## ----------------------
		std::vector<cdata_pair> data;
	};

	// #1
		nf_int ID;
		nf_int format;// == -601
		nf_int color;
		nf_int type;
		nf_int subtype;
		nf_int layer;
		nf_int FunctionCount;
	// #2
		nf_char title[26];
	// #3
		nf_int Bcount; // == 10
	// #4
		nf_bool bval[10];
	// #5
		nf_int Icount;
	// #6,7,8
		nf_int ival[25];
	// #9
		nf_int Mcount;
	// #10-29
		nf_float mval[200];
	// #30
		nf_int Fcount;
	// #31-35
		nf_int fval[50];
	// #36
		nf_int Tcount;
	// #37-43
		nf_int tval[70];
	// # ---------------------
		std::vector<cfunc_rec> func_list;

	// ============================================================================
	nf_float& E(int i) { assert(0<=i&&i<=2); return mval[i];   }
	nf_float& G(int i) { assert(0<=i&&i<=2); return mval[i+3]; }
	nf_float& NU(int i){ assert(0<=i&&i<=2); return mval[i+6]; }
	nf_float& THERMAL_EXPANSION(int i)   { assert(0<=i&&i<=5); return mval[i+36]; }
	nf_float& THERMAL_CONDUCTIVITY(int i){ assert(0<=i&&i<=5); return mval[i+42]; }
	nf_float& THERMAL_CAPACITY()         { return mval[48]; }
	nf_float& DENSITY()	{ return mval[49]; }
	nf_float& DAMPING()	{ return mval[50]; }
	nf_float& TEMPERATURE()	{ return mval[51]; }
};


#endif
