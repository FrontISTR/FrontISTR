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
	conv_neu2hec ver.1.0
*/


#ifndef conv_neu2hecH
#define conv_neu2hecH


#include "CNFData.h"
#include "CNFDataBlock.h"
#include "CHECData.h"
#include "CHECDB.h"

// solution
enum {
	sol_static = 0,
	sol_heat,
	sol_eigen
};

void conv_neu2hec( CNFData& neu, CHECData& hec, int solution=0 );



#endif
