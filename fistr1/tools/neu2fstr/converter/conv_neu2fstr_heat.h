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
	conv_neu2fstr_heat ver.1.0
*/


#ifndef conv_neu2fstr_heatH
#define conv_neu2fstr_heatH

#include "CNFData.h"
#include "CHECData.h"
#include "CHECDB.h"
#include "CFSTRDB.h"

void conv_neu2fstr_heat( CNFData& neu, CHECData& hec );


#endif
