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

/* conv_neu2fstr_dynamic ver.1.0 */

#ifndef conv_neu2fstr_dynamicH
#define conv_neu2fstr_dynamicH

#include "CNFData.h"
#include "CHECData.h"
#include "CHECDB.h"
#include "CFSTRDB.h"

bool conv_neu2fstr_dynamic( CNFData& neu, CHECData& hec );


#endif
