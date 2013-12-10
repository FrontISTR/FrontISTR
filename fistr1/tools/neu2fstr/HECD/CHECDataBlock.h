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
	CHECDataBlock Ver.1.0
*/

#ifndef CHECDataBlockH
#define CHECDataBlockH

#include <stdio.h>


const int hec_name_size = 40;
const int hec_str_size = 256;

class CHECDataBlock {
public:
	int data_type;
	CHECDataBlock(int dtype ) : data_type(dtype){}
	virtual ~CHECDataBlock() {}
	virtual void Clear() = 0;
	virtual void Write( class CHECData* hecd ) = 0;
	virtual bool Read( class CHECData* hecd, char* header_line ) = 0;
	virtual bool IsMesh() { return true; }
};


#endif
