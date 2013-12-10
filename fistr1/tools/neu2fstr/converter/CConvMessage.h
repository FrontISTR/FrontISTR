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
/* CConvMessage class Ver.1.0 */

#ifndef CConvMessageH
#define CConvMessageH

#include <stdarg.h>
#include <stdio.h>
#include <string.h>


enum {
	CONV_NO_ERROR = 0,
	CONV_UNKNOWN_ERROR,
	CONV_COORDINATE_ERROR,
	CONV_NO_SUPPORTED_ELEMENT,
	CONV_INVALID_ELEMENT_PROPERTY,
	CONV_NO_SUPPORTED_ELEM_PROP
};


class CConvMessage {
protected:
	static char msg[256];
public:
	char option_msg[256];
	int no;
	CConvMessage(int No = 0, const char* op_msg="", ... );
	virtual const char* Msg();
};


#endif
