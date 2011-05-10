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


#include "CConvMessage.h"


const char ERROR_MSG[][80] = {
	"No error",
	"Unknown error",
	"Coordinate error",
	"Not supported element",
	"Invalid element property",
	"Not supported property of element"
};




char CConvMessage::msg[256] = "";



CConvMessage::CConvMessage(int No, const char* op_msg, ... )
 : no(No)
{
	if( op_msg[0] == 0 ) {
		option_msg[0] = 0;
		return;
	}

	va_list va;
	va_start(va, op_msg );
	vsprintf( option_msg, op_msg, va );
	va_end( va );
}


const char* CConvMessage::Msg()
{
	if( option_msg[0] != 0 ){
		sprintf( msg, "##Error: %s : %s",  ERROR_MSG[no], option_msg );
	} else {
		sprintf( msg, "##Error: %s", ERROR_MSG[no]);
	}
	return msg;
}
