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
	CNFMessage ver.1.0
*/


#include <stdio.h>
#include "CNFMessage.h"


static char ERROR_MSG[][80] = {
	"No error",
	"Unknown error",
	"Cannot open NEU file",
	"Read data error",
	"Data line required",
	"Data block required",
	"Invalid token",
	"Item requred",
	"A record is required"
};


static char WARNING_MSG[][80] = {
	"Non supported data block"
};



char CNFMessage::msg[256] = "";


const char* CNFError::Msg()
{
	char s[256];

	strcpy( msg, "##Error" );
	if( line >= 0 ){
		sprintf(s,"(line:%d", line );
		strcat( msg, s );
		if( column >0 ){
			sprintf(s,",col:%d", column );
			strcat( msg, s );
		}
		strcat( msg, ")");
	}
	strcat( msg, ": ");
	strcat( msg, ERROR_MSG[no] );
	strcat( msg, option_msg );
	return msg;
}



const char* CNFWarning::Msg()
{
	char s[256];

	strcpy( msg, "##Warning" );
	if( line >= 0 ){
		sprintf(s,"(line:%d", line );
		strcat( msg, s );
		if( column >0 ){
			sprintf(s,",col:%d", column );
			strcat( msg, s );
		}
		strcat( msg, ")");
	}
	strcat( msg, ": ");
	strcat( msg, WARNING_MSG[no] );
	strcat( msg, option_msg );
	return msg;
}
