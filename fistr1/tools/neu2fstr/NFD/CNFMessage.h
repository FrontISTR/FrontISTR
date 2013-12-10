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

#ifndef CNFMessageH
#define CNFMessageH

#include <stdio.h>
#include <string.h>

// error messages
enum {
	NFE_NO_ERROR = 0,
	NFE_UNKNOWN_ERROR,
	NFE_OPEN_ERROR,
	NFE_READDATA_ERROR,
	NFE_LINE_REQUIRED,
	NFE_DATA_BLOCK_REQUIRED,
	NFE_INVALID_TOKEN,
	NFE_ITEM_REQUIRED,
	NFE_RECORD_REQUIRED,

	NFE_WRITEDATA_ERROR
};


// warning messages
enum {
	NFW_NON_SUPPORTED_DATA_BLOCK
};


class CNFMessage {
protected:
	static char msg[256];
public:
	int no;
	int line;
	int column;
	char option_msg[256];

	CNFMessage(int No = 0, const char* opt_msg="", int Line=-1, int Col=-1 )
	 : no(No), line(Line), column(Col) {
		strcpy( option_msg, opt_msg );
	}
	virtual ~CNFMessage() {}
	virtual const char* Msg() = 0;
};


class CNFError : public CNFMessage {
public:
	CNFError(int No, const char* opt_msg="", int Line=-1, int Col=-1 )
	 : CNFMessage(No, opt_msg, Line, Col ) {}
	CNFError(int No, int Line, int Col=-1 )
	 : CNFMessage(No, "", Line, Col ) {}
	virtual ~CNFError() {}
	virtual const char* Msg();
};


class CNFWarning : public CNFMessage {
public:
	CNFWarning(int No, const char* opt_msg="", int Line=-1, int Col=-1 )
	 : CNFMessage(No, opt_msg, Line, Col ) {}
	CNFWarning(int No, int Line, int Col=-1 )
	 : CNFMessage(No, "", Line, Col ) {}
	virtual ~CNFWarning() {}
	virtual const char* Msg();
};


#endif

