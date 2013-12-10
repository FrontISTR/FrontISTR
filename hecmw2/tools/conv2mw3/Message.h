/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Message.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/

#ifndef C5865088_771A_46a1_9911_9D5F41B3734C
#define C5865088_771A_46a1_9911_9D5F41B3734C

#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h> //struct stat
#include <ctime>
using namespace std;

#include "FrontISTR_Type.h"

class CMessage{
private:
	CMessage();
public:
	static CMessage* Instance(){
		static CMessage oMssg;
		return &oMssg;
	}
	virtual ~CMessage();

public:
	void error(string str);
	void warn(string str);
	void info(string str);
	void echo(string str);
	void debug(string str);
	void banner();

	void usage(string str);

	void user_in(string ostr, string& istr);
};
#endif //include_guard




