/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Message.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Message.h"

CMessage::CMessage()
{
}
CMessage::~CMessage()
{
}
//--
//
//--
void CMessage::error(string str)
{
	cout << "Error ---- " << str << endl;
}
void CMessage::warn(string str)
{
	cout << "Warn ----- " << str << endl;
}
void CMessage::info(string str)
{
	cout << "Info ----- " << str << endl;
}
void CMessage::echo(string str)
{
	cout << "Echo ----- " << str << endl;
}
void CMessage::debug(string str)
{
	cout << "Debug ---- " << str << endl;
}
void CMessage::banner()
{
	string sDate;
	struct stat oStat;
	size_t num=26;
	char  *buff= new char[num];
#ifdef MSVC
	if(stat("conv2mw3.exe",&oStat)==0){
		// exe自身のファイル時刻
		errno_t et= ctime_s(buff, num, &oStat.st_mtime);

		sDate = "build date:";
		sDate+= buff;
	}else{
		// 別名で作った場合:現在時刻
		time_t t;
		struct tm timeptr;

		t=time(&t);
		errno_t et= localtime_s(&timeptr,&t);
		asctime_s(buff, num, &timeptr);

		sDate = "      date:";
		sDate+= buff;
	}
	
#elif GNU
	if(stat("conv2mw3",&oStat)==0){
		// exe自身のファイル時刻
		ctime_r(&oStat.st_mtime, buff);

		sDate = "build date:";
		sDate+= buff;
	}else{
		// 別名で作った場合:現在時刻
		time_t t;
		struct tm timeptr;

		t=time(&t);
		localtime_r(&t, &timeptr);
		asctime_r(&timeptr, buff);

		sDate = "      date:";
		sDate+= buff;
	}
#endif
	delete[] buff;

	// bunner display
	cout << "---------------------------------------------------\n"
					"  fistr v4 mesh format >> MW3 mesh format ver " + MW3_Tag::VerNumber() + "  \n"
					"                                                     " << endl;
	cout << "              " << sDate << 
				  "---------------------------------------------------\n" << endl;
}
void CMessage::usage(string str)
{
	cout << "Usage:" << str << endl;
}
//--
//
//--
void CMessage::user_in(string ostr, string& istr)
{
	cout << ostr << flush;
	cin  >> istr;
}


