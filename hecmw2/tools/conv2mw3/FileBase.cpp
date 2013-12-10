/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   FileBase.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FileBase.h"

CFileBase::CFileBase()
{
}
CFileBase::~CFileBase()
{
}
vstring CFileBase::SplitToken(string sLine)
{
	vstring vStr;

	char_separator<char> sep(", \t\r\n");
	tokenizer< char_separator<char> > tokens(sLine, sep);

	typedef tokenizer< char_separator<char> >::iterator Iter;

	for(Iter it=tokens.begin(); it != tokens.end(); ++it){
		vStr.push_back(*it);
	};

	return vStr;
}
void CFileBase::replaceComma(string& sLine)
{
	size_t nLength=sLine.length();

	for(size_t i=0; i < nLength; i++){
		if( ','==sLine[i]) sLine[i]=' ';
	};
}
bool CFileBase::existStr(const string& sLine, string str)
{
	string::size_type pos = sLine.find(str);

	if(string::npos != pos){
		return true;
	}else{
		return false;
	}
}
string CFileBase::getRearStr(string sLine, string symbol)
{
	string sRear;

	string::size_type pos= sLine.find(symbol);
	if( pos!=string::npos ){
		sRear = sLine.substr(pos+1);
	}else{
		sRear = "";//Error
	}

	return sRear;
}



