/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderVersion.cpp
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FileReaderVersion.h"
using namespace FileIO;

CFileReaderVersion::CFileReaderVersion()
{
    msVersion="00";
}
CFileReaderVersion::~CFileReaderVersion()
{
}
string CFileReaderVersion::Name()
{
    return "FileReaderVersion";
}
string& CFileReaderVersion::getVersion()
{
    return msVersion;
}

bool CFileReaderVersion::Read(ifstream& ifs, string& sLine)
{
    if( TagCheck(sLine, FileBlockName::StartVersion()) ) {
        while(!ifs.eof()) {
            sLine = getLineSt(ifs);
            if( TagCheck(sLine, FileBlockName::EndNode()) ) return false;
            //--
            // バージョン(ファイル書式バージョン): バージョン0.4 = "MW3_0.4"
            //--
            string::size_type pos= sLine.find("MW3_");
            if( pos != string::npos ) {
                pos= sLine.find("_");
                string str2= sLine.substr(pos+1);//---- '_'以降の取得(バージョン番号 \n)

                str2= Cleaning(str2);//---- 特殊記号の排除
                msVersion= str2;

                return true;
            } else {
                msVersion="-0";//--------------"MW3_"プリフィックスすら無い、Versionブロックだけ存在.
                return true;
            }
        };
    } else {
        return false;
    }
}
bool CFileReaderVersion::Read_bin(ifstream& ifs)
{
    return true;
}


