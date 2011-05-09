/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderCnt.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FileReaderCnt.h"
#include "FileBlockName.h"
using namespace FileIO;
CFileReaderCnt::CFileReaderCnt()
{
    ;
}
CFileReaderCnt::~CFileReaderCnt()
{
    ;
}
bool CFileReaderCnt::Read(ifstream& ifs, string& sLine)
{
    if(TagCheck(sLine, FileBlockName::StartMeshFileName()) ){
        while(true){
            sLine = getLineSt(ifs);
            if(sLine==FileBlockName::EndMeshFileName()) break;
            istringstream iss(sLine.c_str());
            iss >> msMeshFileBaseName;
            mpLogger->Info(Utility::LoggerMode::Debug,"ベースファイル名=> ",msMeshFileBaseName);
        };
        return true;
    }else{
        return false;
    }
}
