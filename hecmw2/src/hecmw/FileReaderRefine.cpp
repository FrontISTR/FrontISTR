/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderRefine.cpp
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
#include "FileReaderRefine.h"
using namespace FileIO;
CFileReaderRefine::CFileReaderRefine()
{
    ;
}
CFileReaderRefine::~CFileReaderRefine()
{
    ;
}
string CFileReaderRefine::Name()
{
    return  "FileReaderRefine";
}

bool CFileReaderRefine::Read(ifstream& ifs, string& sLine)
{
    uiint numOfRefine;
    if(TagCheck(sLine, FileBlockName::StartRefine()) ) {
        while(true) {
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndRefine()) ) break;
            istringstream iss(sLine.c_str());
            iss >> numOfRefine;
            mpLogger->Info(Utility::LoggerMode::MWDebug, "FileReaderRefine::numOfRefine == ", numOfRefine);
        };
        mpFactory->GeneAssyModel(numOfRefine+1);
        mpFactory->setMGLevel(numOfRefine);
        return true;
    } else {
        return false;
    }
}
bool CFileReaderRefine::Read_bin(ifstream& ifs)
{
    return true;
}
