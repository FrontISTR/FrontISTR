/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderCommMesh2.cpp
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
#include "FileReaderCommMesh2.h"
using namespace FileIO;
CFileReaderCommMesh2::CFileReaderCommMesh2()
{
    ;
}
CFileReaderCommMesh2::~CFileReaderCommMesh2()
{
    ;
}
string CFileReaderCommMesh2::Name()
{
    return  "FileReaderCommMesh2";
}

bool CFileReaderCommMesh2::Read(ifstream& ifs, string& sLine)
{
    uiint  mgLevel(0);
    uiint  nMeshID, numOfCommMesh;
    uiint  nCommMeshID, numOfFace, numOfCommNode, myRank, nTransmitRank;
    istringstream iss;
    if(TagCheck(sLine, FileBlockName::StartCommMesh2()) ) {
        mpLogger->Info(Utility::LoggerMode::MWDebug, "FileReaderCommMesh2", sLine);
        while(true) {
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndCommMesh2()) ) break;
            iss.clear();
            iss.str(sLine);
            iss >> nMeshID >> numOfCommMesh;
            for(uiint i=0; i< numOfCommMesh; i++) {
                sLine = getLineSt(ifs);
                iss.clear();
                iss.str(sLine);
                iss >> nCommMeshID >> numOfFace >> numOfCommNode >> myRank >> nTransmitRank;
                mpFactory->GeneCommMesh2(mgLevel, nMeshID, nCommMeshID,
                                         numOfFace, numOfCommNode, myRank, nTransmitRank);
            };
        };
        return true;
    } else {
        return false;
    }
}
bool CFileReaderCommMesh2::Read_bin(ifstream& ifs)
{
    return true;
}
