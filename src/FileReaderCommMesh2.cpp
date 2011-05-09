/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderCommMesh2.cxx
|
|                     Written by T.Takeda,    2010/06/01
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
bool CFileReaderCommMesh2::Read(ifstream& ifs, string& sLine)
{
    uint  mgLevel(0);
    uint  nMeshID, numOfCommMesh;
    uint  nCommMeshID, numOfFace, numOfCommNode, myRank, nTransmitRank;
    istringstream iss;
    if(TagCheck(sLine, FileBlockName::StartCommMesh2()) ){
        mpLogger->Info(Utility::LoggerMode::MWDebug, "FileReaderCommMesh2", sLine);
        while(true){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndCommMesh2()) ) break;
            iss.clear();
            iss.str(sLine);
            iss >> nMeshID >> numOfCommMesh;
            for(uint i=0; i< numOfCommMesh; i++){
                sLine = getLineSt(ifs);
                iss.clear();
                iss.str(sLine);
                iss >> nCommMeshID >> numOfFace >> numOfCommNode >> myRank >> nTransmitRank;
                mpFactory->GeneCommMesh2(mgLevel, nMeshID, nCommMeshID,
                                            numOfFace, numOfCommNode,myRank, nTransmitRank);
            };
        };
        return true;
    }else{
        return false;
    }
}
