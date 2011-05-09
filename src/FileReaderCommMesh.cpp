/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderCommMesh.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FileReaderCommMesh.h"
using namespace FileIO;
CFileReaderCommMesh::CFileReaderCommMesh()
{
    ;
}
CFileReaderCommMesh::~CFileReaderCommMesh()
{
    ;
}
bool CFileReaderCommMesh::Read(ifstream& ifs, string& sLine)
{
    uint mgLevel(0);
    uint  numOfMesh, maxMeshID, minMeshID;
    uint  nMeshID, numOfCommMesh;
    uint  nCommMeshID, myRank,nTransmitRank;
    string   sLinePara;
    vstring svLinePara; svLinePara.resize(3);
    string   white(" ");
    if(TagCheck(sLine, FileBlockName::StartCommMesh()) ){
        mpLogger->Info(Utility::LoggerMode::MWDebug, "FileReaderCommMesh", sLine);
        sLine = getLineSt(ifs);
        istringstream iss(sLine.c_str());
        iss >> numOfMesh >> maxMeshID >> minMeshID;
        while(true){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndCommMesh()) ) break;
            istringstream iss(sLine.c_str());
            iss >> nMeshID >> numOfCommMesh;
            mpFactory->reserveCommMesh(mgLevel, nMeshID, numOfCommMesh);
            for(uint i=0; i< numOfCommMesh; i++){
                sLine = getLineSt(ifs);
                istringstream issin(sLine.c_str());
                issin >> nCommMeshID >> myRank >> nTransmitRank;
                svLinePara[0]= boost::lexical_cast<string>(nCommMeshID);
                svLinePara[1]= boost::lexical_cast<string>(myRank);
                svLinePara[2]= boost::lexical_cast<string>(nTransmitRank);
                sLinePara = svLinePara[0] + white + svLinePara[1] + white + svLinePara[2];
                mpLogger->Info(Utility::LoggerMode::MWDebug, "CommMesh, myRank, TransmitRank => ", sLinePara);
                mpFactory->GeneCommMesh(mgLevel, nMeshID, nCommMeshID, myRank, nTransmitRank);
            };
        };
        return true;
    }else{
        return false;
    }
}
