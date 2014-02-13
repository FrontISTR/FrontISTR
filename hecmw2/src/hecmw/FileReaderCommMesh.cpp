/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderCommMesh.cpp
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
string CFileReaderCommMesh::Name()
{
    return  "FileReaderCommMesh";
}


bool CFileReaderCommMesh::Read(ifstream& ifs, string& sLine)
{
    uiint mgLevel(0);
    uiint  numOfMesh, maxMeshID, minMeshID;
    uiint  nMeshID, numOfCommMesh;
    uiint  nCommMeshID, myRank,nTransmitRank;
    string   sLinePara;
    vstring svLinePara;
    svLinePara.resize(3);
    string   white(" ");
    if(TagCheck(sLine, FileBlockName::StartCommMesh()) ) {
        mpLogger->Info(Utility::LoggerMode::MWDebug, "FileReaderCommMesh", sLine);
        sLine = getLineSt(ifs);
        istringstream iss(sLine.c_str());
        iss >> numOfMesh >> maxMeshID >> minMeshID;
        while(true) {
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndCommMesh()) ) break;
            istringstream iss(sLine.c_str());
            iss >> nMeshID >> numOfCommMesh;
            mpFactory->reserveCommMesh(mgLevel, nMeshID, numOfCommMesh);

            for(uiint i=0; i< numOfCommMesh; i++) {
                sLine = getLineSt(ifs);
                istringstream issin(sLine.c_str());
                issin >> nCommMeshID >> myRank >> nTransmitRank;

                //////    svLinePara[0]= boost::lexical_cast<string>(nCommMeshID);
                //////    svLinePara[1]= boost::lexical_cast<string>(myRank);
                //////    svLinePara[2]= boost::lexical_cast<string>(nTransmitRank);
                //////    sLinePara = svLinePara[0] + white + svLinePara[1] + white + svLinePara[2];
                //////    mpLogger->Info(Utility::LoggerMode::MWDebug, "CommMesh, myRank, TransmitRank => ", sLinePara);

                mpFactory->GeneCommMesh(mgLevel, nMeshID, nCommMeshID, myRank, nTransmitRank);
            };
        };
        return true;
    } else {
        return false;
    }
}
bool CFileReaderCommMesh::Read_bin(ifstream& ifs)
{
    return true;
}
