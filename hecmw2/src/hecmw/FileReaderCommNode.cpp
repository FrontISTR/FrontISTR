/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderCommNode.cpp
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
#include "FileReaderCommNode.h"
using namespace FileIO;
CFileReaderCommNode::CFileReaderCommNode()
{
    ;
}
CFileReaderCommNode::~CFileReaderCommNode()
{
    ;
}
string CFileReaderCommNode::Name()
{
    return  "FileReaderCommNode";
}

bool CFileReaderCommNode::Read(ifstream& ifs, string& sLine)
{
    uiint numOfNode, nCommMeshID, maxID, minID;
    uiint nCommNodeID, nMeshID, nNodeID, nRank;
    uiint  mgLevel(0);
    if(TagCheck(sLine, FileBlockName::StartCommNode()) ) {
        sLine = getLineSt(ifs);
        istringstream iss(sLine.c_str());
        iss >> numOfNode >> nMeshID >> nCommMeshID >> maxID >> minID;
        mpFactory->reserveCommNode(mgLevel,nMeshID, nCommMeshID, numOfNode);
        while(true) {
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndCommNode()) ) break;
            istringstream iss(sLine.c_str());
            iss >> nCommNodeID >> nNodeID >> nRank;
            mpFactory->GeneCommNode(mgLevel, nCommNodeID, nMeshID, nCommMeshID, nNodeID, nRank);
        };
        return true;
    } else {
        return false;
    }
}
bool CFileReaderCommNode::Read_bin(ifstream& ifs)
{
    return true;
}
