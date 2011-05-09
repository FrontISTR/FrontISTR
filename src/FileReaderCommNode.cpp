/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderCommNode.cxx
|
|                     Written by T.Takeda,    2010/06/01
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
bool CFileReaderCommNode::Read(ifstream& ifs, string& sLine)
{
    uint numOfNode, nCommMeshID, maxID, minID;
    uint nCommNodeID, nMeshID, nNodeID, nRank;
    uint  mgLevel(0); 
    if(TagCheck(sLine, FileBlockName::StartCommNode()) ){
        sLine = getLineSt(ifs);
        istringstream iss(sLine.c_str());
        iss >> numOfNode >> nMeshID >> nCommMeshID >> maxID >> minID;
        mpFactory->reserveCommNode(mgLevel,nMeshID, nCommMeshID, numOfNode);
        while(true){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndCommNode()) ) break;
            istringstream iss(sLine.c_str());
            iss >> nCommNodeID >> nNodeID >> nRank;
            mpFactory->GeneCommNode(mgLevel, nCommNodeID, nMeshID, nCommMeshID, nNodeID, nRank);
        };
        return true;
    }else{
        return false;
    }
}
