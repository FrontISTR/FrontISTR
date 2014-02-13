/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderNode.cpp
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
#include "FileReaderNode.h"
using namespace FileIO;
CFileReaderNode::CFileReaderNode()
{
}
CFileReaderNode::~CFileReaderNode()
{
}
string CFileReaderNode::Name()
{
    return  "FileReaderNode";
}

bool CFileReaderNode::Read(ifstream& ifs, string& sLine)
{
    uiint  nNodeID, nMeshID, nNumOfNode, maxID, minID;
    vdouble vCoord;
    vCoord.resize(3);
    uiint  nType, nNumOfScalarDOF, nNumOfVectorDOF;
    string sType;
    uiint  mgLevel(0);

    if(TagCheck(sLine, FileBlockName::StartNode()) ) {

        //cout << "ReaderNode --- A" << endl;

        sLine = getLineSt(ifs);
        istringstream iss(sLine.c_str());
        iss >> nNumOfNode >> nMeshID >> maxID >> minID;

        //cout << "ReaderNode --- B " << " NumOfNode:" << nNumOfNode
        //        << " nMeshID:" << nMeshID << " maxID:" << maxID << " minID:" << minID << endl;

        mpFactory->reserveNode(mgLevel, nMeshID, nNumOfNode);
        mpFactory->initBucketNode(mgLevel, nMeshID, maxID, minID);

        //cout << "ReaderNode --- C" << endl;

        uiint nCount(0);
        while(!ifs.eof()) {
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndNode()) ) break;

            //cout << "ReaderNode --- D" << endl;

            istringstream iss(sLine.c_str());
            iss  >> sType >> nNumOfScalarDOF >> nNumOfVectorDOF >> nNodeID >> vCoord[0] >> vCoord[1] >> vCoord[2];

            //cout << "ReaderNode --- E" << endl;

            if(sType=="sv"||sType=="SV") {
                nType=pmw::NodeType::ScalarVector;
            } else if(sType=="s"||sType=="S") {
                nType=pmw::NodeType::Scalar;
            } else if(sType=="v"||sType=="V") {
                nType=pmw::NodeType::Vector;
            } else {
                mpLogger->Info(Utility::LoggerMode::Error,"NodeType mismatch, at FileReaderNode");
            }

            if(!mpFactory) mpLogger->Info(Utility::LoggerMode::MWDebug, "Factory => NULL, at FileReaderNode");

            mpFactory->GeneNode(mgLevel, nMeshID, nNodeID, vCoord, nType, nNumOfScalarDOF, nNumOfVectorDOF);
            mpFactory->setIDBucketNode(mgLevel, nMeshID, nNodeID, nCount);
            nCount++;

            //cout << "ReaderNode --- E" << endl;
        };
        mpFactory->setupNode(mgLevel, nMeshID);
        mpFactory->resizeAggregate(mgLevel, nMeshID, nCount);
        mpFactory->GeneAggregate(mgLevel, nMeshID, nCount);

        //cout << "ReaderNode --- F" << endl;

        return true;
    } else {
        return false;
    }
}
bool CFileReaderNode::Read_bin(ifstream& ifs)
{
    return true;
}
