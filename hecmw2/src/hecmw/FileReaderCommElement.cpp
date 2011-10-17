/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileReaderCommElement.cpp
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "HEC_MPI.h"
#include "FileReaderCommElement.h"
using namespace FileIO;
CFileReaderCommElement::CFileReaderCommElement()
{
    ;
}
CFileReaderCommElement::~CFileReaderCommElement()
{
    ;
}
bool CFileReaderCommElement::Read(ifstream& ifs, string& sLine)
{
    uiint mgLevel(0);
    uiint  numOfCommElement, nMeshID, nCommMeshID, nMaxCommID, nMinCommID;
    uiint  nElementID;
    vuint vCommNodeID;
    string sElemType;
    uiint   nElemType;
    if(TagCheck(sLine, FileBlockName::StartCommElement()) ){
        sLine = getLineSt(ifs);
        istringstream iss(sLine.c_str());
        iss >> numOfCommElement >> nMeshID >> nCommMeshID >> nMaxCommID >> nMinCommID;
        mpFactory->reserveCommElement(mgLevel, nMeshID, nCommMeshID, numOfCommElement);
        while(true){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndCommElement()) ) break;
            istringstream iss(sLine.c_str());
            iss >> sElemType >> nElementID;
            nElemType= IntElemType(sElemType);
            vCommNodeID.clear();
            uiint ivert;
            switch(nElemType){
                case(pmw::ElementType::Hexa):
                    vCommNodeID.resize(8);
                    for(ivert=0; ivert< 8; ivert++) iss >> vCommNodeID[ivert];
                    break;
                case(pmw::ElementType::Tetra):
                    vCommNodeID.resize(4);
                    for(ivert=0; ivert< 4; ivert++) iss >> vCommNodeID[ivert];
                    break;
                case(pmw::ElementType::Prism):
                    vCommNodeID.resize(6);
                    for(ivert=0; ivert< 6; ivert++) iss >> vCommNodeID[ivert];
                    break;
                case(pmw::ElementType::Quad):
                    vCommNodeID.resize(4);
                    for(ivert=0; ivert< 4; ivert++) iss >> vCommNodeID[ivert];
                    break;
                case(pmw::ElementType::Triangle):
                    vCommNodeID.resize(3);
                    for(ivert=0; ivert< 3; ivert++) iss >> vCommNodeID[ivert];
                    break;
                case(pmw::ElementType::Beam):
                    vCommNodeID.resize(2);
                    for(ivert=0; ivert< 2; ivert++) iss >> vCommNodeID[ivert];
                    break;
                default:
                    break;
            }
            mpFactory->GeneCommElement(mgLevel, nMeshID, nCommMeshID, nElemType, nElementID, vCommNodeID);
        };
        return true;
    }else{
        return false;
    }
}
bool CFileReaderCommElement::Read_bin(ifstream& ifs)
{
    mpLogger->Info(Utility::LoggerMode::Warn, "FileReaderCommElement Read_bin invalid method");
    return false;
}
