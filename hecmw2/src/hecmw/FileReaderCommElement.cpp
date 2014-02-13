/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderCommElement.cpp
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
#include "FileReaderCommElement.h"
using namespace FileIO;
CFileReaderCommElement::CFileReaderCommElement()
{
    ;
}
CFileReaderCommElement::~CFileReaderCommElement()
{
    //cout << "~CFileReaderCommElement" << endl;
}
string CFileReaderCommElement::Name()
{
    return  "FileReaderCommElement";
}

bool CFileReaderCommElement::Read(ifstream& ifs, string& sLine)
{
    uiint mgLevel(0);
    uiint  numOfCommElement, nMeshID, nCommMeshID, nMaxCommID, nMinCommID;
    uiint  nElementID;
    vuint vCommNodeID;
    string sElemType;
    uiint   nElemType;
    if(TagCheck(sLine, FileBlockName::StartCommElement()) ) {

        //cout << "FileReaderCommElement -- A" << endl;

        sLine = getLineSt(ifs);
        istringstream iss(sLine.c_str());
        iss >> numOfCommElement >> nMeshID >> nCommMeshID >> nMaxCommID >> nMinCommID;

        mpFactory->reserveCommElement(mgLevel, nMeshID, nCommMeshID, numOfCommElement);

        //cout << "FileReaderCommElement -- B" << endl;

        while(true) {
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndCommElement()) ) break;
            istringstream iss(sLine.c_str());
            iss >> sElemType >> nElementID;
            nElemType= IntElemType(sElemType);
            vCommNodeID.clear();
            uiint ivert;
            switch(nElemType) {
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
            //cout << "FileReaderCommElement -- C" << endl;

            mpFactory->GeneCommElement(mgLevel, nMeshID, nCommMeshID, nElemType, nElementID, vCommNodeID);

            //cout << "FileReaderCommElement -- D" << endl;
        };
        return true;
    } else {
        return false;
    }
}
bool CFileReaderCommElement::Read_bin(ifstream& ifs)
{
    mpLogger->Info(Utility::LoggerMode::Warn, "FileReaderCommElement Read_bin invalid method");
    return false;
}
