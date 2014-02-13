/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderElementGroup.cpp
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
#include "FileReaderElementGroup.h"
using namespace FileIO;
CFileReaderElementGroup::CFileReaderElementGroup()
{
    ;
}
CFileReaderElementGroup::~CFileReaderElementGroup()
{
    ;
}
string CFileReaderElementGroup::Name()
{
    return  "FileReaderElementGroup";
}

bool CFileReaderElementGroup::Read(ifstream& ifs, string& sLine)
{
    uiint nGrpID, nMeshID;
    string sGrpName;
    vuint vMeshID;
    map<uiint, vuint, less<uiint> > mavGrpID;
    map<uiint, vstring, less<uiint> > mavGrpName;
    istringstream iss;
    if(TagCheck(sLine, FileBlockName::StartElementGroup()) ) {
        ////cout << "ElementGroup:Block Name : " << sLine << endl;
        while(!ifs.eof()) {
            sLine = getLine(ifs);
            if(TagCheck(sLine, FileBlockName::EndElementGroup()) ) break;
            iss.clear();
            iss.str(sLine.c_str());
            iss >> nGrpID >> sGrpName >> nMeshID;
            ////cout << "GrpID=" << nGrpID << ", GrpName=" << sGrpName << ", MeshID=" << nMeshID << endl;
            vMeshID.push_back(nMeshID);
            mavGrpID[nMeshID].push_back(nGrpID);
            mavGrpName[nMeshID].push_back(sGrpName);
        };
        sort(vMeshID.begin(), vMeshID.end());
        vector<uiint>::iterator unqEnd = unique(vMeshID.begin(), vMeshID.end());
        vMeshID.erase(unqEnd, vMeshID.end());
        vector<uiint>::iterator it;
        for(it=vMeshID.begin(); it != vMeshID.end(); it++) {
            nMeshID = *it;
            mpFactory->GeneElemGrpOBJ(0, nMeshID, mavGrpID[nMeshID], mavGrpName[nMeshID]);
        }
        return true;
    } else {
        return false;
    }
}
bool CFileReaderElementGroup::Read_bin(ifstream& ifs)
{
    return true;
}
