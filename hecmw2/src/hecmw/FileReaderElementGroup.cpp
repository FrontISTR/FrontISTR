//
//  FileReaderElementGroup.cpp
//
//
//          2010.10.22
//          k.Takeda
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

bool CFileReaderElementGroup::Read(ifstream& ifs, string& sLine)
{
    uint nGrpID, nMeshID;
    string sGrpName;
    
    vuint vMeshID;
    map<uint, vuint, less<uint> > mavGrpID;
    map<uint, vstring, less<uint> > mavGrpName;

    istringstream iss;

    if(TagCheck(sLine, FileBlockName::StartElementGroup()) ){
        //debug
        cout << "ElementGroup:Block Name : " << sLine << endl;

        while(!ifs.eof()){

            sLine = getLine(ifs);

            if(TagCheck(sLine, FileBlockName::EndElementGroup()) ) break;
            
            iss.clear();
            iss.str(sLine.c_str());

            iss >> nGrpID >> sGrpName >> nMeshID;

            //debug
            cout << "GrpID=" << nGrpID << ", GrpName=" << sGrpName << ", MeshID=" << nMeshID << endl;

            vMeshID.push_back(nMeshID);
            mavGrpID[nMeshID].push_back(nGrpID);
            mavGrpName[nMeshID].push_back(sGrpName);
        };

        sort(vMeshID.begin(), vMeshID.end());
        vector<uint>::iterator unqEnd = unique(vMeshID.begin(), vMeshID.end());
        vMeshID.erase(unqEnd, vMeshID.end());

        vector<uint>::iterator it;
        for(it=vMeshID.begin(); it != vMeshID.end(); it++){
            nMeshID = *it;
            mpFactory->GeneElemGrpOBJ(0, nMeshID, mavGrpID[nMeshID], mavGrpName[nMeshID]);
        }

        return true;
    }else{
        return false;
    }
}


















