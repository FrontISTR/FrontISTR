//
//  FileReaderElementGroupEntity.cpp
//
//
//              2010.10.22
//              k.Takeda
#include "FileReaderElementGroupEntity.h"
using namespace FileIO;
using namespace boost;

CFileReaderElementGroupEntity::CFileReaderElementGroupEntity()
{
    ;
}
CFileReaderElementGroupEntity::~CFileReaderElementGroupEntity()
{
    ;
}
bool CFileReaderElementGroupEntity::Read(ifstream& ifs, string& sLine)
{
    uint nGrpID;
    uint nElemID, nMeshID;

    vuint vElemID;

    if(TagCheck(sLine, FileBlockName::StartElementGroupEntity()) ){
        //debug
        cout << "ElementGroupEntity:Block Name : " << sLine << endl;

        uint nCount(0);

        while(!ifs.eof()){
            sLine = getLine(ifs);
            if(TagCheck(sLine, FileBlockName::EndElementGroupEntity()) ) break;
            
            if(nCount==0){
                istringstream iss(sLine.c_str());
                iss >> nGrpID >> nMeshID;
            }else{
                // Tokenizer
                char_separator<char> sep(" \t\n");
                tokenizer< char_separator<char> > tokens(sLine, sep);

                typedef tokenizer< char_separator<char> >::iterator Iter;

                for(Iter it=tokens.begin(); it != tokens.end(); ++it){
                    string str = *it;
                    nElemID = atoi(str.c_str());
                    vElemID.push_back(nElemID);
                    //debug
                    cout << " ElemID:" << nElemID;
                };
                //debug
                cout << endl;
            }
            nCount++;
        };
        mpFactory->setElemID_with_ElemGrp(0, nMeshID, nGrpID, vElemID);
        return true;
    }else{
        return false;
    }
}









