/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderAssyModel.cpp
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
#include "FileReaderAssyModel.h"
using namespace FileIO;
CFileReaderAssyModel::CFileReaderAssyModel()
{
    ;
}
CFileReaderAssyModel::~CFileReaderAssyModel()
{
    ;
}
string CFileReaderAssyModel::Name()
{
    return "FileReaderAssyModel";
}

bool CFileReaderAssyModel::Read(ifstream& ifs, string& sLine)
{
    uiint nMeshID, nNumOfMesh, maxID, minID, nProp, mgLevel(0);
    vuint vMeshID, vProp(0);

    uiint nNumGlobalComm;
    vector<pair<uiint,uiint> > vRankPair;
    vuint vMeshID_CommID;

    istringstream iss;
    if(TagCheck(sLine, FileBlockName::StartAssyModel()) ) {
        sLine = getLine(ifs);
        iss.clear();
        iss.str(sLine.c_str());
        iss >> nNumOfMesh >> maxID >> minID;
        mpFactory->setupBucketMesh(mgLevel, maxID, minID);

        //アセンブルデータ
        for(uiint imesh=0; imesh < nNumOfMesh; imesh++) {
            sLine = getLine(ifs);
            if(TagCheck(sLine, FileBlockName::EndAssyModel()) ) break;

            iss.clear();
            iss.str(sLine.c_str());
            uiint nCount(0);
            while(iss) {
                if(nCount==0) {
                    iss >> nMeshID;
                    vMeshID.push_back(nMeshID);
                }
                if(nCount==1) {
                    iss >> nProp;      //-----Prop is CodeType
                    vProp.push_back(nProp);
                }
                nCount++;
                if(nCount > 1) break;
            }
        };
        mpFactory->reserveMesh(mgLevel, nNumOfMesh);

        for(uiint imesh=0; imesh < nNumOfMesh; imesh++) {
            mpFactory->GeneMesh(mgLevel, vMeshID[imesh], imesh, vProp[imesh]);
        };

        //グローバル通信テーブル
        sLine= getLine(ifs);
        iss.clear();
        iss.str(sLine.c_str());
        iss >> nNumGlobalComm;

        vRankPair.resize(nNumGlobalComm);
        vMeshID_CommID.resize(nNumGlobalComm);
        //
        // icomm: 暗黙のグローバル通信テーブル番号
        //
        for(uiint icomm=0; icomm < nNumGlobalComm; icomm++) {
            sLine= getLine(ifs);
            iss.clear();
            iss.str(sLine.c_str());
            iss >> vMeshID_CommID[icomm] >> vRankPair[icomm].first >> vRankPair[icomm].second;
        };
        mpFactory->setGlobalCommData(mgLevel, nNumGlobalComm, vRankPair, vMeshID_CommID);

        return true;
    } else {
        return false;
    }
}
bool CFileReaderAssyModel::Read_bin(ifstream& ifs)
{
    CFileReaderBinCheck *pBinCheck= CFileReaderBinCheck::Instance();
    bool bOrder= pBinCheck->isByteOrderSwap();
    bool b32, bCheck;
    string sClassName("FileReaderAssyModel");
    if( !Check_IntSize(b32, bCheck, sClassName) ) return false;
    char cHead='A';
    if( !TagCheck_Bin(ifs, bCheck, cHead, FileBlockName::StartAssyModel(), FileBlockName::AssyModel_Len())) return false;
    uiint nMeshID, nNumOfMesh, maxID, minID, nProp, mgLevel(0);
    vuint vMeshID;
    vuint vProp(0);
    ifs.read((char*)&nNumOfMesh, sizeof(uiint));
    if(bOrder) pBinCheck->ByteOrderSwap(nNumOfMesh);
    ifs.read((char*)&maxID, sizeof(uiint));
    if(bOrder) pBinCheck->ByteOrderSwap(maxID);
    ifs.read((char*)&minID, sizeof(uiint));
    if(bOrder) pBinCheck->ByteOrderSwap(minID);
    mpFactory->setupBucketMesh(mgLevel, maxID, minID);
    while(!ifs.eof()) {
        if(Check_End(ifs)) break;
        ifs.read((char*)&nMeshID, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(nMeshID);
        vMeshID.push_back(nMeshID);
        ifs.read((char*)&nProp, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(nProp);
        vProp.push_back(nProp);
    };
    mpFactory->reserveMesh(mgLevel, nNumOfMesh);
    uiint imesh;
    for(imesh=0; imesh < nNumOfMesh; imesh++) {
        mpFactory->GeneMesh(mgLevel, vMeshID[imesh], imesh, vProp[imesh]);
    };
    return true;
}
