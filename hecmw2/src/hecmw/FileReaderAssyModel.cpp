//
//  FileReaderAssyModel.cpp
//
//              2009.04.08
//              2009.04.08
//              k.Takeda
#include "FileReaderAssyModel.h"
using namespace FileIO;
//using namespace boost;
//  construct & destruct
//
CFileReaderAssyModel::CFileReaderAssyModel()
{
    ;
}

CFileReaderAssyModel::~CFileReaderAssyModel()
{
    ;
}

// アッセンブリーモデル => Meshの領域確保
//
bool CFileReaderAssyModel::Read(ifstream& ifs, string& sLine)
{
    uiint nMeshID, nNumOfMesh, maxID, minID, nProp, mgLevel(0);// mgLevel=0 ::ファイル入力時のマルチグリッド==0
    vuint vMeshID;
    vuint vProp(0);//属性: 0:構造，1:流体
    istringstream iss;

    // MeshIDデータ for AssyModel
    if(TagCheck(sLine, FileBlockName::StartAssyModel()) ){
        
        // メッシュ数,MaxID,MinID
        sLine = getLine(ifs);
        iss.clear();
        iss.str(sLine.c_str());
        iss >> nNumOfMesh >> maxID >> minID;
        
        // setup to BucketMesh in AssyModel
        mpFactory->setupBucketMesh(mgLevel, maxID, minID);
        
        // MeshID の連続データ
        //
        while(!ifs.eof()){
            sLine = getLine(ifs);
            if(TagCheck(sLine, FileBlockName::EndAssyModel()) ) break;

            iss.clear();
            iss.str(sLine.c_str());

            uiint nCount(0);
            while(iss){
                if(nCount==0){ iss >> nMeshID; vMeshID.push_back(nMeshID);}
                if(nCount==1){ iss >> nProp;   vProp.push_back(nProp);}
                nCount++;
                if(nCount > 1) break;
            }
            //cout << "nMeshID " << nMeshID << " nProp " << nProp << endl;
        };
        // Meshの領域確保
        //
        mpFactory->reserveMesh(mgLevel, nNumOfMesh);//ファイル読み込みなので,mgLevel=0

        // Meshの生成 for AssyModel(at mgLevel=0)
        //
        uiint imesh;
        for(imesh=0; imesh < nNumOfMesh; imesh++){
            mpFactory->GeneMesh(mgLevel, vMeshID[imesh], imesh, vProp[imesh]);
        };
        return true;
    }else{
        return false;
    }
}

bool CFileReaderAssyModel::Read_bin(ifstream& ifs)
{
    CFileReaderBinCheck *pBinCheck= CFileReaderBinCheck::Instance();
    bool bOrder= pBinCheck->isByteOrderSwap();
    
    
    bool b32, bCheck;
    string sClassName("FileReaderAssyModel");
    //BinCheckのサイズ指定との整合性
    if( !Check_IntSize(b32, bCheck, sClassName) ) return false;

    //AssyModelタグ
    char cHead='A';
    if( !TagCheck_Bin(ifs, bCheck, cHead, FileBlockName::StartAssyModel(), FileBlockName::AssyModel_Len())) return false;

    
    uiint nMeshID, nNumOfMesh, maxID, minID, nProp, mgLevel(0);// mgLevel=0
    vuint vMeshID;
    vuint vProp(0);//属性: 0:構造，1:流体

    // メッシュ数,MaxID,MinID
    ifs.read((char*)&nNumOfMesh, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(nNumOfMesh);
    ifs.read((char*)&maxID, sizeof(uiint));      if(bOrder) pBinCheck->ByteOrderSwap(maxID);
    ifs.read((char*)&minID, sizeof(uiint));      if(bOrder) pBinCheck->ByteOrderSwap(minID);

    // setup to BucketMesh in AssyModel
    mpFactory->setupBucketMesh(mgLevel, maxID, minID);

    // MeshID の連続データ
    //
    while(!ifs.eof()){
        if(Check_End(ifs)) break;

        ifs.read((char*)&nMeshID, sizeof(uiint));  if(bOrder) pBinCheck->ByteOrderSwap(nMeshID);
        vMeshID.push_back(nMeshID);

        ifs.read((char*)&nProp, sizeof(uiint));    if(bOrder) pBinCheck->ByteOrderSwap(nProp);
        vProp.push_back(nProp);
    };
    // Meshの領域確保
    //
    mpFactory->reserveMesh(mgLevel, nNumOfMesh);//ファイル読み込みなので,mgLevel=0

    // Meshの生成  AssyModel
    //
    uiint imesh;
    for(imesh=0; imesh < nNumOfMesh; imesh++){
        mpFactory->GeneMesh(mgLevel, vMeshID[imesh], imesh, vProp[imesh]);
    };

    return true;
   
}

