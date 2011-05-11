//
//  FileReaderRes.cpp
//
//
//              2011.02.25
//              k.Takeda
#include "FileReaderRes.h"
#include "AssyVector.h"
using namespace FileIO;

CFileReaderRes::CFileReaderRes()
{
    ;
}
CFileReaderRes::~CFileReaderRes()
{
    ;
}


bool CFileReaderRes::Read(ifstream& ifs, string& sLine)
{
    //全体 階層モデル
    pmw::CGMGModel *pGMGModel= pmw::CGMGModel::Instance();

    uiint nNumOfLevel, nLevel;
    uiint nNumOfAlgEquation, nAlgEquation;
    uiint nNumOfMesh, nMeshID;
    uiint nNumOfNode, nNodeID;
    uiint nNumOfDOF;
    vdouble vValue;
    
    istringstream iss;

    // -----------------------------------------------------
    // Res ブロック (リスタート・データ) : 線形方程式の解ベクトル値
    // -----------------------------------------------------
    if(TagCheck(sLine, FileBlockName::StartRes()) ){
        sLine = getLineSt(ifs);
        iss.clear();
        iss.str(sLine);
        
        iss >> nNumOfLevel >> nNumOfMesh >> nNumOfAlgEquation;// マルチグリッド・レベル数, メッシュ・パーツ数, 線形方程式数
        
        for(uiint iLevel=0; iLevel < nNumOfLevel; iLevel++){
            sLine = getLineSt(ifs);
            iss.clear();
            iss.str(sLine);
            
            iss >> nLevel >> nAlgEquation >> nNumOfDOF;//マルチグリッド番号, 方程式番号, 自由度数
            
            //アセンブル・モデル => 線形方程式取得
            pmw::CAssyModel *pAssyModel= pGMGModel->getAssyModel(nLevel);
            pmw::CAssyVector *pSolAssyVec = pAssyModel->getSolutionAssyVector(nAlgEquation);//アセンブル_解ベクトル(nAlgEquation番目)

            // 各メッシュ・パーツ 読み込み
            for(uiint imesh=0; imesh < nNumOfMesh; imesh++){
                
                sLine = getLineSt(ifs);
                iss.clear();
                iss.str(sLine);
                
                iss >> nMeshID >> nNumOfNode;//メッシュ・パーツID, ノード数
                
                pmw::CVector *pSolVec = pSolAssyVec->getVector(nMeshID);//パーツ別_解ベクトル

                vValue.resize(nNumOfDOF);

                for(uiint inode=0; inode < nNumOfNode; inode++){
                    sLine = getLineSt(ifs);
                    iss.clear();
                    iss.str(sLine);
                    // ノードID,  ベクトル値(自由度数分)
                    iss >> nNodeID;
                    for(uiint idof=0; idof < nNumOfDOF; idof++){
                        iss >> vValue[idof];

                        //線形方程式の解ベクトルに代入
                        pSolVec->setValue(inode, idof, vValue[idof]);
                    };// DOF loop end
                };// Node loop end
            };// Mesh loop end
        };// Level loop end
        
        return true;
    }else{
        return false;
    }
}

bool CFileReaderRes::Read_bin(ifstream& ifs)
{
    CFileReaderBinCheck *pBinCheck= CFileReaderBinCheck::Instance();
    bool bOrder= pBinCheck->isByteOrderSwap();
    
    //BinCheckのサイズ指定との整合性
    bool b32, bCheck;
    string sClassName("FileReaderRes");

    if( !Check_IntSize(b32, bCheck, sClassName) ) return false;

    char cHead='R';
    if( !TagCheck_Bin(ifs, bCheck, cHead, FileBlockName::StartRes(), FileBlockName::Res_Len())) return false;

    
    uiint nNumOfLevel, nLevel;
    uiint nNumOfAlgEquation, nAlgEquation;
    uiint nNumOfMesh, nMeshID;
    uiint nNumOfNode, nNodeID;
    uiint nNumOfDOF;
    vdouble vValue;
    
    // -----------------------------------------------------
    // Res ブロック (リスタート・データ) : 線形方程式の解ベクトル値
    // -----------------------------------------------------
    
    // マルチグリッド・レベル数, メッシュ・パーツ数, 線形方程式数
    ifs.read((char*)&nNumOfLevel, sizeof(uiint));       if(bOrder) pBinCheck->ByteOrderSwap(nNumOfLevel);
    ifs.read((char*)&nNumOfMesh, sizeof(uiint));        if(bOrder) pBinCheck->ByteOrderSwap(nNumOfMesh);
    ifs.read((char*)&nNumOfAlgEquation, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(nNumOfAlgEquation);
    

    //全体 階層モデル
    pmw::CGMGModel *pGMGModel= pmw::CGMGModel::Instance();
    
    for(uiint iLevel=0; iLevel < nNumOfLevel; iLevel++){
        //マルチグリッド番号, 方程式番号, 自由度数
        ifs.read((char*)&nLevel, sizeof(uiint));       if(bOrder) pBinCheck->ByteOrderSwap(nLevel);
        ifs.read((char*)&nAlgEquation, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(nAlgEquation);
        ifs.read((char*)&nNumOfDOF, sizeof(uiint));    if(bOrder) pBinCheck->ByteOrderSwap(nNumOfDOF);


        //アセンブル・モデル => 線形方程式取得
        pmw::CAssyModel  *pAssyModel= pGMGModel->getAssyModel(nLevel);
        pmw::CAssyVector *pSolAssyVec= pAssyModel->getSolutionAssyVector(nAlgEquation);//アセンブル_解ベクトル(nAlgEquation番目)

        // 各メッシュ・パーツ 読み込み
        for(uiint imesh=0; imesh < nNumOfMesh; imesh++){

            //メッシュ・パーツID, ノード数
            ifs.read((char*)&nMeshID, sizeof(uiint));    if(bOrder) pBinCheck->ByteOrderSwap(nMeshID);
            ifs.read((char*)&nNumOfNode, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(nNumOfNode);

            pmw::CVector *pSolVec = pSolAssyVec->getVector(nMeshID);//パーツ別_解ベクトル

            vValue.resize(nNumOfDOF);

            for(uiint inode=0; inode < nNumOfNode; inode++){
                // ノードID,  ベクトル値(自由度数分)
                ifs.read((char*)&nNodeID, sizeof(uiint));  if(bOrder) pBinCheck->ByteOrderSwap(nNodeID);

                for(uiint idof=0; idof < nNumOfDOF; idof++){
                    ifs.read((char*)&vValue[idof], sizeof(double)); if(bOrder) pBinCheck->ByteOrderSwap(vValue[idof]);

                    //線形方程式の解ベクトルに代入
                    pSolVec->setValue(inode, idof, vValue[idof]);

                };// DOF loop end
            };// Node loop end
        };// Mesh loop end
    };// Level loop end
    
    //    uiint nLength;
    //    nLength= FileBlockName::End_Len();
    //    char cTag[nLength+1];
    //    ifs.read(cTag, FileBlockName::End_Len());//End

    return true;
}







