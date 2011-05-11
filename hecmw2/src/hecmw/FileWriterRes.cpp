//
// FileWriterRes.cpp
//
//              2011.03.04
//              k.Takeda
#include "FileWriterRes.h"
#include "FileBlockName.h"
#include "AssyVector.h"
using namespace FileIO;


CFileWriterRes::CFileWriterRes()
{
    ;
}
CFileWriterRes::~CFileWriterRes()
{
    ;
}

void CFileWriterRes::WriteDebug(ofstream& ofs, const uiint& mgLevel)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();

    pLogger->Info(Utility::LoggerMode::Warn, "invalid method, FileWriterRes::Write");
}

//
// Algebraブロック出力(各線形方程式のDOF)
//
void CFileWriterRes::WriteAlgebra(ofstream& ofs)
{
    pmw::CAssyModel *pAssyModel = mpGMGModel->getAssyModel(0);//コース・グリッド

    ofs << FileBlockName::StartAlgebra() << endl;
    {
        uiint nNumOfEquation = pAssyModel->getNumOfEquation();
        ofs << " " << nNumOfEquation << endl;

        for(uiint ieq=0; ieq < nNumOfEquation; ieq++){
            pmw::CAssyVector *pSolAssyVec = pAssyModel->getSolutionAssyVector(ieq);
            uiint nDOF = pSolAssyVec->getDOF();

            ofs << " " << ieq << " " << nDOF << endl;
        };
    }
    ofs << FileBlockName::EndAlgebra() << endl;
}

// Binary
void CFileWriterRes::WriteAlgebra_bin(ofstream& ofs)
{
    pmw::CAssyModel *pAssyModel = mpGMGModel->getAssyModel(0);//コース・グリッド

    ofs.write(FileBlockName::StartAlgebra(), FileBlockName::Algebra_Len());//Block Name

    uiint nNumOfEquation = pAssyModel->getNumOfEquation();
    ofs.write((char*)&nNumOfEquation, sizeof(uiint));

    for(uiint ieq=0; ieq < nNumOfEquation; ieq++){
        pmw::CAssyVector *pSolAssyVec = pAssyModel->getSolutionAssyVector(ieq);
        uiint nDOF= pSolAssyVec->getDOF();

        ofs.write((char*)&ieq, sizeof(uiint));
        ofs.write((char*)&nDOF, sizeof(uiint));
    }
    ofs.write(FileBlockName::End(), FileBlockName::End_Len());//End Block
}


//
// Resブロック出力(SolutionVectorの値)
//
void CFileWriterRes::WriteRes(ofstream& ofs)
{
    // Resブロック
    ofs << FileBlockName::StartRes() << endl;
    
    uiint nNumOfLevel = mpGMGModel->getNumOfLevel();
    pmw::CAssyModel *pAssyModel= mpGMGModel->getAssyModel(0);
    uiint nNumOfMesh= pAssyModel->getNumOfMesh();
    uiint nNumOfEquation= pAssyModel->getNumOfEquation();
    
    ofs << " " << nNumOfLevel << " " << nNumOfMesh << " " << nNumOfEquation << endl;// MG数, パーツ数, 方程式数
    
    pmw::CAssyVector *pFGSolAssyVec;//解アセンブル・ベクトル(ファイン・グリッド)
    pmw::CVector *pFGSolVec;        //解ベクトル(ファイン・グリッド)
    pmw::CAssyVector *pSolAssyVec;  //解アセンブル・ベクトル(カレント・グリッド)
    pmw::CVector *pSolVec;          //解ベクトル(カレント・グリッド)
    
    for(uiint iLevel=0; iLevel < nNumOfLevel; iLevel++){
        for(uiint ieq=0; ieq < nNumOfEquation; ieq++){

            if(iLevel != nNumOfLevel-1){
            pAssyModel= mpGMGModel->getAssyModel(nNumOfLevel-1);  //AssyModel(ファイン・グリッド)
            pFGSolAssyVec= pAssyModel->getSolutionAssyVector(ieq);//解アセンブル・ベクトル(ファイン・グリッド)
            }
            
            pAssyModel= mpGMGModel->getAssyModel(iLevel);//AssyModel(カレント・グリッド)

            pSolAssyVec= pAssyModel->getSolutionAssyVector(ieq);//解アセンブル・ベクトル(カレント・グリッド)
            uiint nNumOfDOF = pSolAssyVec->getDOF();
            
            ofs << " " << iLevel << " " << ieq << " " << nNumOfDOF << endl;// MG番号, 方程式番号, 自由度数

            for(uiint imesh=0; imesh < nNumOfMesh; imesh++){
                if(iLevel != nNumOfLevel-1){
                pFGSolVec= pFGSolAssyVec->getVector(imesh);//解ベクトル(ファイン・グリッド)
                }
                pSolVec= pSolAssyVec->getVector(imesh);//解ベクトル(カレント・グリッド)
                
                pmw::CMesh *pMesh= pAssyModel->getMesh(imesh);
                uiint nMeshID= pMesh->getMeshID();
                uiint nNumOfNode= pMesh->getNumOfNode();
                
                ofs << " " << nMeshID << " " << nNumOfNode << endl;// パーツ番号, ノード数

                for(uiint inode=0; inode < nNumOfNode; inode++){
                    pmw::CNode *pNode = pMesh->getNodeIX(inode);
                    uiint nNodeID = pNode->getID();
                    
                    ofs << " " << nNodeID << " ";

                    // コース・グリッドの解ベクトル出力に,ファイン・グリッドの結果を当てはめる
                    //
                    for(uiint idof=0; idof < nNumOfDOF; idof++){
                        if(iLevel != nNumOfLevel-1){
                            ofs << pFGSolVec->getValue(inode, idof) << " ";
                        }else{
                            ofs <<  pSolVec->getValue(inode, idof)  << " ";
                        }
                    };
                    ofs << endl;
                    
                };// Node loop end
            };// Mesh loop end
        };// Algebra_Equation loop end
    };// MultiGrid loop end
    
    // ブロックエンド
    ofs << FileBlockName::EndRes() << endl;
}

// Binary
void CFileWriterRes::WriteRes_bin(ofstream& ofs)
{
    // Resブロック
    ofs.write(FileBlockName::StartRes(), FileBlockName::Res_Len());
    
    uiint nNumOfLevel = mpGMGModel->getNumOfLevel();
    pmw::CAssyModel *pAssyModel= mpGMGModel->getAssyModel(0);
    uiint nNumOfMesh= pAssyModel->getNumOfMesh();
    uiint nNumOfEquation= pAssyModel->getNumOfEquation();
    
    // MG数, パーツ数, 方程式数
    ofs.write((char*)&nNumOfLevel, sizeof(uiint));
    ofs.write((char*)&nNumOfMesh, sizeof(uiint));
    ofs.write((char*)&nNumOfEquation, sizeof(uiint));
    
    pmw::CAssyVector *pFGSolAssyVec;//解アセンブル・ベクトル(ファイン・グリッド)
    pmw::CVector *pFGSolVec;        //解ベクトル(ファイン・グリッド)
    pmw::CAssyVector *pSolAssyVec;  //解アセンブル・ベクトル(カレント・グリッド)
    pmw::CVector *pSolVec;          //解ベクトル(カレント・グリッド)
    
    for(uiint iLevel=0; iLevel < nNumOfLevel; iLevel++){
        for(uiint ieq=0; ieq < nNumOfEquation; ieq++){

            if(iLevel != nNumOfLevel-1){
            pAssyModel= mpGMGModel->getAssyModel(nNumOfLevel-1);  //AssyModel(ファイン・グリッド)
            pFGSolAssyVec= pAssyModel->getSolutionAssyVector(ieq);//解アセンブル・ベクトル(ファイン・グリッド)
            }
            
            pAssyModel= mpGMGModel->getAssyModel(iLevel);//AssyModel(カレント・グリッド)

            pSolAssyVec= pAssyModel->getSolutionAssyVector(ieq);//解アセンブル・ベクトル(カレント・グリッド)
            uiint nNumOfDOF = pSolAssyVec->getDOF();

            // MG番号, 方程式番号, 自由度数
            ofs.write((char*)&iLevel, sizeof(uiint));
            ofs.write((char*)&ieq, sizeof(uiint));
            ofs.write((char*)&nNumOfDOF, sizeof(uiint));

            for(uiint imesh=0; imesh < nNumOfMesh; imesh++){
                if(iLevel != nNumOfLevel-1){
                pFGSolVec= pFGSolAssyVec->getVector(imesh);//解ベクトル(ファイン・グリッド)
                }
                pSolVec= pSolAssyVec->getVector(imesh);//解ベクトル(カレント・グリッド)
                
                pmw::CMesh *pMesh= pAssyModel->getMesh(imesh);
                uiint nMeshID= pMesh->getMeshID();
                uiint nNumOfNode= pMesh->getNumOfNode();
                
                // パーツ番号, ノード数
                ofs.write((char*)&nMeshID, sizeof(uiint));
                ofs.write((char*)&nNumOfNode, sizeof(uiint));

                for(uiint inode=0; inode < nNumOfNode; inode++){
                    pmw::CNode *pNode = pMesh->getNodeIX(inode);
                    uiint nNodeID = pNode->getID();
                    
                    ofs.write((char*)&nNodeID, sizeof(uiint));

                    // コース・グリッドの解ベクトル出力に,ファイン・グリッドの結果を当てはめる
                    //
                    for(uiint idof=0; idof < nNumOfDOF; idof++){
                        if(iLevel != nNumOfLevel-1){
                            ofs.write((char*)&pFGSolVec->getValue(inode, idof), sizeof(double));
                        }else{
                            ofs.write((char*)&pSolVec->getValue(inode, idof), sizeof(double));
                        }
                    };
                    
                };// Node loop end
            };// Mesh loop end
        };// Algebra_Equation loop end
    };// MultiGrid loop end
    
    // ブロックエンド
    ofs.write(FileBlockName::EndRes(), FileBlockName::End_Len());
}




