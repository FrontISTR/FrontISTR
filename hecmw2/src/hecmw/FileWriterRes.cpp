/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterRes.cpp
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
void CFileWriterRes::WriteAlgebra(ofstream& ofs)
{
    uiint nNumOfLevel= mpGMGModel->getNumOfLevel();

    pmw::CAssyModel *pAssyModel= mpGMGModel->getAssyModel(0);

    ofs << FileBlockName::StartAlgebra() << endl;
    {
        uiint nNumOfEquation= pAssyModel->getNumOfEquation();
        uiint nNumOfMesh= pAssyModel->getNumOfMesh();//---------ローカルMesh数

        // Level数、線形方程式数、Mesh数
        ofs << " " << nNumOfLevel << " " << nNumOfEquation << "  " << nNumOfMesh << endl;

        // 方程式番号、ローカルMesh番号、DOF
        for(uiint ieq=0; ieq < nNumOfEquation; ieq++) {
            pmw::CAssyVector *pSolAssyVec= pAssyModel->getSolutionAssyVector(ieq);
            ofs << " " << ieq;

            // ローカルMeshパーツ別
            for(uiint imesh=0; imesh < nNumOfMesh; imesh++) {
                pmw::CVector *pSolVec= pSolAssyVec->getVector(imesh);

                uiint nDOF = pSolVec->getDOF();
                ofs << " " << imesh << " " << nDOF;
            };//imesh loop
            ofs << endl;

        };//ieq loop
    }
    ofs << FileBlockName::EndAlgebra() << endl;
}
void CFileWriterRes::WriteAlgebra_bin(ofstream& ofs)
{
    uiint nNumOfLevel= mpGMGModel->getNumOfLevel();

    pmw::CAssyModel *pAssyModel= mpGMGModel->getAssyModel(0);

    ofs.write(FileBlockName::StartAlgebra(), FileBlockName::Algebra_Len());

    uiint nNumOfEquation= pAssyModel->getNumOfEquation();
    uiint nNumOfMesh= pAssyModel->getNumOfMesh();//---------ローカルMesh数

    // Level数、線形方程式数、ローカルMesh数
    ofs.write((char*)&nNumOfLevel, sizeof(uiint));
    ofs.write((char*)&nNumOfEquation, sizeof(uiint));
    ofs.write((char*)&nNumOfMesh, sizeof(uiint));

    // 方程式番号、Mesh番号、DOF
    for(uiint ieq=0; ieq < nNumOfEquation; ieq++) {
        pmw::CAssyVector *pSolAssyVec= pAssyModel->getSolutionAssyVector(ieq);
        ofs.write((char*)&ieq, sizeof(uiint));

        // ローカルMeshパーツ別
        for(uiint imesh=0; imesh < nNumOfMesh; imesh++) {
            pmw::CVector *pSolVec= pSolAssyVec->getVector(imesh);
            uiint nDOF= pSolVec->getDOF();

            ofs.write((char*)&imesh, sizeof(uiint));
            ofs.write((char*)&nDOF, sizeof(uiint));
        };
    };
    ofs.write(FileBlockName::End(), FileBlockName::End_Len());
}
void CFileWriterRes::WriteRes(ofstream& ofs)
{
    ofs << FileBlockName::StartRes() << endl;

    uiint nNumOfLevel= mpGMGModel->getNumOfLevel();
    pmw::CAssyModel *pAssyModel= mpGMGModel->getAssyModel(0);

    uiint nNumOfMesh= pAssyModel->getNumOfMesh();
    uiint nNumOfEquation= pAssyModel->getNumOfEquation();

    // MGレベル数, ローカルMesh数, 方程式数
    ofs << " " << nNumOfLevel << " " << nNumOfMesh << " " << nNumOfEquation;

    pmw::CAssyVector *pSolAssyVec;  //コースグリッドの出力もファイングリッドの値を使用するので常にファイングリッド
    pmw::CVector *pSolVec;          //コースグリッドの出力もファイングリッドの値を使用するので常にファイングリッド

    pSolAssyVec= pAssyModel->getSolutionAssyVector(0);

    // ローカルMesh毎: DOF, DOF, DOF, …
    for(uiint imesh=0; imesh < nNumOfMesh; imesh++) {
        pSolVec= pSolAssyVec->getVector(imesh);

        ofs << " " << pSolVec->getDOF();
    };
    ofs << endl;

    for(uiint iLevel=0; iLevel < nNumOfLevel; iLevel++) {
        for(uiint ieq=0; ieq < nNumOfEquation; ieq++) {

            pAssyModel= mpGMGModel->getAssyModel(iLevel);//--- MGレベル指定
            pSolAssyVec= pAssyModel->getSolutionAssyVector(ieq);

            // ローカルMeshパーツ別
            for(uiint imesh=0; imesh < nNumOfMesh; imesh++) {

                pSolVec= pSolAssyVec->getVector(imesh);
                uiint nDOF= pSolVec->getDOF();//Mesh個別のDOF

                pmw::CMesh *pMesh= pAssyModel->getMesh(imesh);
                uiint nNumOfNode= pMesh->getNumOfNode();

                //ノード数
                ofs << nNumOfNode << endl;

                //節点毎: DOF[0]の値、DOF[1]の値、DOF[2]の値、…
                for(uiint inode=0; inode < nNumOfNode; inode++) {
                    for(uiint idof=0; idof < nDOF; idof++)
                        ofs <<  pSolVec->getValue(inode, idof)  << " ";
                    ofs << endl;
                };//inode loop

            };//imesh loop
        };//ieq loop
    };//iLevel loop

    ofs << FileBlockName::EndRes() << endl;
}
void CFileWriterRes::WriteRes_bin(ofstream& ofs)
{
    ofs.write(FileBlockName::StartRes(), FileBlockName::Res_Len());

    uiint nNumOfLevel= mpGMGModel->getNumOfLevel();
    pmw::CAssyModel *pAssyModel= mpGMGModel->getAssyModel(0);

    uiint nNumOfMesh= pAssyModel->getNumOfMesh();
    uiint nNumOfEquation= pAssyModel->getNumOfEquation();

    //MGレベル数, ローカルMesh数, 方程式数
    ofs.write((char*)&nNumOfLevel, sizeof(uiint));
    ofs.write((char*)&nNumOfMesh, sizeof(uiint));
    ofs.write((char*)&nNumOfEquation, sizeof(uiint));

    ////cout << "BinWrite: レベル数:" << nNumOfLevel << " Mesh数:" << nNumOfMesh << " 方程式数:" << nNumOfEquation << endl;

    pmw::CAssyVector *pSolAssyVec;  //コースグリッドの出力もファイングリッドの値を使用するので常にファイングリッド
    pmw::CVector *pSolVec;          //コースグリッドの出力もファイングリッドの値を使用するので常にファイングリッド

    pSolAssyVec= pAssyModel->getSolutionAssyVector(0);

    //ローカルMesh毎: DOF, DOF, DOF, …
    for(uiint imesh=0; imesh < nNumOfMesh; imesh++) {
        pSolVec= pSolAssyVec->getVector(imesh);

        uiint nDOF= pSolVec->getDOF();
        ofs.write((char*)&nDOF, sizeof(uiint));

        ////cout << "BinWrite: imesh:" << imesh << " nDOF:" << nDOF << endl;
    };

    for(uiint iLevel=0; iLevel < nNumOfLevel; iLevel++) {
        for(uiint ieq=0; ieq < nNumOfEquation; ieq++) {

            pAssyModel= mpGMGModel->getAssyModel(iLevel);//--- MGレベル指定
            pSolAssyVec= pAssyModel->getSolutionAssyVector(ieq);

            //ローカルMeshパーツ別
            for(uiint imesh=0; imesh < nNumOfMesh; imesh++) {
                pSolVec= pSolAssyVec->getVector(imesh);
                uiint nDOF= pSolVec->getDOF();

                pmw::CMesh *pMesh= pAssyModel->getMesh(imesh);
                uiint nNumOfNode= pMesh->getNumOfNode();

                //ノード数
                ofs.write((char*)&nNumOfNode, sizeof(uiint));

                ////cout << "BinWrite: imesh:" << imesh << " ノード数:" << nNumOfNode << endl;

                //節点毎: DOF[0]の値、DOF[1]の値、DOF[2]の値、…
                for(uiint inode=0; inode < nNumOfNode; inode++) {
                    ////cout << " inode:" << inode;
                    for(uiint idof=0; idof < nDOF; idof++) {
                        ofs.write((char*)&pSolVec->getValue(inode, idof), sizeof(double));
                        ////cout << " idof:" << idof << " val:" << pSolVec->getValue(inode, idof);
                    };
                    ////cout << endl;
                };
            };
        };
    };
    ofs.write(FileBlockName::EndRes(), FileBlockName::End_Len());
}
