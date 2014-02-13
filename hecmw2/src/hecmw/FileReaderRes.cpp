/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderRes.cpp
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
string CFileReaderRes::Name()
{
    return  "FileReaderRes";
}

bool CFileReaderRes::Read(ifstream& ifs, string& sLine)
{
    pmw::CGMGModel *pGMGModel= pmw::CGMGModel::Instance();

    uiint nNumOfLevel;
    uiint nNumOfAlgEquation, nNumOfMesh, nNumOfNode;
    vvdouble vvValue;

    istringstream iss;
    if(TagCheck(sLine, FileBlockName::StartRes()) ) {
        sLine = getLineSt(ifs);
        iss.clear();
        iss.str(sLine);
        //MGレベル数, ローカルMesh数, 方程式数
        iss >> nNumOfLevel >> nNumOfMesh >> nNumOfAlgEquation;

        //ローカルMesh毎: DOF, DOF, DOF, …
        vuint vDOF;
        vDOF.resize(nNumOfMesh);
        for(uiint imesh=0; imesh < nNumOfMesh; imesh++) {
            iss >> vDOF[imesh];
        };

        for(uiint iLevel=0; iLevel < nNumOfLevel; iLevel++) {
            pmw::CAssyModel *pAssyModel= pGMGModel->getAssyModel(iLevel);//iLevelのアセンブル・モデル

            for(uiint ieq=0; ieq < nNumOfAlgEquation; ieq++) {
                pmw::CAssyVector *pSolAssyVec= pAssyModel->getSolutionAssyVector(ieq);//ieq番の方程式の解ベクトル

                // Meshパーツ別 (ローカルMesh)
                for(uiint imesh=0; imesh < nNumOfMesh; imesh++) {
                    pmw::CVector *pSolVec= pSolAssyVec->getVector(imesh);//imesh番のパーツ解ベクトル

                    sLine = getLineSt(ifs);
                    iss.clear();
                    iss.str(sLine);

                    //ノード数
                    iss >> nNumOfNode;

                    vvValue.resize(nNumOfNode);
                    for(uiint inode=0; inode < nNumOfNode; inode++)
                        vvValue[inode].resize(vDOF[imesh]);

                    //節点毎: DOF[0]の値、DOF[1]の値、DOF[2]の値、…
                    for(uiint inode=0; inode < nNumOfNode; inode++) {
                        sLine = getLineSt(ifs);
                        iss.clear();
                        iss.str(sLine);
                        for(uiint idof=0; idof < vDOF[imesh]; idof++) {
                            iss >> vvValue[inode][idof];
                            pSolVec->setValue(inode, idof, vvValue[inode][idof]);
                        };
                    };//inode loop
                };//imesh loop
            };//ieq loop
        };//iLevel loop

        return true;
    } else {
        return false;
    }
}
bool CFileReaderRes::Read_bin(ifstream& ifs)
{
    CFileReaderBinCheck *pBinCheck= CFileReaderBinCheck::Instance();

    bool bOrder= pBinCheck->isByteOrderSwap();
    bool b32, bCheck;
    string sClassName("FileReaderRes");
    if( !Check_IntSize(b32, bCheck, sClassName) ) return false;
    char cHead='R';
    if( !TagCheck_Bin(ifs, bCheck, cHead, FileBlockName::StartRes(), FileBlockName::Res_Len())) return false;

    uiint nNumOfLevel;
    uiint nNumOfAlgEquation, nNumOfMesh, nNumOfNode;
    vvdouble vvValue;

    //MGレベル数, ローカルMesh数, 方程式数
    ifs.read((char*)&nNumOfLevel, sizeof(uiint));
    if(bOrder) pBinCheck->ByteOrderSwap(nNumOfLevel);
    ifs.read((char*)&nNumOfMesh, sizeof(uiint));
    if(bOrder) pBinCheck->ByteOrderSwap(nNumOfMesh);
    ifs.read((char*)&nNumOfAlgEquation, sizeof(uiint));
    if(bOrder) pBinCheck->ByteOrderSwap(nNumOfAlgEquation);

    //ローカルMesh毎: DOF, DOF, DOF, …
    vuint vDOF;
    vDOF.resize(nNumOfMesh);
    for(uiint imesh=0; imesh < nNumOfMesh; imesh++) {
        ifs.read((char*)&vDOF[imesh], sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(vDOF[imesh]);
    };

    pmw::CGMGModel *pGMGModel= pmw::CGMGModel::Instance();

    for(uiint iLevel=0; iLevel < nNumOfLevel; iLevel++) {
        pmw::CAssyModel  *pAssyModel= pGMGModel->getAssyModel(iLevel);

        for(uiint ieq=0; ieq < nNumOfAlgEquation; ieq++) {
            pmw::CAssyVector *pSolAssyVec= pAssyModel->getSolutionAssyVector(ieq);

            // Meshパーツ別 (ローカルMesh)
            for(uiint imesh=0; imesh < nNumOfMesh; imesh++) {
                pmw::CVector *pSolVec= pSolAssyVec->getVector(imesh);//imesh番のパーツ解ベクトル

                //ノード数
                ifs.read((char*)&nNumOfNode, sizeof(uiint));
                if(bOrder) pBinCheck->ByteOrderSwap(nNumOfNode);

                vvValue.resize(nNumOfNode);
                for(uiint inode=0; inode < nNumOfNode; inode++)
                    vvValue[inode].resize(vDOF[imesh]);

                //節点毎: DOF[0]の値、DOF[1]の値、DOF[2]の値、…
                for(uiint inode=0; inode < nNumOfNode; inode++) {
                    ////cout << " inode:" << inode;
                    for(uiint idof=0; idof < vDOF[imesh]; idof++) {

                        ifs.read((char*)&vvValue[inode][idof], sizeof(double));
                        if(bOrder) pBinCheck->ByteOrderSwap(vvValue[inode][idof]);

                        pSolVec->setValue(inode, idof, vvValue[inode][idof]);

                        ////cout << " idof:" << idof << " val:" << vvValue[inode][idof];
                    };//idof loop
                    ////cout << endl;
                };//inode loop

            };//imesh loop
        };//ieq loop
    };//iLevel loop
    return true;
}
