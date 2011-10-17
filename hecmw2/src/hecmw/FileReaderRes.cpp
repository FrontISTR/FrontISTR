/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileReaderRes.cpp
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "HEC_MPI.h"
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
    pmw::CGMGModel *pGMGModel= pmw::CGMGModel::Instance();
    uiint nNumOfLevel, nLevel;
    uiint nNumOfAlgEquation, nAlgEquation;
    uiint nNumOfMesh, nMeshID;
    uiint nNumOfNode, nNodeID;
    uiint nNumOfDOF;
    vdouble vValue;
    istringstream iss;
    if(TagCheck(sLine, FileBlockName::StartRes()) ){
        sLine = getLineSt(ifs);
        iss.clear();
        iss.str(sLine);
        iss >> nNumOfLevel >> nNumOfMesh >> nNumOfAlgEquation;
        for(uiint iLevel=0; iLevel < nNumOfLevel; iLevel++){
            sLine = getLineSt(ifs);
            iss.clear();
            iss.str(sLine);
            iss >> nLevel >> nAlgEquation >> nNumOfDOF;
            pmw::CAssyModel *pAssyModel= pGMGModel->getAssyModel(nLevel);
            pmw::CAssyVector *pSolAssyVec = pAssyModel->getSolutionAssyVector(nAlgEquation);
            for(uiint imesh=0; imesh < nNumOfMesh; imesh++){
                sLine = getLineSt(ifs);
                iss.clear();
                iss.str(sLine);
                iss >> nMeshID >> nNumOfNode;
                pmw::CVector *pSolVec = pSolAssyVec->getVector(nMeshID);
                vValue.resize(nNumOfDOF);
                for(uiint inode=0; inode < nNumOfNode; inode++){
                    sLine = getLineSt(ifs);
                    iss.clear();
                    iss.str(sLine);
                    iss >> nNodeID;
                    for(uiint idof=0; idof < nNumOfDOF; idof++){
                        iss >> vValue[idof];
                        pSolVec->setValue(inode, idof, vValue[idof]);
                    };
                };
            };
        };
        return true;
    }else{
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
    uiint nNumOfLevel, nLevel;
    uiint nNumOfAlgEquation, nAlgEquation;
    uiint nNumOfMesh, nMeshID;
    uiint nNumOfNode, nNodeID;
    uiint nNumOfDOF;
    vdouble vValue;
    ifs.read((char*)&nNumOfLevel, sizeof(uiint));       if(bOrder) pBinCheck->ByteOrderSwap(nNumOfLevel);
    ifs.read((char*)&nNumOfMesh, sizeof(uiint));        if(bOrder) pBinCheck->ByteOrderSwap(nNumOfMesh);
    ifs.read((char*)&nNumOfAlgEquation, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(nNumOfAlgEquation);
    pmw::CGMGModel *pGMGModel= pmw::CGMGModel::Instance();
    for(uiint iLevel=0; iLevel < nNumOfLevel; iLevel++){
        ifs.read((char*)&nLevel, sizeof(uiint));       if(bOrder) pBinCheck->ByteOrderSwap(nLevel);
        ifs.read((char*)&nAlgEquation, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(nAlgEquation);
        ifs.read((char*)&nNumOfDOF, sizeof(uiint));    if(bOrder) pBinCheck->ByteOrderSwap(nNumOfDOF);
        pmw::CAssyModel  *pAssyModel= pGMGModel->getAssyModel(nLevel);
        pmw::CAssyVector *pSolAssyVec= pAssyModel->getSolutionAssyVector(nAlgEquation);
        for(uiint imesh=0; imesh < nNumOfMesh; imesh++){
            ifs.read((char*)&nMeshID, sizeof(uiint));    if(bOrder) pBinCheck->ByteOrderSwap(nMeshID);
            ifs.read((char*)&nNumOfNode, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(nNumOfNode);
            pmw::CVector *pSolVec = pSolAssyVec->getVector(nMeshID);
            vValue.resize(nNumOfDOF);
            for(uiint inode=0; inode < nNumOfNode; inode++){
                ifs.read((char*)&nNodeID, sizeof(uiint));  if(bOrder) pBinCheck->ByteOrderSwap(nNodeID);
                for(uiint idof=0; idof < nNumOfDOF; idof++){
                    ifs.read((char*)&vValue[idof], sizeof(double)); if(bOrder) pBinCheck->ByteOrderSwap(vValue[idof]);
                    pSolVec->setValue(inode, idof, vValue[idof]);
                };
            };
        };
    };
    return true;
}
