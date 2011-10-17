/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileWriterRes.cpp
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
    pmw::CAssyModel *pAssyModel = mpGMGModel->getAssyModel(0);
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
void CFileWriterRes::WriteAlgebra_bin(ofstream& ofs)
{
    pmw::CAssyModel *pAssyModel = mpGMGModel->getAssyModel(0);
    ofs.write(FileBlockName::StartAlgebra(), FileBlockName::Algebra_Len());
    uiint nNumOfEquation = pAssyModel->getNumOfEquation();
    ofs.write((char*)&nNumOfEquation, sizeof(uiint));
    for(uiint ieq=0; ieq < nNumOfEquation; ieq++){
        pmw::CAssyVector *pSolAssyVec = pAssyModel->getSolutionAssyVector(ieq);
        uiint nDOF= pSolAssyVec->getDOF();
        ofs.write((char*)&ieq, sizeof(uiint));
        ofs.write((char*)&nDOF, sizeof(uiint));
    }
    ofs.write(FileBlockName::End(), FileBlockName::End_Len());
}
void CFileWriterRes::WriteRes(ofstream& ofs)
{
    ofs << FileBlockName::StartRes() << endl;
    uiint nNumOfLevel = mpGMGModel->getNumOfLevel();
    pmw::CAssyModel *pAssyModel= mpGMGModel->getAssyModel(0);
    uiint nNumOfMesh= pAssyModel->getNumOfMesh();
    uiint nNumOfEquation= pAssyModel->getNumOfEquation();
    ofs << " " << nNumOfLevel << " " << nNumOfMesh << " " << nNumOfEquation << endl;
    pmw::CAssyVector *pFGSolAssyVec;
    pmw::CVector *pFGSolVec;        
    pmw::CAssyVector *pSolAssyVec;  
    pmw::CVector *pSolVec;          
    for(uiint iLevel=0; iLevel < nNumOfLevel; iLevel++){
        for(uiint ieq=0; ieq < nNumOfEquation; ieq++){
            if(iLevel != nNumOfLevel-1){
            pAssyModel= mpGMGModel->getAssyModel(nNumOfLevel-1);  
            pFGSolAssyVec= pAssyModel->getSolutionAssyVector(ieq);
            }
            pAssyModel= mpGMGModel->getAssyModel(iLevel);
            pSolAssyVec= pAssyModel->getSolutionAssyVector(ieq);
            uiint nNumOfDOF = pSolAssyVec->getDOF();
            ofs << " " << iLevel << " " << ieq << " " << nNumOfDOF << endl;
            for(uiint imesh=0; imesh < nNumOfMesh; imesh++){
                if(iLevel != nNumOfLevel-1){
                pFGSolVec= pFGSolAssyVec->getVector(imesh);
                }
                pSolVec= pSolAssyVec->getVector(imesh);
                pmw::CMesh *pMesh= pAssyModel->getMesh(imesh);
                uiint nMeshID= pMesh->getMeshID();
                uiint nNumOfNode= pMesh->getNumOfNode();
                ofs << " " << nMeshID << " " << nNumOfNode << endl;
                for(uiint inode=0; inode < nNumOfNode; inode++){
                    pmw::CNode *pNode = pMesh->getNodeIX(inode);
                    uiint nNodeID = pNode->getID();
                    ofs << " " << nNodeID << " ";
                    for(uiint idof=0; idof < nNumOfDOF; idof++){
                        if(iLevel != nNumOfLevel-1){
                            ofs << pFGSolVec->getValue(inode, idof) << " ";
                        }else{
                            ofs <<  pSolVec->getValue(inode, idof)  << " ";
                        }
                    };
                    ofs << endl;
                };
            };
        };
    };
    ofs << FileBlockName::EndRes() << endl;
}
void CFileWriterRes::WriteRes_bin(ofstream& ofs)
{
    ofs.write(FileBlockName::StartRes(), FileBlockName::Res_Len());
    uiint nNumOfLevel = mpGMGModel->getNumOfLevel();
    pmw::CAssyModel *pAssyModel= mpGMGModel->getAssyModel(0);
    uiint nNumOfMesh= pAssyModel->getNumOfMesh();
    uiint nNumOfEquation= pAssyModel->getNumOfEquation();
    ofs.write((char*)&nNumOfLevel, sizeof(uiint));
    ofs.write((char*)&nNumOfMesh, sizeof(uiint));
    ofs.write((char*)&nNumOfEquation, sizeof(uiint));
    pmw::CAssyVector *pFGSolAssyVec;
    pmw::CVector *pFGSolVec;        
    pmw::CAssyVector *pSolAssyVec;  
    pmw::CVector *pSolVec;          
    for(uiint iLevel=0; iLevel < nNumOfLevel; iLevel++){
        for(uiint ieq=0; ieq < nNumOfEquation; ieq++){
            if(iLevel != nNumOfLevel-1){
            pAssyModel= mpGMGModel->getAssyModel(nNumOfLevel-1);  
            pFGSolAssyVec= pAssyModel->getSolutionAssyVector(ieq);
            }
            pAssyModel= mpGMGModel->getAssyModel(iLevel);
            pSolAssyVec= pAssyModel->getSolutionAssyVector(ieq);
            uiint nNumOfDOF = pSolAssyVec->getDOF();
            ofs.write((char*)&iLevel, sizeof(uiint));
            ofs.write((char*)&ieq, sizeof(uiint));
            ofs.write((char*)&nNumOfDOF, sizeof(uiint));
            for(uiint imesh=0; imesh < nNumOfMesh; imesh++){
                if(iLevel != nNumOfLevel-1){
                pFGSolVec= pFGSolAssyVec->getVector(imesh);
                }
                pSolVec= pSolAssyVec->getVector(imesh);
                pmw::CMesh *pMesh= pAssyModel->getMesh(imesh);
                uiint nMeshID= pMesh->getMeshID();
                uiint nNumOfNode= pMesh->getNumOfNode();
                ofs.write((char*)&nMeshID, sizeof(uiint));
                ofs.write((char*)&nNumOfNode, sizeof(uiint));
                for(uiint inode=0; inode < nNumOfNode; inode++){
                    pmw::CNode *pNode = pMesh->getNodeIX(inode);
                    uiint nNodeID = pNode->getID();
                    ofs.write((char*)&nNodeID, sizeof(uiint));
                    for(uiint idof=0; idof < nNumOfDOF; idof++){
                        if(iLevel != nNumOfLevel-1){
                            ofs.write((char*)&pFGSolVec->getValue(inode, idof), sizeof(double));
                        }else{
                            ofs.write((char*)&pSolVec->getValue(inode, idof), sizeof(double));
                        }
                    };
                };
            };
        };
    };
    ofs.write(FileBlockName::EndRes(), FileBlockName::End_Len());
}
