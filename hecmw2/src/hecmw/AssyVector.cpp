/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/AssyVector.cpp
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
#include "HEC_MPI.h"
#include "AssyVector.h"
#include "AssyModel.h"
#include "CommNode.h"
#include "CommMesh2.h"
namespace pmw
{
typedef std::vector<CVector*> CVVec;
typedef CVVec::iterator CVVecIter;
typedef CVVec::const_iterator CVVecConstIter;


CAssyVector::CAssyVector()
{
    ;
}
//
// 引数 vDOF:グローバルMesh番号(MeshID)に対応するDOF
//
void CAssyVector::initAssyVector(CAssyModel* pAssyModel, vuint& vDOF)
{
    ////CHecMPI *pMPI= CHecMPI::Instance();
    ////cout << "-- AssyVector::initAssyVector ------ enter -- rank:" << pMPI->getRank() << endl;

    mpAssyModel= pAssyModel;

    for(uiint imesh=0; imesh < mvVector.size(); imesh++) delete mvVector[imesh];
    mvVector.clear();

    // Mesh別Vectorの生成
    for(uiint imesh=0; imesh < pAssyModel->getNumOfMesh(); imesh++) {
        CMesh *pMesh= pAssyModel->getMesh(imesh);
        uiint nMeshID= pMesh->getMeshID();

        CVector *pVector = new CVector( pAssyModel->getMesh(imesh), vDOF[nMeshID] );//-- グローバルモデルでは,MeshIDは通し番号
        mvVector.push_back(pVector);
    };

    // 通信オーバーラップのカウント
    for(uiint imesh=0; imesh < pAssyModel->getNumOfMesh(); imesh++) {
        CMesh *pMesh= pAssyModel->getMesh(imesh);
        CIndexBucket* pBucket= pMesh->getBucket();

        uiint nNumOfCommMesh= pMesh->getCommMesh2Size();
        for(uiint icom=0; icom < nNumOfCommMesh; icom++) {
            CCommMesh2 *pCommMesh2= pMesh->getCommMesh2IX(icom);

            uiint nNumOfCommNode= pCommMesh2->getCommNodeSize();
            for(uiint icnode=0; icnode < nNumOfCommNode; icnode++) {
                CCommNode *pCommNode= pCommMesh2->getCommNodeIX(icnode);

                uiint node_id= pCommNode->getNodeID();
                uiint inode= pBucket->getIndexNode(node_id);

                mvVector[imesh]->increment_OverlapComm(inode);//---- 通信オーバーラップ数++
            };
        };
    };

    ////cout << "-- AssyVector::initAssyVector ------ exit -- rank:" << pMPI->getRank() << endl;
}

void CAssyVector::initAssyVector(const CAssyVector* pAssyVector)
{
    mpAssyModel = pAssyVector->mpAssyModel;

    for(uiint imesh=0; imesh < mvVector.size(); imesh++) delete mvVector[imesh];
    mvVector.clear();

    for(uiint imesh=0; imesh < pAssyVector->getNumOfVector(); imesh++) {
        CVector *pVector = new CVector( pAssyVector->getVector(imesh) );// 通信オーバーラップ数:CVectorコピーコンストラクターで処理
        mvVector.push_back(pVector);
    };
}

//
// 引数 vDOF:グローバルMesh番号(MeshID)に対応するDOF
//
CAssyVector::CAssyVector(CAssyModel *pAssyModel, vuint& vDOF)
    : mpAssyModel( pAssyModel )
{
    initAssyVector(pAssyModel, vDOF);//--------------- 初期化関数
}

CAssyVector::CAssyVector(const CAssyVector *pAssyVector)
{
    initAssyVector(pAssyVector);//-------------------- 初期化関数
}

CAssyVector::~CAssyVector()
{
    for_each(mvVector.begin(), mvVector.end(), DeleteObject());
}

uiint CAssyVector::size() const
{
    uiint sum = 0;
    for (CVVecConstIter icv = mvVector.begin(); icv != mvVector.end(); icv++) {
        sum += (*icv)->size();
    }
    return sum;
}
uiint CAssyVector::size()
{
    uiint sum = 0;
    for (CVVecConstIter icv = mvVector.begin(); icv != mvVector.end(); icv++) {
        sum += (*icv)->size();
    }
    return sum;
}
uiint CAssyVector::getNumOfMesh() const
{
    return mvVector.size();
}
uiint CAssyVector::getNumOfMesh()
{
    return mvVector.size();
}
uiint& CAssyVector::getDOF(const uiint& imesh) const
{
    return mvVector[imesh]->getDOF();
}
uiint& CAssyVector::getDOF(const uiint& imesh)
{
    return mvVector[imesh]->getDOF();
}

const CVector::ElemType &CAssyVector::operator[](uiint idx) const
{
    for (CVVecConstIter icv = mvVector.begin(); icv != mvVector.end(); icv++) {
        uiint len = (*icv)->size();
        if (idx < len) return (**icv)[idx];
        idx -= len;
    }
}
CVector::ElemType &CAssyVector::operator[](uiint idx)
{
    for (CVVecIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
        uiint len = (*iv)->size();
        if (idx < len) return (**iv)[idx];
        idx -= len;
    }
}
const double &CAssyVector::operator()(uiint meshID, uiint nodeID, uiint dof) const
{
    return (*mvVector[meshID])[nodeID](dof);
}
double &CAssyVector::operator()(uiint meshID, uiint nodeID, uiint dof)
{
    return (*mvVector[meshID])[nodeID](dof);
}
void CAssyVector::Vector_Clear(const uiint& iMesh)
{
    mvVector[iMesh]->Vector_Clear();
}
void CAssyVector::Vector_Clear()
{
    for(uiint i=0; i < mvVector.size(); i++) {
        mvVector[i]->Vector_Clear();
    }
}
void CAssyVector::setZero()
{
    for (CVVecIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
        (*iv)->setZero();
    }
}
void CAssyVector::setValue(uiint imesh, uiint inode, uiint idof, double value)
{
    mvVector[imesh]->setValue(inode, idof, value);
}
void CAssyVector::addValue(const uiint& imesh, const uiint& inode, const uiint& idof, const double& value)
{
    mvVector[imesh]->addValue(inode, idof, value);
}
double& CAssyVector::getValue(uiint imesh, uiint inode, uiint idof)
{
    return( mvVector[imesh]->getValue(inode, idof) );
}
// 代入
void CAssyVector::set(const CAssyVector* pV)
{
    CVVecConstIter icv = pV->mvVector.begin();
    for (CVVecIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
        (*iv)->set((*icv));
        icv++;
    }
}
void CAssyVector::add(const CAssyVector *pV)
{
    CVVecConstIter icv = pV->mvVector.begin();
    for (CVVecIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
        (*iv)->add((*icv));
        icv++;
    }
}
// 代入
void CAssyVector::subst(const CAssyVector *pV)
{
    CVVecConstIter icv = pV->mvVector.begin();
    for (CVVecIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
        (*iv)->subst((*icv));
        icv++;
    }
}
//減算
void CAssyVector::subtrac(const CAssyVector *pV)
{
    for(uiint i=0; i < mvVector.size(); i++) {
        mvVector[i]->subtrac( pV->mvVector[i] );
    };
}

double CAssyVector::norm2() const
{
    double sum = 0;
    for (CVVecConstIter icv = mvVector.begin(); icv != mvVector.end(); icv++) {
        sum += (*icv)->norm2();
    }
    pmw::CHecMPI *pMPI= pmw::CHecMPI::Instance();
    double sum0 = sum;
    pMPI->Allreduce(&sum0, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sum;
}
double CAssyVector::innerProd(const CAssyVector *pX) const
{
    double sum = 0;
    CVVecConstIter icX = pX->mvVector.begin();
    for (CVVecConstIter icv = mvVector.begin(); icv != mvVector.end(); icv++) {
        sum += (*icv)->innerProd((*icX));
        icX++;
    }
    pmw::CHecMPI *pMPI= pmw::CHecMPI::Instance();
    double sum0 = sum;
    pMPI->Allreduce(&sum0, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sum;
}
//--
// Update
//--
void CAssyVector::updateCommBoundary()
{
    //--
    // グローバル通信テーブル
    //--
    uiint bufsize(0);
    for(uiint i=0; i < mvVector.size(); i++) {
        uiint vecDOF= mvVector[i]->getDOF();
        uiint tempsize= mvVector[i]->size() * vecDOF;
        if(bufsize < tempsize) bufsize= tempsize;
    };
    bufsize *= mvVector.size();
    double* buffer = (double*)malloc(sizeof(double) * bufsize);

    pmw::CHecMPI *pMPI= pmw::CHecMPI::Instance();
    MPI_Status stat;
    MPI_Request req;
    uiint myrank= pMPI->getRank();

    //--
    // icomm_id : GlobalCommID Number
    //--
    for(uiint icomm_id=0; icomm_id < mpAssyModel->getNumOfGlobalCommMesh2(); icomm_id++) {

        if( myrank == mpAssyModel->getGlobalPairRank_1st(icomm_id) || myrank==mpAssyModel->getGlobalPairRank_2nd(icomm_id) ) {

            uiint nMeshID= mpAssyModel->getMeshID_with_CommID(icomm_id);
            uiint imesh= mpAssyModel->getIndex_of_Mesh(nMeshID);
            ////////////////debug
            //cout << "AssyVector::updateCommBoundary  MeshID:" << nMeshID << " imesh:" << imesh << " icomm_id:" << icomm_id
            //        << " rank:" << pMPI->getRank() << endl;

            CVector* pVec= mvVector[imesh];
            uiint nDOF = pVec->getDOF();

            CMesh *pMesh= mpAssyModel->getMesh(imesh);
            CIndexBucket *pBucket= pMesh->getBucket();

            CCommMesh2 *pCommMesh2= pMesh->getCommMesh2(icomm_id);
            iint source      = pCommMesh2->getRank();
            iint destination = pCommMesh2->getTrasmitRank();
            uiint nNumOfCommNode= pCommMesh2->getCommNodeSize();

            if( nNumOfCommNode*nDOF > IINT_MAX ) { //---------- サイズ・エラー 2012.01.04
                Utility::CLogger *pLogger= Utility::CLogger::Instance();
                pLogger->Info(Utility::LoggerMode::Error, "AssyVector::updateCommBoundary, over IINT_MAX");
            }

            if( source < destination ) {
                // 送信
                uiint ic= 0;
                for(uiint icnode=0; icnode < nNumOfCommNode; icnode++) {
                    CCommNode *pCommNode= pCommMesh2->getCommNodeIX(icnode);
                    uiint nNodeID= pCommNode->getNode()->getID();
                    uiint inode = pBucket->getIndexNode(nNodeID);

                    ////cout << "update 送信 icom_node:" << icnode << " inode:" << inode << " node_id:" << nNodeID
                    ////        << " rank:" << pMPI->getRank() << endl;

                    for(uiint idof=0; idof < nDOF; idof++) {
                        buffer[ic++] = pVec->getValue(inode, idof);//------------ rank小:値を与える.
                    };
                };//icnode loop end

                //pMPI->Send(buffer, nNumOfCommNode*nDOF, MPI_DOUBLE, destination, 100+icomm_id, MPI_COMM_WORLD);
                pMPI->Isend(buffer, nNumOfCommNode*nDOF, MPI_DOUBLE, destination, 100+icomm_id, MPI_COMM_WORLD, &req);
                pMPI->Wait(&req, &stat);
            } else {
                // 受信
                //pMPI->Recv(buffer, nNumOfCommNode*nDOF, MPI_DOUBLE, destination, 100+icomm_id, MPI_COMM_WORLD, &stat);
                pMPI->Irecv(buffer, nNumOfCommNode*nDOF, MPI_DOUBLE, destination, 100+icomm_id, MPI_COMM_WORLD, &req);
                pMPI->Wait(&req, &stat);

                uiint ic = 0;
                for(uiint icnode=0; icnode< nNumOfCommNode; icnode++) {
                    CCommNode *pCommNode= pCommMesh2->getCommNodeIX(icnode);
                    uiint nNodeID= pCommNode->getNode()->getID();
                    uiint inode = pBucket->getIndexNode(nNodeID);

                    ////cout << "update 受信 icom_node:" << icnode << " inode:" << inode << " node_id:" << nNodeID
                    ////        << " rank:" << pMPI->getRank() << endl;

                    for(uiint idof=0; idof < nDOF; idof++) {
                        pVec->setValue(inode, idof, buffer[ic++]);//------------- rank大:値を代入(更新)
                    };
                };
            }//---- if(source < destination) end

        }// if(myrank==rank_pair) end

    };// NumGlobalComm loop end

    //pMPI->Barrier(MPI_COMM_WORLD);

    if(bufsize > 0) {
        if(buffer) free(buffer);
        buffer=NULL;
    }

}
//--
// Sumup
//--
void CAssyVector::sumupCommBoundary()
{
    //--
    // グローバル通信テーブル
    //--
    uiint bufsize(0);
    for(uiint i=0; i < mvVector.size(); i++) {
        uiint vecDOF= mvVector[i]->getDOF();
        uiint tempsize= mvVector[i]->size() * vecDOF;
        if(bufsize < tempsize) bufsize= tempsize;
    };
    bufsize *= mvVector.size();
    double* buffer = (double*)malloc(sizeof(double) * bufsize);

    pmw::CHecMPI *pMPI= pmw::CHecMPI::Instance();

    MPI_Status stat;
    MPI_Request req;
    uiint myrank = pMPI->getRank();
    //--
    // icomm_id : GlobalCommID Number
    //--
    for(uiint icomm_id=0; icomm_id < mpAssyModel->getNumOfGlobalCommMesh2(); icomm_id++) {

        if(myrank == mpAssyModel->getGlobalPairRank_1st(icomm_id) || myrank==mpAssyModel->getGlobalPairRank_2nd(icomm_id)) {

            uiint nMeshID= mpAssyModel->getMeshID_with_CommID(icomm_id);
            uiint imesh= mpAssyModel->getIndex_of_Mesh(nMeshID);
            ////////////////debug
            ////cout << "AssyVector::sumupCommBoundary  MeshID:" << nMeshID << " imesh:" << imesh << " icomm_id:" << icomm_id
            ////        << " rank:" << pMPI->getRank() << endl;

            CVector* pVec= mvVector[imesh];
            uiint nDOF = pVec->getDOF();

            CMesh *pMesh= mpAssyModel->getMesh(imesh);
            CIndexBucket *pBucket= pMesh->getBucket();

            CCommMesh2 *pCommMesh2= pMesh->getCommMesh2(icomm_id);
            iint source      = pCommMesh2->getRank();
            iint destination = pCommMesh2->getTrasmitRank();
            uiint nNumOfCommNode= pCommMesh2->getCommNodeSize();

            ////cout << "sumup   source:" << source << " destination:" << destination << "  NumOfCommNode:" << nNumOfCommNode << endl;

            if( nNumOfCommNode*nDOF > IINT_MAX ) { //---------- サイズ・エラー 2012.01.04
                Utility::CLogger *pLogger= Utility::CLogger::Instance();
                pLogger->Info(Utility::LoggerMode::Error, "AssyVector::sumupCommBoundary, over IINT_MAX");
            }

            if( source > destination ) {
                // 送信
                uiint ic= 0;
                for(uiint icnode=0; icnode< nNumOfCommNode; icnode++) {
                    CCommNode *pCommNode= pCommMesh2->getCommNodeIX(icnode);
                    uiint nNodeID = pCommNode->getNode()->getID();
                    uiint inode = pBucket->getIndexNode(nNodeID);

                    ////cout << "sumup 送信 icom_node:" << icnode << " inode:" << inode << " node_id:" << nNodeID
                    ////        << " rank:" << pMPI->getRank() << endl;

                    for(uiint idof=0; idof < nDOF; idof++) {
                        buffer[ic++] = pVec->getValue(inode, idof);
                        pVec->setValue(inode, idof, 0.0);//---------------------- rank大:0.0
                    };
                };
                ////cout << " sumup buffer send" << endl;
                ////for(uiint i=0; i < nNumOfCommNode*nDOF; i++) cout << " i:" << buffer[i] ;
                ////cout << endl;

                //pMPI->Send(buffer, nNumOfCommNode*nDOF, MPI_DOUBLE, destination, 100+icomm_id, MPI_COMM_WORLD);
                pMPI->Isend(buffer, nNumOfCommNode*nDOF, MPI_DOUBLE, destination, 100+icomm_id, MPI_COMM_WORLD, &req);
                pMPI->Wait(&req, &stat);
            } else {
                // 受信
                //pMPI->Recv(buffer, nNumOfCommNode*nDOF, MPI_DOUBLE, destination, 100+icomm_id, MPI_COMM_WORLD, &stat);
                pMPI->Irecv(buffer, nNumOfCommNode*nDOF, MPI_DOUBLE, destination, 100+icomm_id, MPI_COMM_WORLD, &req);
                pMPI->Wait(&req, &stat);

                ////cout << " sumup buffer recv" << endl;
                ////for(uiint i=0; i < nNumOfCommNode*nDOF; i++) cout << " i:" << buffer[i] ;
                ////cout << endl;

                uiint ic= 0;
                for(uiint icnode=0; icnode< nNumOfCommNode; icnode++) {
                    CCommNode *pCommNode= pCommMesh2->getCommNodeIX(icnode);
                    uiint nNodeID= pCommNode->getNode()->getID();
                    uiint inode = pBucket->getIndexNode(nNodeID);

                    ////cout << "sumup 受信 icom_node:" << icnode << " inode:" << inode << " node_id:" << nNodeID
                    ////        << " rank:" << pMPI->getRank() << endl;

                    for(uiint idof=0; idof < nDOF; idof++) {
                        pVec->addValue(inode, idof, buffer[ic++]);//------------- rank小:加算
                    };
                };
            }//---- if(source > destination) end

        }// if(myrank==rank_pair) end

    };// NumGlobalComm loop end

    //pMPI->Barrier(MPI_COMM_WORLD);

    if( bufsize > 0 ) {
        if(buffer) free(buffer);
        buffer=NULL;
    }

}



//--
// restrict
//--
uiint CAssyVector::restrictTo(CAssyVector *pVc) const
{
    CVVecIter ivc = pVc->mvVector.begin();
    for (CVVecConstIter icv = mvVector.begin(); icv != mvVector.end(); icv++) {
        (*icv)->restrictTo((*ivc));//------ Vector::restrictTo
        ivc++;
    }
    return 1;
}
//--
// prolongate
//--
uiint CAssyVector::prolongateFrom(const CAssyVector *pcVc)
{
    CVVecConstIter icvc = pcVc->mvVector.begin();
    for (CVVecIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
        (*iv)->prolongateFrom((*icvc));//------ Vector::prolongateFrom
        icvc++;
    }
    return 1;
}
uiint CAssyVector::getNumOfVector() const
{
    return mvVector.size();
}
const CVector *CAssyVector::getVector(uiint index) const
{
    return mvVector[index];
}
CVector *CAssyVector::getVector(uiint index)
{
    return mvVector[index];
}
void CAssyVector::dump() const
{
    uiint nNumOfParts = mvVector.size();// Mesh数
    uiint ipart;
    for(ipart=0; ipart < nNumOfParts; ipart++) {
        cout << " ---- iParts : " << ipart << " Vector ---- " << endl;
        mvVector[ipart]->dump();
    };
    cout << endl;
}
void CAssyVector::dump()
{
    CHecMPI *pMPI= CHecMPI::Instance();

    uiint nNumOfParts = mvVector.size();// Mesh数
    uiint ipart;
    for(ipart=0; ipart < nNumOfParts; ipart++) {
        cout << " ---- iParts : " << ipart << " Vector ---- rank:" << pMPI->getRank() << endl;
        mvVector[ipart]->dump();
    };
    cout << endl;
}
//--
// ベクトル要素を通信オーバーラップ数で除算
//--
void CAssyVector::divis_OverlapNum()
{
    uiint nNumOfMesh= mvVector.size();
    for(uiint imesh=0; imesh < nNumOfMesh; imesh++) {
        mvVector[imesh]->divis_OverlapNum();
    };
}

}
