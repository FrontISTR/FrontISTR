/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/AssyMatrix.cpp
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
#include "AssyMatrix.h"
#include "MatrixBCRS.h"
#include "ContactMatrix.h"
#include "MPCMatrix.h"
#include "AssyVector.h"
#include "Solver.h"
#include "AssyModel.h"
#include "Equation.h"
#include "SolverCG.h"
#include "SolverGMRES.h"

#include "Logger.h"

namespace pmw
{
typedef std::vector<CMatrixBCRS*> CVMat;
typedef CVMat::iterator CVMatIter;
typedef CVMat::const_iterator CVMatConstIter;
//--
//
//--
CAssyMatrix::CAssyMatrix(CAssyModel *pAssyModel, vuint& vDOF, const uiint& ieq)
{

    Utility::CLogger *pLogger= Utility::CLogger::Instance();

    CHecMPI* pMPI= CHecMPI::Instance();

    MPI_Datatype MPI_UIINT= pMPI->MPI_UIINT();
    MPI_Datatype MPI_IINT = pMPI->MPI_IINT();
    iint tag(200);
    MPI_Status stat;
    MPI_Request req[2];

////    //debug file.
////    string fname("debug_assy_matrix_");
////    stringstream ss; ss << pMPI->getRank();
////    fname += ss.str();
////    fname += ".txt";
////    ofstream dofs(fname.c_str());


    mpAssyModel = pAssyModel;
    uiint nNumOfMesh = pAssyModel->getNumOfMesh();

    //--
    // # vDOFは、グローバル番号(MeshID)ごとのDOF
    // # MeshIDは、グローバルインデックス番号
    // # nNumOfMesh は、ローカルMesh数
    //--
    for(uiint i=0; i < nNumOfMesh; i++) {
        CMesh *pMesh = pAssyModel->getMesh(i);
        uiint nMeshID= pMesh->getMeshID();

        CMatrixBCRS *pMatrix = new CMatrixBCRS(pMesh, vDOF[nMeshID]);
        mvMatrix.push_back( pMatrix );
    }
    mMGLevel = pAssyModel->getMGLevel();
    uiint nNumOfContact= pAssyModel->getNumOfContactMesh();
    uiint idof;

    CEquation *vEquation;//--[スレーブ点数*DOF数]確保 : スレーブ点1個に対して、DOF順にEquationが存在


    // ContactMesh 数ぶんのMPCMatrix を生成
    for(uiint icont = 0; icont < nNumOfContact; icont++) {
        CContactMesh* pConMesh= pAssyModel->getContactMesh(icont);
        uiint nNumOfSPoint = pConMesh->getNumOfSlavePoint();

        ////pLogger->Info_format(Utility::LoggerMode::Info, "%s%d%s%d", "CAssyMatrix::CAssyMatrix, rank:", pMPI->getRank()," NumOfSPoint:", nNumOfSPoint);

        CMPCMatrix* mpc = new CMPCMatrix();//--------- MPCMatrix生成 (ContactMesh数ぶん)

        //--
        // Ax=b方程式のマスターとスレーブの小さいほうのDOF : 2つのMeshパーツの通し番号を元にvDOFから自由度を取得
        //--
        uiint nMasterDOF, nSlaveDOF, nLessDOF;
        CContactNode* psCNode = pConMesh->getSlaveConNode(0);
        CContactNode* pmCNode = pConMesh->getMasterConNode(0);
        uiint sMeshID= psCNode->getMeshID();// スレーブのMeshグローバルインデックス:MeshID
        uiint mMeshID= pmCNode->getMeshID();// マスターのMeshグローバルインデックス:MeshID

        nSlaveDOF = vDOF[sMeshID]; // スレーブAlgebraEquation(線形方程式)のDOF
        nMasterDOF= vDOF[mMeshID]; // マスターAlgebraEquation(線形方程式)のDOF

        bool bSurplusDOF(false);
        if(nMasterDOF < nSlaveDOF) {
            // マスターの方がDOFが少ない
            nLessDOF= nMasterDOF;
            bSurplusDOF= true;
        } else {
            // スレーブの方がDOFが少ない, もしくは両者等しい
            nLessDOF= nSlaveDOF;
            if(nMasterDOF > nSlaveDOF) {
                pLogger->Info_format(Utility::LoggerMode::Info, "%s%d", "CAssyMatrix::CAssyMatrix, master_dof >= slave_dof,  rank:", pMPI->getRank());
            }
        }
        // ---- マスター面の法線方向に固定になっていないので、コメントアウト ----
        //////    //--
        //////    // スレーブのDOF数が多い場合に、余った自由度を固定:(対象: ieq番目のAx=b:線形方程式)
        //////    //--
        //////    if(bSurplusDOF){
        //////        CAssyVector *pRHSAssyVec= pAssyModel->getRHSAssyVector(ieq);
        //////        double fixValue= 0.0, diagValue= 1.0;
        //////
        //////        for(uiint islave=0; islave < nNumOfSPoint; islave++){
        //////            CContactNode* pSlaveNode = pConMesh->getSlaveConNode(islave);
        //////            //--
        //////            // 自身Rank内のスレーブ点のみ処理
        //////            //--
        //////            if(pSlaveNode->getRank()==pMPI->getRank()){
        //////                //--
        //////                // マスター面に載っているスレーブ点
        //////                //--
        //////                if(pSlaveNode->have_MasterFaceID(mMGLevel)){
        //////                    uiint sMeshID= pSlaveNode->getMeshID();
        //////                    //--
        //////                    // 自身ランクにMeshが存在する場合に余分DOFを拘束
        //////                    //--
        //////                    bool bSelfRank= pAssyModel->isSelfMesh(sMeshID);
        //////
        //////                    if( bSelfRank ){
        //////                        CMesh* pSMesh = pAssyModel->getMesh_ID(sMeshID);
        //////                        CIndexBucket *pSBucket = pSMesh->getBucket();
        //////
        //////                        uiint node_id= pSlaveNode->getNodeID();
        //////                        uiint inode  = pSBucket->getIndexNode(node_id);
        //////                        uiint sMeshIX= pAssyModel->getIndex_of_Mesh(sMeshID);
        //////
        //////                        for(uiint idof=nLessDOF; idof < nSlaveDOF; idof++){
        //////                            //--
        //////                            // ## CMW::Set_BC_Mat_RHS2 ## Dirichlet Boundary
        //////                            //--
        //////                            pRHSAssyVec->setValue(sMeshIX, inode, idof, fixValue);
        //////                            this->setZero_NonDiag(sMeshIX, inode, idof, pRHSAssyVec, fixValue);
        //////                            this->setValue(sMeshIX, inode, idof, diagValue, pRHSAssyVec, fixValue);
        //////                        };
        //////                    }//if(bSelfRank) end
        //////                }//if(pSlave::have_MasterFaceID) end
        //////            }//if(myRank) end
        //////        };
        //////    }//if(bSurplusDOF) end

        // --
        // MPC_Equation配列：スレーブ点1個に対して、DOF順にEquationを生成
        // --
        vEquation= new CEquation[nNumOfSPoint * nLessDOF];//------------ Equation生成 (スレーブ点数 * DOF)

        for(uiint islave = 0; islave< nNumOfSPoint; islave++) {
            CContactNode* pSlaveNode = pConMesh->getSlaveConNode(islave);

            if(pSlaveNode->have_MasterFaceID(mMGLevel)) {
                uiint smesh=IINT_MAX, snode=IINT_MAX;
                uiint s_commrank= pSlaveNode->getRank();//---スレーブRank(通信)
                uiint s_solvrank= pSlaveNode->getRank();//---スレーブRank(計算)
                uiint s_overlap = 1;//-------------------- 初期値:スレーブ点オーバーラップ数
                CMesh* pSMesh;
                CIndexBucket *pSBucket;
                uiint smesh_ix=IINT_MAX, snode_ix=IINT_MAX;
                pSlaveNode->sort_OverlapRank();//--------------- vRankソート
                vuint s_vRank= pSlaveNode->getOverlapRank();
                uiint nsVal[4];//--通信用テンポラリー変数
                //--
                // スレーブ点
                //--
                // オーバーラップ点のrank更新
                if( pMPI->getRank()==s_solvrank && pSlaveNode->isOverlap() ) {
                    if( pSlaveNode->getRank() >= pSlaveNode->getOverlapMinRank() ) {
                        s_commrank= pSlaveNode->getOverlapMinRank();
                    }
                    s_overlap= pSlaveNode->getNumOfOverlapRank();//--- オーバーラップ数

                    ////debug
                    //dofs << "スレーブ オーバーラップ MeshID:" << pSlaveNode->getMeshID() << " NodeID:" << pSlaveNode->getNodeID() << " rank:" << pMPI->getRank() << endl;
                }
                //if(pMPI->getRank()==s_commrank){
                if(pMPI->getRank()==s_solvrank) {
                    smesh = pSlaveNode->getMeshID();
                    pSMesh= pAssyModel->getMesh_ID(smesh);//----- スレーブMesh
                    pSBucket= pSMesh->getBucket();        //----- スレーブMeshのBucket
                    smesh_ix= pAssyModel->getIndex_of_Mesh(smesh);//----ローカルMeshインデックス
                    snode = pSlaveNode->getNodeID();
                    snode_ix = pSBucket->getIndexNode(snode);//-- スレーブ節点Index番号

                    nsVal[0]=snode_ix;
                    nsVal[1]=smesh_ix;
                    nsVal[2]=snode;
                    nsVal[3]=smesh;

                    ////debug
                    //dofs << "スレーブsolvrank MeshID:" << smesh << " NodeID:" << snode
                    //        << " solvrank:" << s_solvrank << " commrank:" << s_commrank << endl;
                }
                pMPI->Bcast((void*)nsVal, 4, MPI_UIINT, s_commrank, MPI_COMM_WORLD);

                snode_ix= nsVal[0];
                smesh_ix= nsVal[1];
                snode= nsVal[2];
                smesh= nsVal[3];

                ////debug
                //dofs << "スレーブBcast MeshID:" << smesh << " NodeID:" << snode
                //        << " solvrank:" << s_solvrank << " commrank:" << s_commrank << endl;

                //--
                // #Term(スレーブ): スレーブ点1個  : islaveのidofのMPC_Equationにセット
                //--
                for(idof=0; idof < nLessDOF; idof++)
                    vEquation[islave*nLessDOF + idof].setSlave(smesh, smesh_ix, snode, snode_ix, idof, s_commrank, s_solvrank, s_overlap, s_vRank);

                //--
                // マスター面
                //--
                uiint masterFaceID = pSlaveNode->getMasterFaceID(mMGLevel);
                CSkinFace* pMasterFace = pConMesh->getMasterFace_ID(masterFaceID);
                uiint nNumOfVert = pMasterFace->getNumOfNode();

////                            //debug
////                            if(pMasterFace->isSelf()){
////                                dofs << "--------- self master_face id:" << pMasterFace->getID() << endl;
////                            }else{
////                                dofs << "----- non-self master_face id:" << pMasterFace->getID() << endl;
////                            }

                // マスター面の頂点 ループ
                for(uiint ivert=0; ivert< nNumOfVert; ivert++) {
                    CContactNode *pMasterCNode= pMasterFace->getNode(ivert);
                    uiint mmesh=IINT_MAX, mnode=IINT_MAX;
                    uiint m_commrank= pMasterCNode->getRank();//--- マスター点Rank(通信)
                    uiint m_solvrank= pMasterCNode->getRank();//--- マスター点Rank(計算)
                    uiint m_overlap = 1;//--------------------- 初期値:マスター点オーバーラップ数
                    CMesh* pMMesh;
                    CIndexBucket *pMBucket;
                    uiint mmesh_ix=IINT_MAX, mnode_ix=IINT_MAX;//-- 初期値
                    pMasterCNode->sort_OverlapRank();//------------ vRankソート
                    vuint m_vRank= pMasterCNode->getOverlapRank();
                    uiint nmVal[4];//--通信用テンポラリー変数
                    //--
                    // マスター点
                    //--
                    // オーバーラップ点のrank更新
                    if( pMPI->getRank()==m_solvrank && pMasterCNode->isOverlap() ) {
                        if( pMasterCNode->getRank() >= pMasterCNode->getOverlapMinRank() ) {
                            m_commrank= pMasterCNode->getOverlapMinRank();
                        }
                        m_overlap= pMasterCNode->getNumOfOverlapRank();//--- オーバーラップ数

                        ////debug
                        //dofs << "マスター オーバーラップ MeshID:" << pMasterCNode->getMeshID() << " NodeID:" << pMasterCNode->getNodeID() << " rank:" << pMPI->getRank() << endl;
                    }
                    //if(pMPI->getRank()==m_commrank){
                    if(pMPI->getRank()==m_solvrank) {
                        mmesh = pMasterCNode->getMeshID();//----- マスターMeshID
                        pMMesh= pAssyModel->getMesh_ID(mmesh);//----- マスターMesh
                        pMBucket= pMMesh->getBucket();        //----- マスターMeshのBucket
                        mmesh_ix= pAssyModel->getIndex_of_Mesh(mmesh);//---ローカルMeshインデックス
                        mnode = pMasterCNode->getNodeID();
                        mnode_ix = pMBucket->getIndexNode(mnode);//-- マスター節点のIndex番号

                        nmVal[0]= mnode_ix;
                        nmVal[1]= mmesh_ix;
                        nmVal[2]= mnode;
                        nmVal[3]= mmesh;

////                                        //debug
////                                        dofs << "マスターsolvrank MeshID:" << mmesh << " NodeID:" << mnode
////                                                << " solvrank:" << m_solvrank << " commrank:" << m_commrank << endl;
                    }
                    pMPI->Bcast((void*)nmVal, 4, MPI_UIINT, m_commrank, MPI_COMM_WORLD);

                    mnode_ix= nmVal[0];
                    mmesh_ix= nmVal[1];
                    mnode= nmVal[2];
                    mmesh= nmVal[3];

////                                    //debug
////                                    dofs << "マスターBcast MeshID:" << mmesh << " NodeID:" << mnode
////                                            << " solvrank:" << m_solvrank << " commrank:" << m_commrank << endl;

                    uiint slaveID = pSlaveNode->getID();//----------------- スレーブ点ID
                    double coef = pMasterFace->getCoef(slaveID, ivert);//-- マスター面の各頂点のCoef
                    coef *= pConMesh->getTransCoeff(ieq);//---------------- MPC係数に伝達率

                    //--
                    // #Term(マスター)の生成:マスター面の各頂点の係数  : islaveのidofのMPC_Equationにセット
                    //--
                    for(idof=0; idof < nLessDOF; idof++)
                        vEquation[ islave*nLessDOF + idof].addMaster(mmesh, mmesh_ix, mnode, mnode_ix, idof, coef, m_commrank, m_solvrank, m_overlap, m_vRank);

                };
                // --
                // 接合Mesh毎のMPCMatrixに、全スレーブ点(全DOF)のMPC_Equation配列を保存
                // --
                // MPC_Equation単体の内容:
                // - Term(スレーブ) : スレーブ点 1点のスレーブCoef(=1.0)
                // - Term(マスター) : スレーブ点 1点に対して、マスター面頂点数 のマスターCoef
                for(idof=0; idof < nLessDOF; idof++)
                    mpc->addEquation(&vEquation[islave*nLessDOF + idof]);

            }//end if(マスター面を所有している点)
        };//end スレーブ点 ループ

        mvMPCMatrix.push_back(mpc);//----- MPCMatrix 保存

    };//end 接合Mesh ループ

////    dofs.close();//debug file.
////
////    pLogger->Info_format(Utility::LoggerMode::Info, "%s%d%s%d", "AssyMatrix exit rank:", pMPI->getRank(), " mgLevel:", pAssyModel->getMGLevel());

    mpCoarseMatrix=NULL;
    mpSolver=NULL;
}
CAssyMatrix::~CAssyMatrix()
{
    for_each(mvMatrix.begin(), mvMatrix.end(), DeleteObject());
    for_each(mvMPCMatrix.begin(), mvMPCMatrix.end(), DeleteObject());

    if(mpSolver) delete mpSolver;
    mpSolver=NULL;
    mpCoarseMatrix=NULL;
}


uiint CAssyMatrix::Matrix_Add_Nodal(const uiint& iMesh, const uiint& inode, const uiint& jnode, double* NodalMatrix)
{
    mvMatrix[iMesh]->Matrix_Add_Nodal(inode, jnode, NodalMatrix);
    return 1;
}
uiint CAssyMatrix::Matrix_Add_Elem(CAssyModel *pAssyModel, const uiint& iMesh, const uiint& iElem, double *ElemMatrix)
{
    CMesh *pMesh = pAssyModel->getMesh(iMesh);
    mvMatrix[iMesh]->Matrix_Add_Elem(pMesh, iElem, ElemMatrix);

    return 1;
}
void CAssyMatrix::Matrix_Clear(const uiint& iMesh)
{
    mvMatrix[iMesh]->Matrix_Clear();
}
void CAssyMatrix::Matrix_Clear()
{
    for(uiint i=0; i < mvMatrix.size(); i++) {
        mvMatrix[i]->Matrix_Clear();
    };
}

void CAssyMatrix::setValue_D(const uiint& imesh, const uiint& inode, const uiint& idof, const double& value)
{
    mvMatrix[imesh]->setValue_D(inode, idof, value);
}
void CAssyMatrix::setValue(const uiint& imesh, const uiint& inode, const uiint& idof, const double& dDiag, CAssyVector *pRHS, const double& dRHS)
{
    mvMatrix[imesh]->setValue(inode, idof, dDiag, pRHS->getVector(imesh), dRHS);
}
void CAssyMatrix::setZero_NonDiag(const uiint& imesh, const uiint& inode, const uiint& idof, CAssyVector *pRHS, const double& dRHS)
{
    mvMatrix[imesh]->setZero_NonDiag(inode, idof, pRHS->getVector(imesh), dRHS);
}


//--
// Av=p 並列通信あり : 各種ソルバーが使用するmultVector
//--
uiint CAssyMatrix::multVector(CAssyVector *pV, CAssyVector *pP, CAssyVector *pW) const
{
    //
    // 各種ソルバーが使用するmultVector
    //
    bool bFlagW(false);

////    //debug
////    CHecMPI *pMPI= CHecMPI::Instance();
////    Utility::CLogger *pLogger= Utility::CLogger::Instance();

    pV->updateCommBoundary();//------通信: 更新 x-ベクトル

//    //--
//    //debug
//    //--
//    cout << "------ AssyMatrix::multVector update-X 左辺ベクトル ------" << endl;
//    pV->dump();//左辺ベクトル
//    cout << "------ AssyMatrix::multVector update-X 右辺ベクトル ------" << endl;
//    pP->dump();//結果

////    pLogger->Info_format(Utility::LoggerMode::Info, "%s%d%s%d%s", "AssyMatrix::multVector  rank:", pMPI->getRank(), " mgLevel:", mpAssyModel->getMGLevel(), " ----- updateComm X");

    if(mvMPCMatrix.size()==0) {
        //標準
        for (uiint i = 0; i < mvMatrix.size(); i++) {
            mvMatrix[i]->multVector(pV->getVector(i), pP->getVector(i));// A*v=P
        };
    } else {
        if(pW == NULL) {
            pW = new CAssyVector(pP);
            bFlagW=true;
        }

        //MPC
        for (uiint i = 0; i < mvMPCMatrix.size(); i++) {
            mvMPCMatrix[i]->multVector(pV, pP);//------------------------------Tu  => pP
        };
        for (uiint i = 0; i < mvMatrix.size(); i++) {
            mvMatrix[i]->multVector(pP->getVector(i), pW->getVector(i));//-----KTu => pW
        };
        for (uiint i = 0; i < mvMPCMatrix.size(); i++) {
            mvMPCMatrix[i]->transMultVector(pW, pP);//-------------------(転置T)KTu => pP
        };
    }
////    pLogger->Info_format(Utility::LoggerMode::Info, "%s%d%s%d%s", "AssyMatrix::multVector  rank:", pMPI->getRank(), " mgLevel:", mpAssyModel->getMGLevel(), " ----- multVector");

    pP->sumupCommBoundary(); //------通信: 加算 y-ベクトル

//    //--
//    //debug
//    //--
//    cout << "------ AssyMatrix::multVector sumup 左辺ベクトル ------" << endl;
//    pV->dump();//左辺ベクトル
//    cout << "------ AssyMatrix::multVector sumup 右辺ベクトル ------" << endl;
//    pP->dump();//結果


////    pLogger->Info_format(Utility::LoggerMode::Info, "%s%d%s%d%s", "AssyMatrix::multVector  rank:", pMPI->getRank(), " mgLevel:", mpAssyModel->getMGLevel(), " ----- updateComm Y");

    if(bFlagW) delete pW;

//    //--
//    //debug : 行列・ベクトルチェック
//    //--
//    cout << "------ AssyMatrix::multVector 行列        ------" << endl;
//    dump();
//    cout << "------ AssyMatrix::multVector 左辺ベクトル ------" << endl;
//    pV->dump();//左辺ベクトル
//    cout << "------ AssyMatrix::multVector 右辺ベクトル ------" << endl;
//    pP->dump();//結果


    return 1;
}
//--
// solver iter 終了後にコール => MPCの最後の処理:  u = Tu' + ug
//--
uiint CAssyMatrix::multMPC(CAssyVector *pV, CAssyVector *pP) const
{
    if(mvMPCMatrix.size() == 0) {
        ;
    } else {
        for(uiint i = 0; i < mvMPCMatrix.size(); i++) {
            mvMPCMatrix[i]->multVector(pV, pP);
        };
    }

    return 1;
}

//
// 残差
//
uiint CAssyMatrix::residual(CAssyVector *pV, const CAssyVector *pF, CAssyVector *pR) const
{
    ////CHecMPI *pMPI= CHecMPI::Instance();

    multVector(pV, pR);//----- r=K*u

    for(uiint i = 0; i < pR->size(); i++) {
        (*pR)[i] = (*pF)[i] - (*pR)[i];//----- r = f - K*u
    };

    return 1;
}
//uiint CAssyMatrix::setupSolver(iint type){ return 1;}
uiint CAssyMatrix::setupSolver(iint iter_max, double tolerance, iint method, iint pre, bool b_iterlog,bool b_timelog)
{
    switch(method) {
    case(1):
        mpSolver = new CSolverCG(iter_max, tolerance, method, pre, b_iterlog, b_timelog);
        break;
    case(2):
        mpSolver = new CSolverBiCGSTAB(iter_max, tolerance, method, pre, b_iterlog, b_timelog);
        break;
    case(3):
        mpSolver = new CSolverGPBiCG(iter_max, tolerance, method, pre, b_iterlog, b_timelog);
        break;
    case(4):
        mpSolver = new CSolverGMRES(iter_max, tolerance, method, pre, b_iterlog, b_timelog);
        break;
    }
    return 1;
}
void CAssyMatrix::deleteSolver()
{
    if(mpSolver) delete mpSolver;
}
//
// preconditionの設定：コースグリッドマトリックスへもセットされていく
//
uiint CAssyMatrix::setupPreconditioner(iint type) const
{
    for (CVMatConstIter im = mvMatrix.begin(); im != mvMatrix.end(); im++) {
        (*im)->setupPreconditioner(type);
    }
    if( mpCoarseMatrix ) mpCoarseMatrix->setupPreconditioner(type);

    return 1;
}
uiint CAssyMatrix::setupSmoother(iint type)
{
    for(CVMatIter im = mvMatrix.begin(); im != mvMatrix.end(); im++) {
        (*im)->setupSmoother(type);
    }
    return 1;
}
//uiint CAssyMatrix::solve(const CAssyVector *pF, CAssyVector *pV, iint iter = 0) const
//{
//    int iter_save = mpSolver->getIterMax();
//    if (iter > 0) mpSolver->setIterMax(iter);
//    mpSolver->solve(this, pF, pV);
//    if (iter > 0) mpSolver->setIterMax(iter_save);
//    return 1;
//}

//
// 前処理MG smoother & solver:  pR:右辺、pZ:解
//
uiint CAssyMatrix::precond(const CAssyVector *pR, CAssyVector *pZ, iint iter) const
{
    CAssyVector z_resid(pR);

    ////////////////
    pZ->setZero();//------ 解ベクトル初期化:”0”クリア
    ////////////////

    for(iint i=0; i < iter; i++) {
        uiint index = 0;
        for(CVMatConstIter im = mvMatrix.begin(); im != mvMatrix.end(); im++) {
            (*im)->precond(pR->getVector(index), z_resid.getVector(index));
            index++;//--- パーツ番号
        }
        pZ->add(&z_resid);

        if (i == iter - 1) break;
        residual(pZ, pR, &z_resid);//---- 残差  Δu = R-K・u
    }
    return 1;
}
//
// MG_Solver本体のsmoother : pF:右辺、pV:解
//
uiint CAssyMatrix::relax(const CAssyVector *pF, CAssyVector *pV, iint iter) const
{
    for(iint i = 0; i < iter; i++) {
        uiint index = 0;
        for(CVMatConstIter im= mvMatrix.begin(); im != mvMatrix.end(); im++ ) {
            //---- SOR法 MatrixBCRS::relax ----
            (*im)->relax(pF->getVector(index), pV->getVector(index));
            index++;//--- パーツ番号
        };
    }
    return 1;
}
//--
// Restriction 行列 : coarse matrix 生成
//--
void CAssyMatrix::restrict()
{
    //-- mpCoarseMatrix --//
    if( mpCoarseMatrix == NULL ) return;
}

//--
// MG_Solver本体
//--
uiint CAssyMatrix::MGCycle(const CAssyVector *pRHSAssyVec, CAssyVector *pSolAssyVec, iint mg_cycle, iint alpha1, iint alpha2)// const
{
    //
    // アセンブル・ベクトルの各DOFは方程式内では共通なので解ベクトルからDOFを取得
    //
    uiint nNumOfMesh= pSolAssyVec->getNumOfMesh();
    vuint vDOF;
    vDOF.resize(nNumOfMesh);
    for(uiint imesh=0; imesh < nNumOfMesh; imesh++) {
        vDOF[imesh]= pSolAssyVec->getDOF(imesh);
    }

    if (mpCoarseMatrix != NULL) {
        //
        // restriction matrix
        //
        for(uiint i=0; i < mvMatrix.size(); i++){
            CMatrixBCRS *pCrsMat= mpCoarseMatrix->getMatrix(i);
            mvMatrix[i]->restrict(pCrsMat);//-------------------- コースグリッド行列
        };
        
        CAssyVector uf(mpAssyModel, vDOF);//--------------------- 修正量ファイングリッド
        CAssyVector rf(mpAssyModel, vDOF);//--------------------- 残差ファイングリッド

        CAssyVector uc(mpCoarseMatrix->mpAssyModel, vDOF);//----- 修正量コースグリッド
        CAssyVector rc(mpCoarseMatrix->mpAssyModel, vDOF);//----- 残差コースグリッド

        relax(pRHSAssyVec, pSolAssyVec, alpha1);//<<<<<<<<<<<<<<<<<<< pre_smoothing

        //---- 残差
        residual(pSolAssyVec, pRHSAssyVec, &rf);//---- ファインrf(残差)


        rf.restrictTo(&rc);//------- rfからrcにリストリクト:コースrc(残差)

        uc.setZero();//---uc = init_Δu: 修正量(コース)の初期値
        
        // -----------------------------------------------------
        // Kc*uc=rc : Kc:コース行列, uc:コース修正量, rc:コース残差
        // -----------------------------------------------------
        for(iint i=0; i < mg_cycle; i++) {
            mpCoarseMatrix->MGCycle(&rc, &uc, mg_cycle, alpha1, alpha2);//----- コースグリッド(修正方程式)
        }

        //---- 修正量
        uf.prolongateFrom(&uc);//---- ファイングリッドの修正量:コースグリッドの修正量をプロロンゲート


        pSolAssyVec->add(&uf);//----- ファイングリッドの残差の解(修正量)を加算
        
        relax(pRHSAssyVec, pSolAssyVec, alpha2);//<<<<<<<<<<<<<<<<<<<<<<<<< post_smoothing

    } else {
        //コースグリッド (Level=0) ソルバー
        if(mpSolver) mpSolver->solve(this, pRHSAssyVec, pSolAssyVec);
    }

    return 1;
}
//--
//  前処理MG : 各ソルバーのdoSolveから呼ばれる(MG_Solverとは処理が異なる)
//--
uiint CAssyMatrix::MGInitialGuess(const CAssyVector *pRHSAssyVec, CAssyVector *pSolAssyVec)// const
{
    // printf(" enter CAssyMatrix::MGInitialGuess %ld \n", mMGLevel);

    //
    // アセンブル・ベクトルの各DOFは方程式内では共通なので解ベクトルからDOFを取得
    //
    uiint nNumOfMesh= pSolAssyVec->getNumOfMesh();
    vuint vDOF;
    vDOF.resize(nNumOfMesh);
    for(uiint imesh=0; imesh < nNumOfMesh; imesh++) {
        vDOF[imesh]= pSolAssyVec->getDOF(imesh);
    }

    if (mpCoarseMatrix != NULL) {
        //////
        ////// restriction matrix
        //////
        ////for(uiint i=0; i < mvMatrix.size(); i++){
        ////    CMatrixBCRS *pCrsMat= mpCoarseMatrix->getMatrix(i);
        ////    mvMatrix[i]->restrict(pCrsMat);//-------------------- コースグリッド行列
        ////};

        CAssyVector uf(mpAssyModel, vDOF);//------------------- uf:修正量(ファイン)
        CAssyVector uc(mpCoarseMatrix->mpAssyModel, vDOF);//--- uc:修正量(コース)
        CAssyVector rf(mpAssyModel, vDOF);//------------------- rf:残差(ファイン)
        CAssyVector rc(mpCoarseMatrix->mpAssyModel, vDOF);//--- rc:残差(コース)

        residual(pSolAssyVec, pRHSAssyVec, &rf);//---------- 残差 rf: Δr = f - K・u

        rf.restrictTo(&rc);//---------- rfをrcに,リストリクト:残差(コース)
        uc.setZero();//---------------- 修正量(コース) 初期値

        mpCoarseMatrix->MGInitialGuess(&rc, &uc);//---------- コースグリッド修正量方程式:Kc*Δuc=rc

        uf.prolongateFrom(&uc);//------------ uc(修正量コース)をuf(修正量ファイン)にプロロンゲート

        pSolAssyVec->add(&uf);//--------------------- u = u + prolongate(Δuc)
        
        precond(pRHSAssyVec, pSolAssyVec, 1);//---- smoothing
        
    } else {
        //--
        // コースグリッド
        //--
        precond(pRHSAssyVec, pSolAssyVec, 1);
    }
    precond(pRHSAssyVec, pSolAssyVec, 1);//---- smoothing

    // printf(" exit CAssyMatrix::MGInitialGuess %ld \n", mMGLevel);
    return 1;
}
void CAssyMatrix::dump()
{
    uiint i;
    for(i=0; i < mvMatrix.size(); i++) {
        mvMatrix[i]->dump();
    };
    for(i=0; i < mvMPCMatrix.size(); i++) {
        mvMPCMatrix[i]->dump();
    };
}
void CAssyMatrix::dump() const
{
    uiint i;
    for(i=0; i < mvMatrix.size(); i++) {
        mvMatrix[i]->dump();
    };
    for(i=0; i < mvMPCMatrix.size(); i++) {
        mvMPCMatrix[i]->dump();
    };
}
}//namespace pmw;
