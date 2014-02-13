/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/MPCMatrix.cpp
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
#include "MPCMatrix.h"
#include "AssyVector.h"
#include "Equation.h"

#include "Logger.h"

namespace pmw
{
CMPCMatrix::CMPCMatrix():mnEquation(0)
{
    ;
}
CMPCMatrix::~CMPCMatrix()
{
    for_each(mvEquation.begin(), mvEquation.end(), DeleteObject());
}

typedef std::vector<CEquation*> CVEqn;
typedef CVEqn::iterator CVEqnIter;
typedef CVEqn::const_iterator CVEqnConstIter;
typedef std::vector<CEqnTerm> CVEqTrm;
typedef CVEqTrm::iterator CVEqTrmIter;
typedef CVEqTrm::const_iterator CVEqTrmConstIter;

//--
// Tu: Coef*pV = pP
//--
void CMPCMatrix::multVector(CAssyVector *pV, CAssyVector *pP) const
{
    CHecMPI *pMPI= CHecMPI::Instance();

    ////    //debug file.
    ////    string mfile("mpc_multvec_");
    ////    stringstream ss; ss << pMPI->getRank();
    ////    mfile += ss.str();
    ////    mfile += ".txt";
    ////    ofstream ofs(mfile.c_str(), ios::out | ios::app);


    pP->subst(pV);//-------------------- pV値をpPへ代入

    //--
    // MPC_CoefEquation数 == スレーブ点*DOF数
    // ## スレーブ点一個につき、DOF数のEquationが並んでいる.
    //--
    for(CVEqnConstIter i_eqn = mvEquation.begin(); i_eqn != mvEquation.end(); i_eqn++) {

        const CVEqTrm &vTerm = (*i_eqn)->getTerm();
        CVEqTrmConstIter term= vTerm.begin();//---- 先頭Term:スレーブ点
        uiint nSlaveCommRank= term->comm_rank();// スレーブ点Rank(通信)
        uiint nSlaveSolvRank= term->solv_rank();// スレーブ点Rnak(計算)

        uiint nSlaveMeshIX= term->meshIX(), nSlaveNodeIX= term->nodeIX(), nSlaveDOF= term->dof();


        double sum(0.0), val(0.0), val_c(0.0), val_v(0.0);
        term++;
        //--
        // 次からのTerm : マスター・ベクトルとの演算
        //--
        for( ; term != vTerm.end(); term++) {
            //--
            // 行列ベクトル積の前にupdateされているので、値は同一:comm_rankの値を利用
            //--
            // マスター "pVベクトル*Coef" 演算
            if( term->comm_rank()==pMPI->getRank() ) {
                val = term->coef() * (*pV)(term->meshIX(), term->nodeIX(), term->dof());
            }
            pMPI->Bcast((void*)&val, 1, MPI_DOUBLE, term->comm_rank(), MPI_COMM_WORLD);

            sum += val;//------------------------------------------- sum

        };//--- vTerm ループ End

        ////cout << "MPCMatrix::multVector ----------- nSlaveCommRank:" << nSlaveCommRank << "  rank:" << pMPI->getRank() << endl;

        //--
        // Pベクトル・スレーブ点 に代入
        // # 最小rankに加算値を代入
        //--
        if( nSlaveCommRank==pMPI->getRank() ) {

            ////cout << "MPCMatrix::multVector ----------- sum:" << sum
            ////        << " nSlaveMeshIX:" << nSlaveMeshIX << " nSlaveNodeIX:" << nSlaveNodeIX << " nSlaveDOF:" << nSlaveDOF
            ////        << "  rank:" << pMPI->getRank() << endl;

            (*pP)( nSlaveMeshIX, nSlaveNodeIX, nSlaveDOF ) = sum;// スレーブpPベクトルに"pVベクトル*Coef"を加算 : sumup
        } else {
            // スレーブ点が相手Rank
            // # 加算すべきベクトルは無い -> なにもしない.
        }
    };
}


//--
// T転置 行列の積 : pPベクトルへのSumup   (pVとpPの関係 : Av=p)
//--
void CMPCMatrix::transMultVector(CAssyVector *pV, CAssyVector *pP) const
{
    CHecMPI *pMPI= CHecMPI::Instance();

    pP->subst(pV);//----------- pV値をpPへ代入

    //--
    // MPC_CoefEquation数 == スレーブ点*DOF数
    // ## スレーブ点一個につき、DOF数のEquationが並んでいる.
    //--
    for(CVEqnConstIter i_eqn = mvEquation.begin(); i_eqn != mvEquation.end(); i_eqn++) {

        const CVEqTrm &vTerm = (*i_eqn)->getTerm();
        CVEqTrmConstIter term= vTerm.begin();       // 先頭Term : スレーブ点
        uiint nSlaveCommRank= term->comm_rank();    // スレーブ点Rank(通信)
        uiint nSlaveSolvRank= term->solv_rank();    // スレーブ点Rank(計算)

        //--
        // ベクトル"P"のスレーブ点
        //--
        if(nSlaveSolvRank==pMPI->getRank()) {
            (*pP)(term->meshIX(), term->nodeIX(), term->dof()) = 0.0;
        } else {
            // スレーブ点が相手Rank
            // # 操作するベクトルが無い -> なにもしない.
        }
        //--
        // ベクトル"V"のスレーブ点の値
        //--
        double slave_val(0.0);
        if( nSlaveCommRank==pMPI->getRank() ) {
            slave_val = (*pV)(term->meshIX(), term->nodeIX(), term->dof());
        }
        pMPI->Bcast((void*)&slave_val, 1, MPI_DOUBLE, nSlaveCommRank, MPI_COMM_WORLD);


        ////cout << "MPCMatrix::transMultVector ----------- nSlaveCommRank:" << nSlaveCommRank << "  rank:" << pMPI->getRank() << endl;


        term++;// 次からのTerm: マスター点
        //--
        //  次からのTerm : マスター面の頂点のCoef
        //--
        for( ; term != vTerm.end(); term++) {
            //--
            // マスター pPベクトルに "slave_val * Coef" を加算 : sumup
            // # 最小rankに加算
            //--
            if(term->comm_rank()==pMPI->getRank()) {
                (*pP)(term->meshIX(), term->nodeIX(), term->dof()) += term->coef() * slave_val;
            } else {
                // マスター頂点が相手Rank
                // # 操作対象となるベクトルが無い -> なにもしない.
            }
        };
    };
}


//--
// dump
//--
void CMPCMatrix::dump()
{
    uiint nNumOfEqu = mvEquation.size();
    uiint i,ii;
    for(i=0; i < nNumOfEqu; i++) {
        uiint nNumOfTerm;
        nNumOfTerm = mvEquation[i]->numTerm();
        cout << " " ;
        for(ii=0; ii < nNumOfTerm; ii++) {
            cout << " NodeID=" << mvEquation[i]->getTerm(ii).nodeID();
            cout << " DOF=" << mvEquation[i]->getTerm(ii).dof();
            cout << " Coef=" << mvEquation[i]->getTerm(ii).coef();
        };
        cout << endl;
    };
}
}
