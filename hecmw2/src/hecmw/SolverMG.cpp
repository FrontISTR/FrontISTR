/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/SolverMG.cpp
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
#include "SolverMG.h"
#include "AssyMatrix.h"
#include "AssyVector.h"
namespace pmw
{
CSolverMG::CSolverMG(iint iter_max = 100,
                     double tolerance = 1.0e-8,
                     iint method = 1,
                     iint precondition = 1,
                     bool flag_iter_log = false,
                     bool flag_time_log = false)
    : CSolver(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log)
{
}
CSolverMG::~CSolverMG()
{
}
//uiint CSolverMG::doSolve_(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
//		iint iter_max, double tolerance,
//		bool flag_iter_log, bool flag_time_log)
//{
//
//}
uiint CSolverMG::doSolve( CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
                          iint iter_max, double tolerance, bool flag_iter_log, bool flag_time_log)
{

    setCGridSolver( pA, iter_max, tolerance, flag_iter_log, flag_time_log);//--- コースグリッド線形ソルバー生成・セット

    pA->setupPreconditioner( mPrecondition );//前処理種類の設定:CoarseMatrixへもセットされる

    double bnrm2_inv, rnrm2;
    double resid;
    iint iter;

    CAssyVector R(pX);
    bnrm2_inv = 1.0 / pB->norm2();

    for( iter=0; iter < mIterMax; iter++) {

        iint mg_cycle = 1;// 1:V-Cycle
        iint alpha1(10);//sor iter :pre smooth
        iint alpha2(10);//sor iter :post smooth

        pA->MGCycle(pB, pX, mg_cycle, alpha1, alpha2);

        pA->residual(pX, pB, &R);

        rnrm2 = R.norm2();
        resid = sqrt(rnrm2 * bnrm2_inv);

        cout << " mg_cycle:" << iter+1 << "  tolerance:" << resid << endl;  // << "  R.norm2:" << rnrm2 << endl;

        if( getFlagIterLog() ) {
            ;
        }
        if( resid < getTolerance() ) break;

    };//iter ループ

    if( iter == getIterMax() ) {
        deleteCGridSolver();//--- コースグリッドソルバー削除
        return 0;
    }

    ////pX->updateCommBoundary();//---- 通信界面 update

    deleteCGridSolver();//--- コースグリッドソルバー削除

    return 1;
}
}
