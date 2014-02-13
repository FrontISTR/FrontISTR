/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Solver.cpp
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
#include "Solver.h"
#include "AssyMatrix.h"
#include "AssyVector.h"
namespace pmw
{
CSolver::CSolver(iint iter_max = 100, double tolerance = 1.0e-8, iint method = 1, iint precondition = 1,
                 bool flag_iter_log = false, bool flag_time_log = false)
    : mIterMax(iter_max), mTolerance(tolerance), mMethod(method), mPrecondition(precondition),
      mFlagIterLog(flag_iter_log), mFlagTimeLog(flag_time_log)
{
    mFlagNewMesh = true;
    mFlagNewCoef = true;

    mpCrsMatA = NULL;
}
CSolver::~CSolver()
{
}
void CSolver::setFlagNewMesh(bool flag_new)
{
    mFlagNewMesh = flag_new;
}
void CSolver::setFlagNewCoef(bool flag_new)
{
    mFlagNewCoef = flag_new;
}
void CSolver::setIterMax(iint iter_new)
{
    mIterMax = iter_new;
}
iint CSolver::getIterMax()
{
    return mIterMax;
}
void CSolver::setTolerance(double tol_new)
{
    mTolerance = tol_new;
}
double CSolver::getTolerance()
{
    return mTolerance;
}
void CSolver::setMethod(uiint method)
{
    mMethod = method;
}
uiint CSolver::getMethod()
{
    return mMethod;
}
void CSolver::setPrecondition(iint precondition)
{
    mPrecondition = precondition;
}
iint CSolver::getPrecondition()
{
    return mPrecondition;
}
bool CSolver::getFlagIterLog()
{
    return mFlagIterLog;
}
void CSolver::setFlagIterLog(bool flag_new)
{
    mFlagIterLog = true;
}
void CSolver::setFlagTimeLog(bool flag_new)
{
    mFlagTimeLog = true;
}
bool CSolver::getFlagTimeLog()
{
    return mFlagTimeLog;
}
//uiint CSolver::solve_(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX)
//{
//}
uiint CSolver::solve(CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX)
{
    if (mFlagNewMesh) {
        mFlagNewMesh = false;
        mFlagNewCoef = false;
    } else if (mFlagNewCoef) {
        mFlagNewCoef = false;
    }
    ////CHecMPI *pMPI= CHecMPI::Instance();
    ////pMPI->Barrier(MPI_COMM_WORLD);
    ////cout << "CSolver::solve Start------------------------ PE:" << pMPI->getRank() << endl;

    doSolve(pA, pB, pX, mIterMax, mTolerance, mFlagIterLog, mFlagTimeLog);

    ////pMPI->Barrier(MPI_COMM_WORLD);
    ////cout << "CSolver::solve End-------------------------- PE:" << pMPI->getRank() << endl;

    return 1;
}
// コースグリッド・ソルバーのセット: コースマトリックスに線形ソルバーを生成
uiint CSolver::setCGridSolver(CAssyMatrix* pA, iint iter_max, double tolerance, bool flag_iter_log, bool flag_time_log)
{
    iint method;
    if(mMethod > 4 ) {
        method = mMethod-4;// 5-7:MGSolver == 5:CG-MG, 6:BiCGSTAB-MG, 7:GPBiCG-MG
    } else {
        method = mMethod;  // 1-4:MG前処理  == 1:CG, 2:BiCBSTAB, 3:GPBiCG, 4:GMRES
    }

    //----'11.08----//コースグリッドのMatrixにSolverセット
    CAssyMatrix *pCrsA_curr= pA->getCoarseMatrix();
    CAssyMatrix *pCrsA_next;

    if(pCrsA_curr) { //コースグリッドのpAの場合にスルー
        while(pCrsA_curr) {
            pCrsA_next = pCrsA_curr->getCoarseMatrix();

            if(pCrsA_next) {
                pCrsA_curr=pCrsA_next;
            } else {
                break;
            }
        };
        pCrsA_curr->setupSolver(iter_max, tolerance, method, 3, false, false);//引数:3 == 前処理SSOR
    }
    mpCrsMatA = pCrsA_curr;//delete処理のため保存

    return 1;
}
// コースマトリックスの線形ソルバーを削除
uiint CSolver::deleteCGridSolver()
{
    if(mpCrsMatA) mpCrsMatA->deleteSolver();//コースグリッドソルバー削除

    return 0;
}
}
