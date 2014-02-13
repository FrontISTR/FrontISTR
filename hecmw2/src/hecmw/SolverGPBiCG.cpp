/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/SolverGPBiCG.cpp
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
#include "SolverGPBiCG.h"
#include "AssyMatrix.h"
#include "AssyVector.h"
#include <cmath>

namespace pmw
{
CSolverGPBiCG::CSolverGPBiCG(iint iter_max = 100,
                             double tolerance = 1.0e-8,
                             iint method = 1,
                             iint precondition = 1,
                             bool flag_iter_log = false,
                             bool flag_time_log = false)
    : CSolver(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log)
{
}
CSolverGPBiCG::~CSolverGPBiCG()
{
}
//uiint CSolverGPBiCG::doSolve_(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
//		iint iter_max, double tolerance,
//		bool flag_iter_log, bool flag_time_log)
//{
//
//}
uiint CSolverGPBiCG::doSolve( CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
                              iint iter_max, double tolerance,
                              bool flag_iter_log, bool flag_time_log)
{
    double bnrm2_inv;
    double alpha, beta, rho, rho_prev, resid;
    iint iter, iter_precond;
    printf(" --- start of GPBiCG solver --- \n");
    CAssyVector R (pX);
    CAssyVector Rt(pX);
    CAssyVector Rp(pX);
    CAssyVector T (pX);
    CAssyVector Tt(pX);
    CAssyVector T0(pX);
    CAssyVector P (pX);
    CAssyVector Pt(pX);
    CAssyVector U (pX);
    CAssyVector W1(pX);
    CAssyVector Y (pX);
    CAssyVector Z (pX);
    CAssyVector WK(pX);
    CAssyVector W2(pX);
    uiint len = pB->size();
    pA->multVector(pX, &P);
    for (int i = 0; i < len; i++) {
        R[i] = (*pB)[i] - P[i];
        Rt[i] = R[i];
    }
    bnrm2_inv = 1.0 / pB->norm2();
    rho = Rt.innerProd(&R);
    iint itype = getPrecondition();
    pA->setupPreconditioner( itype );

    setCGridSolver( pA, iter_max, tolerance, flag_iter_log, flag_time_log);//コースグリッド線形ソルバー生成・セット

    for (iter = 0; iter < getIterMax(); iter++) {
        iter_precond = 1;
        if( itype == 2 ) {
            Rp.subst(pX);
            pA->MGInitialGuess(&R, &Rp);//--------- MultiGrid前処理
        } else {
            pA->precond(&R, &Rp, iter_precond);
        }
        if (iter > 0) {
            for (uiint i = 0; i < len; i++) {
                P[i] = Rp[i] + beta * (P[i] - U[i]);
            }
        } else {
            for (uiint i = 0; i < len; i++) {
                P[i] = Rp[i];
            }
        }
        pA->multVector(&P, &Pt);
        alpha = rho / Rt.innerProd(&Pt);
        for (uiint i = 0; i < len; i++) {
            Y[i] = T[i] - R[i] + alpha * (- W1[i] + Pt[i]);
            T[i] = R[i] - alpha * Pt[i];
        }
        iter_precond = 1;
        if( itype == 2 ) {
            Tt.subst(pX);
            pA->MGInitialGuess(&T, &Tt);// ------- MultiGrid前処理
        } else {
            pA->precond(&T, &Tt, iter_precond);
        }
        iter_precond = 1;
        if( itype == 2 ) {
            W2.subst(pX);
            pA->MGInitialGuess(&T0, &W2);// ------- MultiGrid前処理
        } else {
            pA->precond(&T0, &W2, iter_precond);
        }
        for (int i = 0; i < len; i++) {
            T0[i] = W2[i];
        }
        iter_precond = 1;
        if( itype == 2 ) {
            W2.subst(pX);
            pA->MGInitialGuess(&Pt, &W2);// ------- MultiGrid前処理
        } else {
            pA->precond(&Pt, &W2, iter_precond);
        }
        pA->multVector(&Tt, &WK);
        for (uiint i = 0; i < len; i++) {
            Tt[i] = WK[i];
        }
        double c[5];
        c[0] = Y.norm2();
        c[1] = Tt.innerProd(&T);
        c[2] = Y.innerProd(&T);
        c[3] = Tt.innerProd(&Y);
        c[4] = Tt.norm2();
        double qsi, eta;
        if (iter > 0) {
            double denom = c[4] * c[0] - c[3] * c[3];
            qsi = (c[0] * c[1] - c[2] * c[3]) / denom;
            eta = (c[4] * c[2] - c[3] * c[1]) / denom;
        } else {
            qsi = c[1] / c[4];
            eta = 0.0;
        }
        for (uiint i = 0; i < len; i++) {
            U[i] = qsi * W2[i] + eta * (T0[i] - Rp[i] + beta * U[i]);
            Z[i] = qsi * Rp[i] + eta * Z[i] - alpha * U[i];
        }
        for (uiint i = 0; i < len; i++) {
            (*pX)[i] += alpha * P[i] + Z[i];
            R[i] = T[i] - eta * Y[i] - qsi * Tt[i];
            T0[i] = T[i];
        }
        double rnrm2 = R.norm2();
        rho_prev = rho;
        rho = Rt.innerProd(&R);
        beta = alpha * rho / (qsi * rho_prev);
        for (uiint i = 0; i < len; i++) {
            W1[i] = Tt[i] + beta * Pt[i];
        }
        resid = sqrt(rnrm2 * bnrm2_inv);
        if (getFlagIterLog()) {
            printf("%5ld %16.6e\n", iter + 1, resid);
        }
        printf("iteration:%5ld, residue:%16.6e \n", iter + 1, resid);
        if (resid < getTolerance()) break;
    }//iter ループ

    if( iter == getIterMax() ) {
        deleteCGridSolver();//コースグリッドソルバー削除
        return 0;
    }

    P.subst(pX);
    pA->multMPC(&P, pX);

    printf(" --- end of GPBiCG solver --- \n");

    deleteCGridSolver();//コースグリッドソルバー削除

    return 1;
}
}
