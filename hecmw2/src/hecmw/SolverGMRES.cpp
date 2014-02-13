/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/SolverGMRES.cpp
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
#include "SolverGMRES.h"
#include "AssyMatrix.h"
#include "AssyVector.h"
#include <cmath>

namespace pmw
{
CSolverGMRES::CSolverGMRES(iint iter_max = 100,
                           double tolerance = 1.0e-8,
                           iint method = 1,
                           iint precondition = 1,
                           bool flag_iter_log = false,
                           bool flag_time_log = false)
    : CSolver(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log)
{
}
CSolverGMRES::~CSolverGMRES()
{
}
//uiint CSolverGMRES::doSolve_(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
//		iint iter_max, double tolerance,
//		bool flag_iter_log, bool flag_time_log)
//{
//
//}
uiint CSolverGMRES::doSolve( CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
                             iint iter_max, double tolerance,
                             bool flag_iter_log, bool flag_time_log)
{
    double bnrm2_inv, rnrm2;
    double alpha, beta, rho, rho_prev, resid;
    iint iter, iter_precond;
#ifdef ADVANCESOFT_DEBUG
    printf(" enter CSolverGMRES::doSolve \n");
#endif
    printf(" --- start of GMRES solver --- \n");
    uiint len = pB->size();
    CAssyVector R(pX);
    CAssyVector Z(pX);
    CAssyVector P(pX);
    CAssyVector &Q = Z;
    pA->residual(pX, pB, &R);
    bnrm2_inv = 1.0 / pB->norm2();
    iint itype = getPrecondition();
    pA->setupPreconditioner( itype );
    Z.setZero();
    for (iter = 0; iter < getIterMax(); iter++) {
        iter_precond = 1;
        pA->precond(&R, &Z, iter_precond);
        rho = R.innerProd(&Z);
        if (iter > 0) {
            beta = rho / rho_prev;
            for (int i = 0; i < len; i++) {
                P[i] = Z[i] + beta * P[i];
            }
        } else {
            P.subst(&Z);
        }
        pA->multVector(&P, &Q);
        alpha = rho / P.innerProd(&Q);
        for (uiint i = 0; i < len; i++) {
            (*pX)[i] += alpha * P[i];
            R[i]     -= alpha * Q[i];
        }
        rnrm2 = R.norm2();
        resid = sqrt(rnrm2 * bnrm2_inv);
        if (getFlagIterLog()) {
            printf("%5ld %16.6e\n", iter + 1, resid);
        }
        printf("iteration:%5ld, residue:%16.6e \n", iter + 1, resid);
        if (resid < getTolerance()) break;
        rho_prev = rho;
    }
    if( iter == getIterMax() ) return 0;
    P.subst(pX);
    pA->multMPC(&P, pX);
    printf(" --- end of GMRES solver --- \n");
    return 1;
}
}
