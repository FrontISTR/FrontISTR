/*
 * SolverCG.cpp
 *
 *  Created on: Jul 23, 2009
 *      Author: goto
 */

#include "SolverCG.h"
#include "AssyMatrix.h"
#include "AssyVector.h"
#include <cmath>
#include <cstdio>

namespace pmw
{

CSolverCG::CSolverCG(iint iter_max = 100,
		double tolerance = 1.0e-8,
		iint method = 1,
		iint precondition = 1,
		bool flag_iter_log = false,
		bool flag_time_log = false)
	: CSolver(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log)
{
	// TODO Auto-generated constructor stub
}

CSolverCG::~CSolverCG()
{
	// TODO Auto-generated destructor stub
}

//
// Ax=b :  pA, pX, pB
//
uiint CSolverCG::doSolve(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
		iint iter_max, double tolerance,
		bool flag_iter_log, bool flag_time_log)
{
	double bnrm2_inv, rnrm2;
	double alpha, beta, rho, rho_prev, resid;
	iint iter, iter_precond; // TODO: class variable??

#ifdef ADVANCESOFT_DEBUG
   	printf(" enter CSolverCG::doSolve1111 \n");
#endif
   	printf(" --- start of CG solver --- %ld \n",getPrecondition());

	uiint len = pB->size();

	CAssyVector R(pX); // same type as pB and pX
	CAssyVector Z(pX); //
	CAssyVector P(pX); //
	CAssyVector &Q = Z;

	// {R} = {B} - [A] {X}
	pA->residual(pX, pB, &R);
	
	// calc 1 / |{B}|^2
	bnrm2_inv = 1.0 / pB->norm2();

	iint itype = getPrecondition();// 前処理
	pA->setupPreconditioner( itype ); // TODO: rewrite to pA->precond(pR, pZ);
	Z.setZero(); //K.Matsubara

	for (iter = 0; iter < getIterMax(); iter++) {

		// {Z} = [Minv] {R}
		iter_precond = 1;

                // 前処理選択
		if( itype == 2 ) {
			Z.subst(pX);
			pA->MGInitialGuess(&R, &Z);
		} else {
			pA->precond(&R, &Z, iter_precond); // TODO: rewrite to pA->precond(pR, pZ);
		}

		// rho = {R} {Z}
		rho = R.innerProd(&Z);


		// iter > 0:
		//   beta = rho / rho_prev
		//   {P} = {Z} + beta {P}
		// iter == 1:
		//   {P} = {Z}
		if (iter > 0) {
                    beta = rho / rho_prev;

                    // Z.sumSV(beta, &P, &P);
                    for (uiint i = 0; i < len; i++) {
                            P[i] = Z[i] + beta * P[i];
                    }
		} else {
			P.subst(&Z);
		}

		// {Q} = [A] {P}
                //-------------------------------------------------------------------------//
                //		P.updateCommBoundary();// parallel processing : multVector内
		pA->multVector(&P, &Q);
                //		Q.sumupCommBoundary(); // parallel processing : multVector内
                //-------------------------------------------------------------------------//

                
		// alpha = rho / ({P} {Q})
		alpha = rho / P.innerProd(&Q);
					
		// {X} = {X} + alpha {P}
		// {R} = {R} - alpha {Q}
		for (uiint i = 0; i < len; i++) {
                    (*pX)[i] += alpha * P[i];
                    R[i]     -= alpha * Q[i];
		}

		// calc |{R}|^2
		rnrm2 = R.norm2();

		// resid = |{R}| / |{B}|
		resid = sqrt(rnrm2 * bnrm2_inv);

		// iteration history
		if (getFlagIterLog()) {
                    // TODO: replace with logger
                    printf("%5ld %16.6e\n", iter + 1, resid);
		}
		printf("iteration:%5ld, residue:%16.6e \n", iter + 1, resid);

		// check convergence
		if (resid < getTolerance()) break;

		// rho_prev = rho
		rho_prev = rho;

#ifdef ADVANCESOFT_DEBUG
//	printf(" ----- the solution vector in loop ----- \n");
//	for(int i=0; i < len; i++)
//		printf(" %d : %e %e %e \n", i, (*pX)[i][0], (*pX)[i][1], (*pX)[i][2]);
#endif

	}
	
	if( iter == getIterMax() ) return 0;

	// K.Matsubara 2010.03.31
	P.subst(pX);
	pA->multMPC(&P, pX);

	// interface data exchange
        //	pX->updateCommBoundary(); // pA->multVector内

   	printf(" --- end of CG solver --- \n");
	
#ifdef ADVANCESOFT_DEBUG
	printf(" ----- the solution vector of solver ----- \n");
	for(uiint i=0; i < 5; i++)
		printf(" %d : %e %e %e \n", i, (*pX)[i][0], (*pX)[i][1], (*pX)[i][2]);
	for(uiint i=len-5; i < len; i++)
//	for(uint i=0; i < len; i++)
		printf(" %d : %e %e %e \n", i, (*pX)[i][0], (*pX)[i][1], (*pX)[i][2]);
   	printf(" exit CSolverCG::doSolve \n");
#endif

	return 1;
}

}
