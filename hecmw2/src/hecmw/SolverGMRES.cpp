/*
 * SolverGMRES.cpp
 *
 *  Created on: Jul 24, 2009
 *      Author: goto
 */

#include "SolverGMRES.h"
#include "AssyMatrix.h"
#include "AssyVector.h"
#include <cmath>
#include <cstdio>

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
	// TODO Auto-generated constructor stub

}

CSolverGMRES::~CSolverGMRES()
{
	// TODO Auto-generated destructor stub
}

uiint CSolverGMRES::doSolve(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
		iint iter_max, double tolerance,
		bool flag_iter_log, bool flag_time_log)
{
	double bnrm2_inv, rnrm2;
	double alpha, beta, rho, rho_prev, resid;
	iint iter, iter_precond; // TODO: class variable??

#ifdef ADVANCESOFT_DEBUG
   	printf(" enter CSolverGMRES::doSolve \n");
#endif
   	printf(" --- start of GMRES solver --- \n");

	uiint len = pB->size();

	CAssyVector R(pX); // same type as pB and pX
	CAssyVector Z(pX); //
	CAssyVector P(pX); //
	CAssyVector &Q = Z;

	// {R} = {B} - [A] {X}
	pA->residual(pX, pB, &R);
	
	// calc 1 / |{B}|^2
	bnrm2_inv = 1.0 / pB->norm2();

	iint itype = getPrecondition();
	pA->setupPreconditioner( itype ); // TODO: rewrite to pA->precond(pR, pZ);
	Z.setZero(); //K.Matsubara

	for (iter = 0; iter < getIterMax(); iter++) {

		// {Z} = [Minv] {R}
		iter_precond = 1;
		pA->precond(&R, &Z, iter_precond); // TODO: rewrite to pA->precond(pR, pZ);
//		Z.subst(pX);
//		pA->MGInitialGuess(&R, &Z);

		// rho = {R} {Z}
		rho = R.innerProd(&Z);

		if (iter > 0) {
			beta = rho / rho_prev;
			for (int i = 0; i < len; i++) {
				P[i] = Z[i] + beta * P[i];
			}
		} else {
			P.subst(&Z);
		}

		// {Q} = [A] {P}
		pA->multVector(&P, &Q);

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

	}
	
	if( iter == getIterMax() ) return 0;

	// K.Matsubara 2010.03.31
	P.subst(pX);
	pA->multMPC(&P, pX);

   	printf(" --- end of GMRES solver --- \n");
	return 1;
}

}
