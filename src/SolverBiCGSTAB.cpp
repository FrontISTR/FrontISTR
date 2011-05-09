/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   SolverBiCGSTAB.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
/*
 * SolverBiCGSTAB.cpp
 *
 *  Created on: Jul 24, 2009
 *      Author: goto
 */
#include "SolverBiCGSTAB.h"
#include "AssyMatrix.h"
#include "AssyVector.h"
#include <cmath>
#include <cstdio>
namespace pmw
{
CSolverBiCGSTAB::CSolverBiCGSTAB(int iter_max = 100,
		double tolerance = 1.0e-8,
		uint method = 1,
		uint precondition = 1,
		bool flag_iter_log = false,
		bool flag_time_log = false)
	: CSolver(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log)
{
}
CSolverBiCGSTAB::~CSolverBiCGSTAB()
{
}
int CSolverBiCGSTAB::doSolve(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
		int iter_max, double tolerance,
		bool flag_iter_log, bool flag_time_log)
{
	double bnrm2_inv, snrm2;
	double alpha, beta, omega, rho, rho_prev, resid;
	int iter_precond; 
#ifdef ADVANCESOFT_DEBUG
   	printf(" enter CSolverBiCGSTAB::doSolve1111 \n");
#endif
   	printf(" --- start of BiCGSTAB solver --- \n");
	CAssyVector R (pX); 
	CAssyVector Rt(pX); 
	CAssyVector P (pX); 
	CAssyVector Pt(pX); 
	CAssyVector S (pX);
	CAssyVector &St= R;  
	CAssyVector T (pX); 
	CAssyVector V (pX); 
	int len = pB->size();
	pA->multVector(pX, &P);
	for (int i = 0; i < len; i++) {
		R[i] = (*pB)[i] - P[i];
		Rt[i] = R[i];
	}
	bnrm2_inv = 1.0 / pB->norm2();
	int itype = getPrecondition();
	pA->setupPreconditioner( itype ); 
	for (int iter = 0; iter < getIterMax(); iter++) {
		rho = Rt.innerProd(&R);
		if (iter > 0) {
			beta = (rho * alpha) / (rho_prev * omega);
			for (int i = 0; i < len; i++) {
				P[i] = R[i] + beta * (P[i] - omega * V[i]);
			}
		} else {
			for (int i = 0; i < len; i++) {
				P[i] = R[i];
			}
		}
		iter_precond = 1;
		if( itype == 2 ) {
			Pt.subst(pX);
			pA->MGInitialGuess(&P, &Pt);
		} else {
		  pA->precond(&P, &Pt, iter_precond);
		}
		pA->multVector(&Pt, &V);
		alpha = rho / Rt.innerProd(&V);
		for (int i = 0; i < len; i++) {
			S[i] = R[i] - alpha * V[i];
		}
		iter_precond = 1;
		if( itype == 2 ) {
			St.subst(pX);
			pA->MGInitialGuess(&S, &St);
		} else {
		  pA->precond(&S, &St, iter_precond);
		}
		pA->multVector(&St, &T);
		omega = T.innerProd(&S) / T.norm2();
		for (int i = 0; i < len; i++) {
			(*pX)[i] += alpha * Pt[i] + omega * St[i];
			R[i] = S[i] - omega * T[i];
		}
		snrm2 = S.norm2();
		resid = sqrt(snrm2 * bnrm2_inv);
		if (getFlagIterLog()) {
			printf("%5d %16.6e\n", iter + 1, resid);
		}
		printf("iteration:%5d, residue:%16.6e \n", iter + 1, resid);
		if (resid < getTolerance()) break;
		rho_prev = rho;
	}
	P.subst(pX);
	pA->multMPC(&P, pX);
   	printf(" --- end of BiCGSTAB solver --- \n");
	return 1;
}
}
