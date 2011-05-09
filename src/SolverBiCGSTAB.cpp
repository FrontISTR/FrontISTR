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
	// TODO Auto-generated constructor stub
}

CSolverBiCGSTAB::~CSolverBiCGSTAB()
{
	// TODO Auto-generated destructor stub
}

int CSolverBiCGSTAB::doSolve(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
		int iter_max, double tolerance,
		bool flag_iter_log, bool flag_time_log)
{
	double bnrm2_inv, snrm2;
	double alpha, beta, omega, rho, rho_prev, resid;
	int iter_precond; // TODO: class variable???

#ifdef ADVANCESOFT_DEBUG
   	printf(" enter CSolverBiCGSTAB::doSolve1111 \n");
#endif
   	printf(" --- start of BiCGSTAB solver --- \n");

	CAssyVector R (pX); // same type as pB and pX
	CAssyVector Rt(pX); //
	CAssyVector P (pX); //
	CAssyVector Pt(pX); //
	CAssyVector S (pX);
	CAssyVector &St= R;  // can be shared with R
	CAssyVector T (pX); //
	CAssyVector V (pX); //

	int len = pB->size();

	// {R} = {B} - [A] {X}
	// {Rt} = {R}
	pA->multVector(pX, &P);
	for (int i = 0; i < len; i++) {
		R[i] = (*pB)[i] - P[i];
		Rt[i] = R[i];
	}

	// calc 1 / |{B}|^2
	bnrm2_inv = 1.0 / pB->norm2();

	int itype = getPrecondition();
	pA->setupPreconditioner( itype ); // TODO: rewrite to pA->precond(pR, pZ);

	for (int iter = 0; iter < getIterMax(); iter++) {

		// rho = {Rt} {R}
		rho = Rt.innerProd(&R);

		// iter > 0:
		//   beta = (rho / rho_prev) * (alpha / omega)
		//   {P} = {R} + beta ( {P} - omega {V} )
		// iter == 1:
		//   {P} = {R}
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

		// {Pt} = [Minv] {P}
		iter_precond = 1;
		if( itype == 2 ) {
			Pt.subst(pX);
			pA->MGInitialGuess(&P, &Pt);
		} else {
		  pA->precond(&P, &Pt, iter_precond);
		}

		// {V} = [A] {Pt}
		pA->multVector(&Pt, &V);

		// alpha = rho / ({Rt} {V})
		alpha = rho / Rt.innerProd(&V);

		// {S} = {R} - alpha {V}
		for (int i = 0; i < len; i++) {
			S[i] = R[i] - alpha * V[i];
		}

		// {St} = [Minv] {S}
		iter_precond = 1;

		if( itype == 2 ) {
			St.subst(pX);
			pA->MGInitialGuess(&S, &St);
		} else {
		  pA->precond(&S, &St, iter_precond);
		}
		// {T} = [A] {St}
		pA->multVector(&St, &T);

		// omega = ({T} {S}) / ({T} {T})
		omega = T.innerProd(&S) / T.norm2();
		//printf(" - omega %e \n", omega);

		// {X} = {X} + alpha {Pt} + omega {St}
		// {R} = {S} - omega {T}
		for (int i = 0; i < len; i++) {
			(*pX)[i] += alpha * Pt[i] + omega * St[i];
			R[i] = S[i] - omega * T[i];
		}

		// calc |{S}|^2
		snrm2 = S.norm2();

		// resid = |{S}| / |{B}|
		resid = sqrt(snrm2 * bnrm2_inv);

		// iteration history
		if (getFlagIterLog()) {
			// TODO: replace with logger
			printf("%5d %16.6e\n", iter + 1, resid);
		}
		printf("iteration:%5d, residue:%16.6e \n", iter + 1, resid);
		
		// check convergence
		if (resid < getTolerance()) break;

		// rho_prev = rho
		rho_prev = rho;
	}

	// K.Matsubara 2010.03.31
	P.subst(pX);
	pA->multMPC(&P, pX);

	// interface data exchange
	//pX->updateCommBoundary();
	
   	printf(" --- end of BiCGSTAB solver --- \n");

	return 1;
}

}
