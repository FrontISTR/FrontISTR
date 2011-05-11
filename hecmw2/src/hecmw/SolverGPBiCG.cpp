/*
 * SolverGPBiCG.cpp
 *
 *  Created on: Jul 24, 2009
 *      Author: goto
 */

#include "SolverGPBiCG.h"
#include "AssyMatrix.h"
#include "AssyVector.h"
#include <cmath>
#include <cstdio>

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
	// TODO Auto-generated constructor stub

}

CSolverGPBiCG::~CSolverGPBiCG()
{
	// TODO Auto-generated destructor stub
}

uiint CSolverGPBiCG::doSolve(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
		iint iter_max, double tolerance,
		bool flag_iter_log, bool flag_time_log)
{
	double bnrm2_inv;
	double alpha, beta, rho, rho_prev, resid;
	iint iter, iter_precond; // TODO: class variable???

   	printf(" --- start of GPBiCG solver --- \n");

	CAssyVector R (pX); // same type as pB and pX
	CAssyVector Rt(pX); //
	CAssyVector Rp(pX); //
	CAssyVector T (pX); //
	CAssyVector Tt(pX); //
	CAssyVector T0(pX); //
	CAssyVector P (pX); //
	CAssyVector Pt(pX); //
	CAssyVector U (pX); //
	CAssyVector W1(pX); //
	CAssyVector Y (pX); //
	CAssyVector Z (pX); //
	CAssyVector WK(pX); //
	CAssyVector W2(pX); //

	uiint len = pB->size();

	// {R} = {B} - [A] {X}
	// {Rt} = {R}
	pA->multVector(pX, &P);
	for (int i = 0; i < len; i++) {
		R[i] = (*pB)[i] - P[i];
		Rt[i] = R[i];
	}

	// calc 1 / |{B}|^2
	bnrm2_inv = 1.0 / pB->norm2();

	// rho = {Rt} {R}
	rho = Rt.innerProd(&R);

	iint itype = getPrecondition();
	pA->setupPreconditioner( itype ); // TODO: rewrite to pA->precond(pR, pZ);

	for (iter = 0; iter < getIterMax(); iter++) {

		// {Rp} = [Minv] {R}
		iter_precond = 1;

		if( itype == 2 ) {
		  Rp.subst(pX);
		  pA->MGInitialGuess(&R, &Rp);
		} else {
		  pA->precond(&R, &Rp, iter_precond);
		}
		// iter > 0:
		//   {P} = {Rp} + beta ( {P} - {U} )
		// iter == 1:
		//   {P} = {Rp}
		if (iter > 0) {
			for (uiint i = 0; i < len; i++) {
				P[i] = Rp[i] + beta * (P[i] - U[i]);
			}
		} else {
			for (uiint i = 0; i < len; i++) {
				P[i] = Rp[i];
			}
		}

		// {Pt} = [A] {P}
		pA->multVector(&P, &Pt);

		// alpha = rho / ({Rt} {Pt})
		alpha = rho / Rt.innerProd(&Pt);

		// {Y} = {T} - {R} + alpha (- {W1} + {Pt})
		// {T} = {R} - alpha {Pt}
		for (uiint i = 0; i < len; i++) {
			Y[i] = T[i] - R[i] + alpha * (- W1[i] + Pt[i]);
			T[i] = R[i] - alpha * Pt[i];
		}

		// {Tt} = [Minv] {T}
		iter_precond = 1;

		if( itype == 2 ) {
		  Tt.subst(pX);
		  pA->MGInitialGuess(&T, &Tt);
		} else {
		  pA->precond(&T, &Tt, iter_precond);
		}
		// {T} = [Minv] {T}
		iter_precond = 1;

		if( itype == 2 ) {
		  W2.subst(pX);
		  pA->MGInitialGuess(&T0, &W2);
		} else {
		  pA->precond(&T0, &W2, iter_precond);
		}
		for (int i = 0; i < len; i++) {
			T0[i] = W2[i];
		}

		// {W2} = [Minv] {Pt}
		iter_precond = 1;

		if( itype == 2 ) {
		  W2.subst(pX);
		  pA->MGInitialGuess(&Pt, &W2);
		} else {
		  pA->precond(&Pt, &W2, iter_precond);
		}
		// {Tt} = [A] {Tt}
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
		// iter > 0:
		//   qsi =
		//   eta =
		// iter == 1:
		//   qsi =
		//   eta =
		if (iter > 0) {
			double denom = c[4] * c[0] - c[3] * c[3];
			qsi = (c[0] * c[1] - c[2] * c[3]) / denom;
			eta = (c[4] * c[2] - c[3] * c[1]) / denom;
		} else {
			qsi = c[1] / c[4];
			eta = 0.0;
		}
		// {U} = qsi {W2} + eta ({T0} + {Rp} + beta {U})
		// {Z} = qsi {Rp} - eta {Z} - alpha {U}
		for (uiint i = 0; i < len; i++) {
			U[i] = qsi * W2[i] + eta * (T0[i] - Rp[i] + beta * U[i]);
			Z[i] = qsi * Rp[i] + eta * Z[i] - alpha * U[i];
		}

		// update {X}, {R}, {W}
		for (uiint i = 0; i < len; i++) {
			(*pX)[i] += alpha * P[i] + Z[i];
			R[i] = T[i] - eta * Y[i] - qsi * Tt[i];
			T0[i] = T[i];
		}

		// calc |{R}|^2
		double rnrm2 = R.norm2();

		// rho_prev = rho
		rho_prev = rho;

		// rho = {Rt} {R}
		rho = Rt.innerProd(&R);

		// beta = (alpha * rho) / (qsi * rho_prev)
		beta = alpha * rho / (qsi * rho_prev);

		// {W1} = {Tt} + beta {Pt}
		for (uiint i = 0; i < len; i++) {
			W1[i] = Tt[i] + beta * Pt[i];
		}

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
	}
	
	if( iter == getIterMax() ) return 0;

	// K.Matsubara 2010.03.31
	P.subst(pX);
	pA->multMPC(&P, pX);

	// interface data exchange
	//pX->updateCommBoundary();

   	printf(" --- end of GPBiCG solver --- \n");

	return 1;
}

}
