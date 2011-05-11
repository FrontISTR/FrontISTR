/*
 * SolverMG.cpp
 *
 *  Created on: Dec 29, 2009
 *      Author: goto
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
	// TODO Auto-generated constructor stub

}

CSolverMG::~CSolverMG()
{
	// TODO Auto-generated destructor stub
}

uiint CSolverMG::doSolve(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
		iint iter_max, double tolerance,
		bool flag_iter_log, bool flag_time_log)
{
	double bnrm2_inv, rnrm2;
	double resid;
	iint iter;

	CAssyVector R(pX); // same type as pB and pX

	// calc 1 / |{B}|^2
	bnrm2_inv = 1.0 / pB->norm2();

	for (iter = 0; iter < getIterMax(); iter++) {
		iint mg_iter = 1; // 1: V-Cycle, 2: W-Cycle
		iint alpha1; // number of pre smoothing
		iint alpha2; // number of post smoothing

		pA->MGCycle(pB, pX, mg_iter, alpha1, alpha2);

		// {R} = {B} - [A] {X}
		pA->residual(pX, pB, &R);  // can't skip [A] {X} ???

		// calc |{R}|^2
		rnrm2 = R.norm2();

		// resid = |{R}| / |{B}|
		resid = sqrt(rnrm2 * bnrm2_inv);

		// iteration history
		if (getFlagIterLog()) {
			// TODO: replace with logger
		  //			printf("%5d %16.6e\n", iter + 1, resid);
		}

		// check convergence
		if (resid < getTolerance()) break;
	}
	
	if( iter == getIterMax() ) return 0;

	// interface data exchange
	pX->updateCommBoundary();
	return 1;
}

}
