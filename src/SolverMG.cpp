/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   SolverMG.cxx
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
CSolverMG::CSolverMG(int iter_max = 100,
		double tolerance = 1.0e-8,
		uint method = 1,
		uint precondition = 1,
		bool flag_iter_log = false,
		bool flag_time_log = false)
	: CSolver(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log)
{
}
CSolverMG::~CSolverMG()
{
}
int CSolverMG::doSolve(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
		int iter_max, double tolerance,
		bool flag_iter_log, bool flag_time_log)
{
	double bnrm2_inv, rnrm2;
	double resid;
	CAssyVector R(pX); 
	bnrm2_inv = 1.0 / pB->norm2();
	for (int iter = 0; iter < getIterMax(); iter++) {
		int mg_iter = 1; 
		int alpha1; 
		int alpha2; 
		pA->MGCycle(pB, pX, mg_iter, alpha1, alpha2);
		pA->residual(pX, pB, &R);  
		rnrm2 = R.norm2();
		resid = sqrt(rnrm2 * bnrm2_inv);
		if (getFlagIterLog()) {
		}
		if (resid < getTolerance()) break;
	}
	pX->updateCommBoundary();
	return 1;
}
}
