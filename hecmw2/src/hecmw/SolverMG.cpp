/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/SolverMG.cpp
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "HEC_MPI.h"
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
}
CSolverMG::~CSolverMG()
{
}
uiint CSolverMG::doSolve(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
		iint iter_max, double tolerance,
		bool flag_iter_log, bool flag_time_log)
{
	double bnrm2_inv, rnrm2;
	double resid;
	iint iter;
	CAssyVector R(pX); 
	bnrm2_inv = 1.0 / pB->norm2();
	for (iter = 0; iter < getIterMax(); iter++) {
		iint mg_iter = 1; 
		iint alpha1; 
		iint alpha2; 
		pA->MGCycle(pB, pX, mg_iter, alpha1, alpha2);
		pA->residual(pX, pB, &R);  
		rnrm2 = R.norm2();
		resid = sqrt(rnrm2 * bnrm2_inv);
		if (getFlagIterLog()) {
		}
		if (resid < getTolerance()) break;
	}
	if( iter == getIterMax() ) return 0;
	pX->updateCommBoundary();
	return 1;
}
}
