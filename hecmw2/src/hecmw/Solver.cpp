/*
 * Solver.cpp
 *
 *  Created on: Jul 23, 2009
 *      Author: goto
 */

#include "Solver.h"
#include "AssyMatrix.h"
#include "AssyVector.h"

namespace pmw
{

CSolver::CSolver(iint iter_max = 100,
		double tolerance = 1.0e-8,
		iint method = 1,
		iint precondition = 1,
		bool flag_iter_log = false,
		bool flag_time_log = false)
	: mIterMax(iter_max),
	  mTolerance(tolerance),
	  mMethod(method),
	  mPrecondition(precondition),
	  mFlagIterLog(flag_iter_log),
	  mFlagTimeLog(flag_time_log)
{
	// TODO Auto-generated constructor stub
	mFlagNewMesh = true;
	mFlagNewCoef = true;
}

CSolver::~CSolver()
{
	// TODO Auto-generated destructor stub
}

void CSolver::setFlagNewMesh(bool flag_new)
{
	mFlagNewMesh = flag_new;
}

void CSolver::setFlagNewCoef(bool flag_new)
{
	mFlagNewCoef = flag_new;
}

// IterMax

void CSolver::setIterMax(iint iter_new)
{
	mIterMax = iter_new;
}

iint CSolver::getIterMax()
{
	return mIterMax;
}

// Tolerance

void CSolver::setTolerance(double tol_new)
{
	mTolerance = tol_new;
}

double CSolver::getTolerance()
{
	return mTolerance;
}

// Method

void CSolver::setMethod(uiint method)
{
	mMethod = method;
}

uiint CSolver::getMethod()
{
	return mMethod;
}

// Precondition

void CSolver::setPrecondition(iint precondition)
{
	mPrecondition = precondition;
}

iint CSolver::getPrecondition()
{
	return mPrecondition;
}

// FlagIterLog

bool CSolver::getFlagIterLog()
{
	return mFlagIterLog;
}

void CSolver::setFlagIterLog(bool flag_new)
{
	mFlagIterLog = true;
}

// FlagTimeLog

void CSolver::setFlagTimeLog(bool flag_new)
{
	mFlagTimeLog = true;
}

bool CSolver::getFlagTimeLog()
{
	return mFlagTimeLog;
}

// solve

uiint CSolver::solve(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX)
{
	// double time_start;
	// double time_comm_setup;
	// double time_comm_solve;
	// double time_setup;
	// double time_solve;

//	if (mFlagTimeLog) {
//		time_start = pmw::Wtime();
//	}

	//TODO: error check
	// ZERO HRS norm
	// ZERO DIAGONAL component

	if (mFlagNewMesh) {
		// pA->setupPrecond(); <-- arguments???
		mFlagNewMesh = false;
		mFlagNewCoef = false;
	} else if (mFlagNewCoef) {
		// pA->setupPrecondCoef(); <-- arguments???
		mFlagNewCoef = false;
	}

//	if (mFlagTimeLog) {
//		time_setup = pmw::Wtime() - time_start;
//	}

	doSolve(pA, pB, pX,
		mIterMax, mTolerance,
		mFlagIterLog, mFlagTimeLog);

//	if (mFlagTimeLog) {
//		time_solve = pmw::Wtime() - time_start;
//		double work_ratio = (time_solve - time_comm_solve) / (time_solve + 1.0e-24) * 100.0;
//		// TODO: replace with logger
//		printf("### summary of linear solver\n");
//		printf("%10d iterations %16.6e\n", iter, resid);
//		printf("set-up time     : %16.6e\n", time_setup);
//		printf("solver time     : %16.6e\n", time_solve);
//		printf("solver/comm time: %16.6e\n", time_comm_solve);
//		printf("work ratio(%)   : %16.6g", work_ratio);
//	}
	return 1;
}

}
