/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Solver.cxx
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
CSolver::CSolver(int iter_max = 100,
		double tolerance = 1.0e-8,
		uint method = 1,
		uint precondition = 1,
		bool flag_iter_log = false,
		bool flag_time_log = false)
	: mIterMax(iter_max),
	  mTolerance(tolerance),
	  mMethod(method),
	  mPrecondition(precondition),
	  mFlagIterLog(flag_iter_log),
	  mFlagTimeLog(flag_time_log)
{
	mFlagNewMesh = true;
	mFlagNewCoef = true;
}
CSolver::~CSolver()
{
}
void CSolver::setFlagNewMesh(bool flag_new)
{
	mFlagNewMesh = flag_new;
}
void CSolver::setFlagNewCoef(bool flag_new)
{
	mFlagNewCoef = flag_new;
}
void CSolver::setIterMax(int iter_new)
{
	mIterMax = iter_new;
}
int CSolver::getIterMax()
{
	return mIterMax;
}
void CSolver::setTolerance(double tol_new)
{
	mTolerance = tol_new;
}
double CSolver::getTolerance()
{
	return mTolerance;
}
void CSolver::setMethod(uint method)
{
	mMethod = method;
}
uint CSolver::getMethod()
{
	return mMethod;
}
void CSolver::setPrecondition(uint precondition)
{
	mPrecondition = precondition;
}
uint CSolver::getPrecondition()
{
	return mPrecondition;
}
bool CSolver::getFlagIterLog()
{
	return mFlagIterLog;
}
void CSolver::setFlagIterLog(bool flag_new)
{
	mFlagIterLog = true;
}
void CSolver::setFlagTimeLog(bool flag_new)
{
	mFlagTimeLog = true;
}
bool CSolver::getFlagTimeLog()
{
	return mFlagTimeLog;
}
int CSolver::solve(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX)
{
	if (mFlagNewMesh) {
		mFlagNewMesh = false;
		mFlagNewCoef = false;
	} else if (mFlagNewCoef) {
		mFlagNewCoef = false;
	}
	doSolve(pA, pB, pX,
		mIterMax, mTolerance,
		mFlagIterLog, mFlagTimeLog);
	return 1;
}
}
