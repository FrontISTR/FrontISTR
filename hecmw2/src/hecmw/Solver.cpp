/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/Solver.cpp
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
void CSolver::setIterMax(iint iter_new)
{
	mIterMax = iter_new;
}
iint CSolver::getIterMax()
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
void CSolver::setMethod(uiint method)
{
	mMethod = method;
}
uiint CSolver::getMethod()
{
	return mMethod;
}
void CSolver::setPrecondition(iint precondition)
{
	mPrecondition = precondition;
}
iint CSolver::getPrecondition()
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
uiint CSolver::solve(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX)
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
