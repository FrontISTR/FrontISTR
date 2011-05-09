/*
 * Solver.h
 *
 *  Created on: Jul 23, 2009
 *      Author: goto
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include "TypeDef.h"

namespace pmw
{

class CAssyMatrix;
class CAssyVector;

class CSolver
{
public:
	CSolver(int iter_max, double tolerance,
			uint method, uint precondition,
			bool flag_iter, bool flag_time);
	virtual ~CSolver();

	void setFlagNewMesh(bool flag_new);
	void setFlagNewCoef(bool flag_new);

	void setIterMax(int iter_new);
	int getIterMax();

	void setTolerance(double tol_new);
	double getTolerance();

	void setMethod(uint method);
	uint getMethod();

	void setPrecondition(uint precondition);
	uint getPrecondition();

	void setFlagIterLog(bool flag_new);
	bool getFlagIterLog();

	void setFlagTimeLog(bool flag_new);
	bool getFlagTimeLog();

	int solve(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX);

private:
	virtual int doSolve(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
			int iter_max, double tolerance,
			bool flag_iter_log, bool flag_time_log) = 0;

	int mIterMax;
	double mTolerance;
	uint mMethod;
	uint mPrecondition;
	bool mFlagIterLog;
	bool mFlagTimeLog;
	bool mFlagNewMesh;
	bool mFlagNewCoef;
};

}

#endif /* SOLVER_H_ */
