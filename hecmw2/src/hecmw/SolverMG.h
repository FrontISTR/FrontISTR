/*
 * SolverMG.h
 *
 *  Created on: Dec 29, 2009
 *      Author: goto
 */

#ifndef SOLVERMG_H_
#define SOLVERMG_H_

#include "Solver.h"

namespace pmw
{

class CSolverMG: public pmw::CSolver
{
public:
	CSolverMG(iint iter_max, double tolerance,
			iint method, iint precondition,
			bool flag_iter, bool flag_time);
	virtual ~CSolverMG();
private:
	uiint doSolve(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
			iint iter_max, double tolerance,
			bool flag_iter_log, bool flag_time_log);
};

}

#endif /* SOLVERMG_H_ */
