/*
 * SolverGMRES.h
 *
 *  Created on: Jul 24, 2009
 *      Author: goto
 */

#ifndef SOLVERGMRES_H_
#define SOLVERGMRES_H_

#include "Solver.h"

namespace pmw
{

class CSolverGMRES: public pmw::CSolver
{
public:
	CSolverGMRES(int iter_max, double tolerance,
			uint method, uint precondition,
			bool flag_iter, bool flag_time);
	virtual ~CSolverGMRES();
private:
	int doSolve(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
			int iter_max, double tolerance,
			bool flag_iter_log, bool flag_time_log);
};

}

#endif /* SOLVERGMRES_H_ */
