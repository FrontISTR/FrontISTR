/*
 * SolverBiCGSTAB.h
 *
 *  Created on: Jul 24, 2009
 *      Author: goto
 */

#ifndef SOLVERBICGSTAB_H_
#define SOLVERBICGSTAB_H_

#include "Solver.h"

namespace pmw
{

class CSolverBiCGSTAB: public pmw::CSolver
{
public:
	CSolverBiCGSTAB(int iter_max, double tolerance,
			uint method, uint precondition,
			bool flag_iter, bool flag_time);
	virtual ~CSolverBiCGSTAB();
private:
	int doSolve(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
			int iter_max, double tolerance,
			bool flag_iter_log, bool flag_time_log);
};

}

#endif /* SOLVERBICGSTAB_H_ */
