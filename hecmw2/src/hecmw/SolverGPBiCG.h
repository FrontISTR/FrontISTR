/*
 * SolverGPBiCG.h
 *
 *  Created on: Jul 24, 2009
 *      Author: goto
 */

#ifndef SOLVERGPBICG_H_
#define SOLVERGPBICG_H_

#include "Solver.h"

namespace pmw
{

class CSolverGPBiCG: public pmw::CSolver
{
public:
	CSolverGPBiCG(iint iter_max, double tolerance,
			iint method, iint precondition,
			bool flag_iter, bool flag_time);
	virtual ~CSolverGPBiCG();
private:
	uiint doSolve(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
			iint iter_max, double tolerance,
			bool flag_iter_log, bool flag_time_log);
};

}

#endif /* SOLVERGPBICG_H_ */
