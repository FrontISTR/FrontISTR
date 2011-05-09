/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   SolverMG.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef SOLVERMG_H_
#define SOLVERMG_H_
#include "Solver.h"
namespace pmw
{
class CSolverMG: public pmw::CSolver
{
public:
	CSolverMG(int iter_max, double tolerance,
			uint method, uint precondition,
			bool flag_iter, bool flag_time);
	virtual ~CSolverMG();
private:
	int doSolve(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
			int iter_max, double tolerance,
			bool flag_iter_log, bool flag_time_log);
};
}
#endif /* SOLVERMG_H_ */
