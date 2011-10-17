/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/SolverMG.h
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
