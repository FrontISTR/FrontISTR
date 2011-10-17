/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/SolverGPBiCG.h
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
