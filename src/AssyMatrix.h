/*
 * AssyMatrix.h
 *
 *  Created on: Oct 16, 2009
 *      Author: goto
 */

#ifndef ASSYMATRIX_H_
#define ASSYMATRIX_H_

#include <vector>

namespace pmw
{

class CMatrixBCRS;
class CContactMatrix;
class CMPCMatrix;
class CAssyVector;
class CSolver;
class CAssyModel;

class CAssyMatrix
{
public:
	CAssyMatrix(/* const */ CAssyModel *pAssyModel);
	virtual ~CAssyMatrix();

	int Matrix_Add_Elem(CAssyModel *pAssyModel, int iMesh, int iElem, double *ElemMatrix);
	void setValue(int imesh, int inode, int idof, double value);
	int multVector(CAssyVector *pV, CAssyVector *pP, CAssyVector *pW = 0) const;
	int multMPC(CAssyVector *pV, CAssyVector *pP) const;
	int residual(CAssyVector *pV, const CAssyVector *pF, CAssyVector *pR) const;

	int setupSolver(int type);
	int setupPreconditioner(int type) const; 
	int setupSmoother(int type);
	
	void setCoarseMatrix(CAssyMatrix *coarse) {mpCoarseMatrix = coarse;};
	CAssyMatrix* getCoarseMatrix() { return mpCoarseMatrix;};

	int solve(const CAssyVector *pF, CAssyVector *pV, int iter) const;
	int precond(const CAssyVector *pR, CAssyVector *pZ, int iter) const;
	int relax(const CAssyVector *pF, CAssyVector *pV, int iter) const;

	int MGCycle(const CAssyVector *pF, CAssyVector *pV, int iter, int alpha1, int alpha2) const;
	int MGInitialGuess(const CAssyVector *pF, CAssyVector *pV) const;

private:
	CAssyModel *mpAssyModel;

	std::vector<CMatrixBCRS*> mvMatrix;
	std::vector<CContactMatrix*> mvContactMatrix;
	std::vector<CMPCMatrix*> mvMPCMatrix;

	CSolver *mpSolver;

	CAssyMatrix *mpCoarseMatrix;
	unsigned int mMGLevel;
};

}

#endif /* ASSYMATRIX_H_ */
