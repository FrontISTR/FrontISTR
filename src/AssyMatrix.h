/*
 * AssyMatrix.h
 *
 *  Created on: Oct 16, 2009
 *      Author: goto
 */

#ifndef ASSYMATRIX_H_
#define ASSYMATRIX_H_

#include "CommonStd.h"
#include "TypeDef.h"

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
	CAssyMatrix(CAssyModel *pAssyModel, const uint& nDOF);//, const vuint& vPartsID);//第二引数：アセンブルするメッシュパーツ指定
	virtual ~CAssyMatrix();

        int Matrix_Add_Nodal(const uint& iMesh, const uint& iNodeID, const uint& jNodeID, double* NodalMatrix);
	int Matrix_Add_Elem(CAssyModel *pAssyModel, const uint& iMesh, const uint& iElem, double *ElemMatrix);
        void Matrix_Clear(const uint& iMesh);// Matrix 0 clear

        uint& getDOF(){ return mnDOF;}

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
        uint mnDOF;//アセンブル方程式に対して、一つのDOF

	std::vector<CMatrixBCRS*> mvMatrix;
	std::vector<CContactMatrix*> mvContactMatrix;
	std::vector<CMPCMatrix*> mvMPCMatrix;

	CSolver *mpSolver;

	CAssyMatrix *mpCoarseMatrix;
	uint mMGLevel;
};

}

#endif /* ASSYMATRIX_H_ */
