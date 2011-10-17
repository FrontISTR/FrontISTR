/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/AssyMatrix.h
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
	CAssyMatrix(CAssyModel *pAssyModel, const uiint& nDOF);
	virtual ~CAssyMatrix();
        uiint Matrix_Add_Nodal(const uiint& iMesh, const uiint& iNodeID, const uiint& jNodeID, double* NodalMatrix);
	uiint Matrix_Add_Elem(CAssyModel *pAssyModel, const uiint& iMesh, const uiint& iElem, double *ElemMatrix);
        void Matrix_Clear(const uiint& iMesh);
        uiint& getDOF(){ return mnDOF;}
        void setValue_D(const uiint& imesh, const uiint& inode, const uiint& idof, const double& value);
	void setValue(const uiint& imesh, const uiint& inode, const uiint& idof, const double& dDiag, CAssyVector *pRHS, const double& dRHS);
        void setZero_NonDiag(const uiint& imesh, const uiint& inode, const uiint& idof, CAssyVector *pRHS, const double& dRHS);
        void multVector(const uiint& imesh, CAssyVector *pX, CAssyVector *pB);
	uiint multVector(CAssyVector *pV, CAssyVector *pP, CAssyVector *pW = 0) const;
	uiint multMPC(CAssyVector *pV, CAssyVector *pP) const;
	uiint residual(CAssyVector *pV, const CAssyVector *pF, CAssyVector *pR) const;
	uiint setupSolver(iint type);
	uiint setupPreconditioner(iint type) const;
	uiint setupSmoother(iint type);
	void setCoarseMatrix(CAssyMatrix *coarse) {mpCoarseMatrix = coarse;};
	CAssyMatrix* getCoarseMatrix() { return mpCoarseMatrix;};
	uiint solve(const CAssyVector *pF, CAssyVector *pV, iint iter) const;
	uiint precond(const CAssyVector *pR, CAssyVector *pZ, iint iter) const;
	uiint relax(const CAssyVector *pF, CAssyVector *pV, iint iter) const;
	uiint MGCycle(const CAssyVector *pF, CAssyVector *pV, iint iter, iint alpha1, iint alpha2) const;
	uiint MGInitialGuess(const CAssyVector *pF, CAssyVector *pV) const;
        void dump();
private:
	CAssyModel *mpAssyModel;
        uiint mnDOF;
	std::vector<CMatrixBCRS*> mvMatrix;
	std::vector<CContactMatrix*> mvContactMatrix;
	std::vector<CMPCMatrix*> mvMPCMatrix;
	CSolver *mpSolver;
	CAssyMatrix *mpCoarseMatrix;
	uiint mMGLevel;
};
}
#endif /* ASSYMATRIX_H_ */
