/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/AssyMatrix.h
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "CommonStd.h"
#include "TypeDef.h"

#include "SolverCG.h"
#include "SolverBiCGSTAB.h"
#include "SolverGPBiCG.h"
#include "SolverGMRES.h"

#include "HEC_MPI.h"
namespace pmw
{
#ifndef ASSYMATRIX_H_
#define ASSYMATRIX_H_
class CMatrixBCRS;
class CContactMatrix;
class CMPCMatrix;
class CAssyVector;
class CAssyModel;

class CAssyMatrix
{
public:
    CAssyMatrix(CAssyModel *pAssyModel, vuint& vDOF, const uiint& ieq);
    virtual ~CAssyMatrix();
    uiint Matrix_Add_Nodal(const uiint& iMesh, const uiint& iNodeID, const uiint& jNodeID, double* NodalMatrix);
    uiint Matrix_Add_Elem(CAssyModel *pAssyModel, const uiint& iMesh, const uiint& iElem, double *ElemMatrix);
    void Matrix_Clear(const uiint& iMesh);
    void Matrix_Clear();

    void setValue_D(const uiint& imesh, const uiint& inode, const uiint& idof, const double& value);
    void setValue(const uiint& imesh, const uiint& inode, const uiint& idof, const double& dDiag, CAssyVector *pRHS, const double& dRHS);
    void setZero_NonDiag(const uiint& imesh, const uiint& inode, const uiint& idof, CAssyVector *pRHS, const double& dRHS);

    uiint multVector(CAssyVector *pV, CAssyVector *pP, CAssyVector *pW = 0) const;

    uiint multMPC(CAssyVector *pV, CAssyVector *pP) const;
    uiint residual(CAssyVector *pV, const CAssyVector *pF, CAssyVector *pR) const;

    uiint setupSolver(iint iter_max, double tolerance, iint method, iint precondition, bool flag_iter_log,bool flag_time_log);
    void deleteSolver();

    uiint setupPreconditioner(iint type) const;
    uiint setupSmoother(iint type);

    //--
    // MG
    //--
    void setCoarseMatrix(CAssyMatrix *coarse) {
        mpCoarseMatrix = coarse;
    };
    CAssyMatrix* getCoarseMatrix() {
        return mpCoarseMatrix;
    };
    CAssyMatrix* getCoarseMatrix() const {
        return mpCoarseMatrix;
    };

    void restrict();//---- coarse matrix 生成 => Restriction 行列

    uiint precond(const CAssyVector *pR, CAssyVector *pZ, iint iter) const;
    uiint relax(const CAssyVector *pF, CAssyVector *pV, iint iter) const;
    uiint MGCycle(const CAssyVector *pF, CAssyVector *pV, iint mg_cycle, iint alpha1, iint alpha2);// const;

    uiint MGInitialGuess(const CAssyVector *pF, CAssyVector *pV);// const;

    void dump();
    void dump() const;

    CMatrixBCRS* getMatrix(const uiint& index) {
        return mvMatrix[index];
    }

private:
    CAssyModel *mpAssyModel;

    std::vector<CMatrixBCRS*> mvMatrix;
    std::vector<CContactMatrix*> mvContactMatrix;
    std::vector<CMPCMatrix*> mvMPCMatrix;

    CSolver *mpSolver;
    CAssyMatrix *mpCoarseMatrix;//---- コースグリッドAssyMatrix

    uiint mMGLevel;
};
#endif /* ASSYMATRIX_H_ */
}

