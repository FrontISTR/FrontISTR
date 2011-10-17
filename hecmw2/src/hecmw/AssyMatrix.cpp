/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/AssyMatrix.cpp
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
#include "HEC_MPI.h"
#include "AssyMatrix.h"
#include "MatrixBCRS.h"
#include "ContactMatrix.h"
#include "MPCMatrix.h"
#include "AssyVector.h"
#include "Solver.h"
#include "AssyModel.h"
#include "Equation.h"
namespace pmw
{
typedef std::vector<CMatrixBCRS*> CVMat;
typedef CVMat::iterator CVMatIter;
typedef CVMat::const_iterator CVMatConstIter;
CAssyMatrix::CAssyMatrix(CAssyModel *pAssyModel, const uiint& nDOF)
{
	mpAssyModel = pAssyModel;
        mnDOF = nDOF;
	for (int i = 0; i < pAssyModel->getNumOfMesh(); i++) {
		CMesh *pMesh = pAssyModel->getMesh(i);
		CMatrixBCRS *pMatrix = new CMatrixBCRS(pMesh, mnDOF);
		mvMatrix.push_back( pMatrix );
	}
	mMGLevel = pAssyModel->getMGLevel();
	uiint numOfContact= pAssyModel->getNumOfContactMesh();
        uiint idof;
        CEquation **vEquation;
        vEquation = new CEquation*[mnDOF];
	for(uiint icont = 0; icont < numOfContact; icont++){
		CContactMesh* pConMesh= pAssyModel->getContactMesh(icont);
		uiint numOfSPoint = pConMesh->getNumOfSlavePoint();
		CMPCMatrix* mpc = new CMPCMatrix();
		for(uiint islave = 0; islave< numOfSPoint; islave++){
			CContactNode* pSlaveNode = pConMesh->getSlaveConNode(islave);
			uiint mgLevel = mMGLevel;
			if(pSlaveNode->have_MasterFaceID(mgLevel)){
                                for(idof=0; idof < mnDOF; idof++) vEquation[idof]= new CEquation();
				uiint masterFaceID = pSlaveNode->getMasterFaceID(mgLevel);
				uiint smesh = pSlaveNode->getMeshID(); 
                                uiint snode = pSlaveNode->getNodeID();
                                CMesh* pSMesh= pAssyModel->getMesh_ID(smesh);
                                CIndexBucket *pSBucket= pSMesh->getBucket();
                                uiint snode_ix = pSBucket->getIndexNode(snode);
                                for(idof=0; idof < mnDOF; idof++) vEquation[idof]->setSlave(smesh, snode_ix, idof);
				CSkinFace* pMasterFace = pConMesh->getMasterFace_ID(masterFaceID);
				uiint numOfVert = pMasterFace->getNumOfNode();
				for(uiint ivert=0; ivert< numOfVert; ivert++){
					uiint mmesh = pMasterFace->getMeshID(); 
					uiint node = ( pMasterFace->getNode(ivert) )->getNodeID();
                                        CMesh* pMMesh= pAssyModel->getMesh_ID(mmesh);
                                        CIndexBucket *pMBucket= pMMesh->getBucket();
                                        uiint node_ix = pMBucket->getIndexNode(node);
                                        snode = pSlaveNode->getID();
                                        snode_ix = pSBucket->getIndexNode(snode);
					double coef = pMasterFace->getCoef(snode_ix, ivert);
                                        for(idof=0; idof < mnDOF; idof++) vEquation[idof]->addMaster(mmesh, node_ix, idof, coef);
				};
                                for(idof=0; idof < mnDOF; idof++) mpc->addEquation(vEquation[idof]);
			};
		};
		mvMPCMatrix.push_back(mpc);
	};
}
CAssyMatrix::~CAssyMatrix()
{
	for (CVMatIter im = mvMatrix.begin(); im != mvMatrix.end(); im++) {
		delete *im;
	};
}
uiint CAssyMatrix::Matrix_Add_Nodal(const uiint& iMesh, const uiint& inode, const uiint& jnode, double* NodalMatrix)
{
    mvMatrix[iMesh]->Matrix_Add_Nodal(inode, jnode, NodalMatrix);
    return 1;
}
uiint CAssyMatrix::Matrix_Add_Elem(CAssyModel *pAssyModel, const uiint& iMesh, const uiint& iElem, double *ElemMatrix)
{
#ifdef ADVANCESOFT_DEBUG
   	printf(" enter CAssyMatrix::Matrix_Add_Elemm %d %d %e \n", iMesh, iElem, ElemMatrix[0]);
#endif
	CMesh *pMesh = pAssyModel->getMesh(iMesh);
	mvMatrix[iMesh]->Matrix_Add_Elem(pMesh, iElem, ElemMatrix);
#ifdef ADVANCESOFT_DEBUG
   	printf(" exit CAssyMatrix::Matrix_Add_Elemm \n");
#endif
    return 1;
}
void CAssyMatrix::Matrix_Clear(const uiint& iMesh)
{
    mvMatrix[iMesh]->Matrix_Clear();
}
void CAssyMatrix::setValue_D(const uiint& imesh, const uiint& inode, const uiint& idof, const double& value)
{
    mvMatrix[imesh]->setValue_D(inode, idof, value);
}
void CAssyMatrix::setValue(const uiint& imesh, const uiint& inode, const uiint& idof, const double& dDiag, CAssyVector *pRHS, const double& dRHS)
{
#ifdef ADVANCESOFT_DEBUG
    printf(" enter CAssyMatrix::setValue %d %d %d %e \n", imesh, inode, idof, dDiag);
#endif
    mvMatrix[imesh]->setValue(inode, idof, dDiag, pRHS->getVector(imesh), dRHS);
#ifdef ADVANCESOFT_DEBUG
    printf(" exit CAssyMatrix::setValue %d %d %d %e \n", imesh, inode, idof, dDiag);
#endif
}
void CAssyMatrix::setZero_NonDiag(const uiint& imesh, const uiint& inode, const uiint& idof, CAssyVector *pRHS, const double& dRHS)
{
    mvMatrix[imesh]->setZero_NonDiag(inode, idof, pRHS->getVector(imesh), dRHS);
}
void CAssyMatrix::multVector(const uiint& imesh, CAssyVector* pX, CAssyVector* pB)
{
    mvMatrix[imesh]->multVector(pX->getVector(imesh), pB->getVector(imesh));
}
uiint CAssyMatrix::multVector(CAssyVector *pV, CAssyVector *pP, CAssyVector *pW) const
{
	pV->updateCommBoundary(); 
	if (mvMPCMatrix.size() == 0) {
		for (uiint i = 0; i < mvMatrix.size(); i++) {
			mvMatrix[i]->multVector(pV->getVector(i), pP->getVector(i));
		}
		if (mvContactMatrix.size() > 0) {
		}
	} else {
		if (pW == 0) {
			pW = new CAssyVector(pP);
		}
		for (uiint i = 0; i < mvMPCMatrix.size(); i++) {
			mvMPCMatrix[i]->multVector(pV, pP);
		}
		for (uiint i = 0; i < mvMatrix.size(); i++) {
			mvMatrix[i]->multVector(pP->getVector(i), pW->getVector(i));
		}
		if (mvContactMatrix.size() > 0) {
			for (uiint i = 0; i < mvContactMatrix.size(); i++) {
			}
		}
		for (uiint i = 0; i < mvMPCMatrix.size(); i++) {
			mvMPCMatrix[i]->transMultVector(pW, pP);
		}
	}
	pP->sumupCommBoundary(); 
	return 1;
}
uiint CAssyMatrix::multMPC(CAssyVector *pV, CAssyVector *pP) const
{
	if (mvMPCMatrix.size() == 0) {
	} else {
		for (uiint i = 0; i < mvMPCMatrix.size(); i++) {
			mvMPCMatrix[i]->multVector(pV, pP);
		}
	}
	return 1;
}
uiint CAssyMatrix::residual(CAssyVector *pV, const CAssyVector *pF, CAssyVector *pR) const
{
    multVector(pV, pR);
    for (uiint i = 0; i < pR->size(); i++) {
            (*pR)[i] = (*pF)[i] - (*pR)[i];
    }
    return 1;
}
uiint CAssyMatrix::setupSolver(iint type)
{
    return 1;
}
uiint CAssyMatrix::setupPreconditioner(iint type) const
{
#ifdef ADVANCESOFT_DEBUG
	printf(" enter CAssyMatrix::setupPreconditioner %d \n", type);
#endif
	for (CVMatConstIter im = mvMatrix.begin(); im != mvMatrix.end(); im++) {
		(*im)->setupPreconditioner(type);
	}
	if( mpCoarseMatrix ) mpCoarseMatrix->setupPreconditioner(type);
#ifdef ADVANCESOFT_DEBUG
	printf(" exit CAssyMatrix::setupPreconditioner \n");
#endif
	return 1;
}
uiint CAssyMatrix::setupSmoother(iint type)
{
	for (CVMatIter im = mvMatrix.begin(); im != mvMatrix.end(); im++) {
		(*im)->setupSmoother(type);
	}
	return 1;
}
uiint CAssyMatrix::solve(const CAssyVector *pF, CAssyVector *pV, iint iter = 0) const
{
	int iter_save = mpSolver->getIterMax();
	if (iter > 0) mpSolver->setIterMax(iter);
	mpSolver->solve(this, pF, pV);
	if (iter > 0) mpSolver->setIterMax(iter_save);
	return 1;
}
uiint CAssyMatrix::precond(const CAssyVector *pR, CAssyVector *pZ, iint iter) const
{
	CAssyVector z_resid(pR);
	pZ->setZero();
	for (iint i = 0; i < iter; i++) {
		uiint index = 0;
		for (CVMatConstIter im = mvMatrix.begin(); im != mvMatrix.end(); im++) {
			(*im)->precond(pR->getVector(index), z_resid.getVector(index));
			index++;
		}
		pZ->add(&z_resid);
		if (i == iter - 1) break;
		residual(pZ, pR, &z_resid);
	}
	return 1;
}
uiint CAssyMatrix::relax(const CAssyVector *pF, CAssyVector *pV, iint iter) const
{
	CAssyVector v_resid(pF);
	pV->setZero();
	for (iint i = 0; i < iter; i++) {
		uiint index = 0;
		for (CVMatConstIter im = mvMatrix.begin(); im != mvMatrix.end(); im++) {
			(*im)->relax(pF->getVector(index), v_resid.getVector(index));
			index++;
		}
		pV->add(&v_resid);
		if (i == iter - 1) break;
		residual(pV, pF, &v_resid);
	}
	return 1;
}
uiint CAssyMatrix::MGCycle(const CAssyVector *pF, CAssyVector *pV, iint iter, iint alpha1, iint alpha2) const
{
	relax(pF, pV, alpha1);
	if (mpCoarseMatrix != 0) {
		CAssyVector w(mpAssyModel, mnDOF);
		CAssyVector fc(mpCoarseMatrix->mpAssyModel, mnDOF);
		CAssyVector vc(mpCoarseMatrix->mpAssyModel, mnDOF);
		residual(pV, pF, &w);
		w.restrictTo(&fc);
		vc.setZero();
		for (iint i = 0; i < iter; i++) {
                    mpCoarseMatrix->MGCycle(&vc, &fc, iter, alpha1, alpha2);
		}
		w.prolongateFrom(&vc);
		pV->add(&w);
	} else {
	}
	relax(pF, pV, alpha2);
	return 1;
}
uiint CAssyMatrix::MGInitialGuess(const CAssyVector *pF, CAssyVector *pV) const
{
	printf(" enter CAssyMatrix::MGInitialGuess %ld \n", mMGLevel);
	int mAlpha; 
        if (mpCoarseMatrix != NULL) {
		CAssyVector w(mpAssyModel, mnDOF);
		CAssyVector fc(mpCoarseMatrix->mpAssyModel, mnDOF);
		CAssyVector vc(mpCoarseMatrix->mpAssyModel, mnDOF);
		residual(pV, pF, &w);
		w.restrictTo(&fc);
		vc.setZero();
		mpCoarseMatrix->MGInitialGuess(&fc, &vc);
		w.prolongateFrom(&vc);
		pV->add(&w);
	} else {
		precond(pF, pV, 1);
	}
	precond(pF, pV, 1);
	printf(" exit CAssyMatrix::MGInitialGuess %ld \n", mMGLevel);
	return 1;
}
void CAssyMatrix::dump()
{
    uiint i;
    for(i=0; i < mvMatrix.size(); i++){
        mvMatrix[i]->dump();
    };
    for(i=0; i < mvMPCMatrix.size(); i++){
        mvMPCMatrix[i]->dump();
    };
}
}
