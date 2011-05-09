/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   AssyMatrix.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
/*
 * AssyMatrix.cpp
 *
 *  Created on: Oct 16, 2009
 *      Author: goto
 */
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
CAssyMatrix::CAssyMatrix(/* const */ CAssyModel *pAssyModel)
{
	mpAssyModel = pAssyModel;
	for (int i = 0; i < pAssyModel->getNumOfMesh(); i++) {
		CMesh *pMesh = pAssyModel->getMesh(i);
		CMatrixBCRS *pMatrix = new CMatrixBCRS(pMesh);
		mvMatrix.push_back( pMatrix );
	}
	mMGLevel = pAssyModel->getMGLevel();
	printf("mMGLevel %d \n",mMGLevel);
	uint numOfContact= pAssyModel->getNumOfContactMesh();
	printf("numOfContact %d \n",numOfContact);
	for(uint icont = 0; icont < numOfContact; icont++){
		CContactMesh* pConMesh= pAssyModel->getContactMesh(icont);
		uint numOfSPoint = pConMesh->getNumOfSlavePoint();
		printf("numOfSPoint, icont %d %d \n", numOfSPoint, icont);
		CMPCMatrix* mpc = new CMPCMatrix();
		for(uint islave = 0; islave< numOfSPoint; islave++){
			CContactNode* pSlaveNode = pConMesh->getSlaveConNode(islave);
			uint mgLevel = mMGLevel;
			if(pSlaveNode->have_MasterFaceID(mgLevel)){
				CEquation *equation0 = new CEquation();
				CEquation *equation1 = new CEquation();
				CEquation *equation2 = new CEquation();
				int masterFaceID = pSlaveNode->getMasterFaceID(mgLevel);
				int smesh = pSlaveNode->getMeshID(); 
				int snode = pSlaveNode->getNodeID();
				printf("masterFaceID, islave, smesh, node %d %d %d %d\n"
					,masterFaceID, islave, smesh, snode);
				equation0->setSlave(smesh, snode, 0);
				equation1->setSlave(smesh, snode, 1);
				equation2->setSlave(smesh, snode, 2);
				CSkinFace* pMasterFace = pConMesh->getMasterFace_ID(masterFaceID);
				int numOfVert= (int) pMasterFace->getNumOfNode();
				printf("numOfVert %d \n",numOfVert);
				for(uint ivert=0; ivert< numOfVert; ivert++){
					int mmesh = pMasterFace->getMeshID(); 
					int node = ( pMasterFace->getNode(ivert) )->getNodeID();
					double coef = pMasterFace->getCoef(pSlaveNode->getID(),ivert);
					equation0->addMaster(mmesh, node, 0, coef);
					equation1->addMaster(mmesh, node, 1, coef);
					equation2->addMaster(mmesh, node, 2, coef);
					printf("ivert, mmesh, node, coef %d %d %d %e \n", ivert, mmesh, node, coef);
				};
				mpc->addEquation(equation0);
				mpc->addEquation(equation1);
				mpc->addEquation(equation2);
			};
		};
		mvMPCMatrix.push_back(mpc);
		printf("mvMPCMatrix.size() %d \n", mvMPCMatrix.size());
	};
}
CAssyMatrix::~CAssyMatrix()
{
	for (CVMatIter im = mvMatrix.begin(); im != mvMatrix.end(); im++) {
		delete *im;
	}
}
int CAssyMatrix::Matrix_Add_Elem(CAssyModel *pAssyModel, int iMesh, int iElem, double *ElemMatrix)
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
void CAssyMatrix::setValue(int imesh, int inode, int idof, double value)
{
#ifdef ADVANCESOFT_DEBUG
	printf(" enter CAssyMatrix::setValue %d %d %d %e \n", imesh, inode, idof, value);
#endif
	mvMatrix[imesh]->setValue(inode, idof, value);
#ifdef ADVANCESOFT_DEBUG
	printf(" exit CAssyMatrix::setValue %d %d %d %e \n", imesh, inode, idof, value);
#endif
}
int CAssyMatrix::multVector(CAssyVector *pV, CAssyVector *pP, CAssyVector *pW) const
{
	pV->updateCommBoundary(); 
	if (mvMPCMatrix.size() == 0) {
		for (int i = 0; i < mvMatrix.size(); i++) {
			mvMatrix[i]->multVector(pV->getVector(i), pP->getVector(i));
		}
		if (mvContactMatrix.size() > 0) {
		}
	} else {
		if (pW == 0) {
			pW = new CAssyVector(pP);
		}
		for (int i = 0; i < mvMPCMatrix.size(); i++) {
			mvMPCMatrix[i]->multVector(pV, pP);
		}
		for (int i = 0; i < mvMatrix.size(); i++) {
			mvMatrix[i]->multVector(pP->getVector(i), pW->getVector(i));
		}
		if (mvContactMatrix.size() > 0) {
			for (int i = 0; i < mvContactMatrix.size(); i++) {
			}
		}
		for (int i = 0; i < mvMPCMatrix.size(); i++) {
			mvMPCMatrix[i]->transMultVector(pW, pP);
		}
	}
	pP->sumupCommBoundary(); 
	return 1;
}
int CAssyMatrix::multMPC(CAssyVector *pV, CAssyVector *pP) const
{
	if (mvMPCMatrix.size() == 0) {
	} else {
		for (int i = 0; i < mvMPCMatrix.size(); i++) {
			mvMPCMatrix[i]->multVector(pV, pP);
		}
	}
	return 1;
}
int CAssyMatrix::residual(CAssyVector *pV, const CAssyVector *pF, CAssyVector *pR) const
{
	multVector(pV, pR);
	for (int i = 0; i < pR->size(); i++) {
		(*pR)[i] = (*pF)[i] - (*pR)[i];
	}
	return 1;
}
int CAssyMatrix::setupSolver(int type)
{
	return 1;
}
int CAssyMatrix::setupPreconditioner(int type) const
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
int CAssyMatrix::setupSmoother(int type)
{
	for (CVMatIter im = mvMatrix.begin(); im != mvMatrix.end(); im++) {
		(*im)->setupSmoother(type);
	}
	return 1;
}
int CAssyMatrix::solve(const CAssyVector *pF, CAssyVector *pV, int iter = 0) const
{
	int iter_save = mpSolver->getIterMax();
	if (iter > 0) mpSolver->setIterMax(iter);
	mpSolver->solve(this, pF, pV);
	if (iter > 0) mpSolver->setIterMax(iter_save);
	return 1;
}
int CAssyMatrix::precond(const CAssyVector *pR, CAssyVector *pZ, int iter) const
{
	CAssyVector z_resid(pR);
	pZ->setZero();
	for (int i = 0; i < iter; i++) {
		int index = 0;
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
int CAssyMatrix::relax(const CAssyVector *pF, CAssyVector *pV, int iter) const
{
	CAssyVector v_resid(pF);
	pV->setZero();
	for (int i = 0; i < iter; i++) {
		int index = 0;
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
int CAssyMatrix::MGCycle(const CAssyVector *pF, CAssyVector *pV, int iter, int alpha1, int alpha2) const
{
	relax(pF, pV, alpha1);
	if (mpCoarseMatrix != 0) {
		CAssyVector w(mpAssyModel);
		CAssyVector fc(mpCoarseMatrix->mpAssyModel);
		CAssyVector vc(mpCoarseMatrix->mpAssyModel);
		residual(pV, pF, &w);
		w.restrictTo(&fc);
		vc.setZero();
		for (int i = 0; i < iter; i++) {
			mpCoarseMatrix->MGCycle(&vc, &fc, iter, alpha1, alpha2);
		}
		w.prolongateFrom(&vc);
		pV->add(&w);
	} else {
	}
	relax(pF, pV, alpha2);
	return 1;
}
int CAssyMatrix::MGInitialGuess(const CAssyVector *pF, CAssyVector *pV) const
{
	printf(" enter CAssyMatrix::MGInitialGuess %d \n", mMGLevel);
	int mAlpha; 
	if (mpCoarseMatrix != 0) {
		CAssyVector w(mpAssyModel);
		CAssyVector fc(mpCoarseMatrix->mpAssyModel);
		CAssyVector vc(mpCoarseMatrix->mpAssyModel);
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
	printf(" exit CAssyMatrix::MGInitialGuess %d \n", mMGLevel);
	return 1;
}
}
