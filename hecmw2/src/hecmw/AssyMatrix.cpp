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

CAssyMatrix::CAssyMatrix(CAssyModel *pAssyModel, const uint& nDOF)//, const vuint& vPartsID)
{
	mpAssyModel = pAssyModel;
        mnDOF = nDOF;//アセンブル方程式のDOF

        //printf("mnDOF %d \n", mnDOF);//'10.12.27 追加
	
	for (int i = 0; i < pAssyModel->getNumOfMesh(); i++) {
		CMesh *pMesh = pAssyModel->getMesh(i);
		CMatrixBCRS *pMatrix = new CMatrixBCRS(pMesh, mnDOF);
		mvMatrix.push_back( pMatrix );
	}

	mMGLevel = pAssyModel->getMGLevel();
	//printf("mMGLevel %d \n",mMGLevel);

	uint numOfContact= pAssyModel->getNumOfContactMesh();/////////////////////////
	//printf("numOfContact %d \n",numOfContact);
	
        uint idof;
        //DOF数ぶんの関係式ポインターを準備
        CEquation **vEquation;
        vEquation = new CEquation*[mnDOF];
        
        

	for(uint icont = 0; icont < numOfContact; icont++){
		CContactMesh* pConMesh= pAssyModel->getContactMesh(icont);
		uint numOfSPoint = pConMesh->getNumOfSlavePoint();
		//printf("numOfSPoint, icont %d %d \n", numOfSPoint, icont);
		
		CMPCMatrix* mpc = new CMPCMatrix();

		for(uint islave = 0; islave< numOfSPoint; islave++){
			CContactNode* pSlaveNode = pConMesh->getSlaveConNode(islave);
			uint mgLevel = mMGLevel;
			if(pSlaveNode->have_MasterFaceID(mgLevel)){
                                ////////////////// DOF数ぶんのMPC関係式 生成 ///////////////////////
                                for(idof=0; idof < mnDOF; idof++) vEquation[idof]= new CEquation();

				uint masterFaceID = pSlaveNode->getMasterFaceID(mgLevel);
				uint smesh = pSlaveNode->getMeshID(); //pConMesh->getSlaveMeshID(islave);
				uint snode = pSlaveNode->getNodeID();
				//printf("masterFaceID, islave, smesh, node %d %d %d %d\n", masterFaceID, islave, smesh, snode);
                                ///////////////////////////////////////////////////////////////////////////////
                                for(idof=0; idof < mnDOF; idof++) vEquation[idof]->setSlave(smesh, snode, idof);
				
				CSkinFace* pMasterFace = pConMesh->getMasterFace_ID(masterFaceID);
				uint numOfVert = pMasterFace->getNumOfNode();
				//printf("numOfVert %d \n",numOfVert);
				
				for(uint ivert=0; ivert< numOfVert; ivert++){
					int mmesh = pMasterFace->getMeshID(); //pConMesh->getMasterMeshID(islave);
					int node = ( pMasterFace->getNode(ivert) )->getNodeID();
					double coef = pMasterFace->getCoef(pSlaveNode->getID(),ivert);
                                        /////////////////////////////////////////////////////////////////////////////////////
                                        for(idof=0; idof < mnDOF; idof++) vEquation[idof]->addMaster(mmesh, node, idof, coef);

					//printf("ivert, mmesh, node, coef %d %d %d %e \n", ivert, mmesh, node, coef);
				};
                                ///////////////////////////////////////////////////////////////////
                                for(idof=0; idof < mnDOF; idof++) mpc->addEquation(vEquation[idof]);
			};
		};
                //////////////////////////
		mvMPCMatrix.push_back(mpc);

		//printf("mvMPCMatrix.size() %d \n", (uint)mvMPCMatrix.size());
	};
}

CAssyMatrix::~CAssyMatrix()
{
	for (CVMatIter im = mvMatrix.begin(); im != mvMatrix.end(); im++) {
		delete *im;
	};
}

int CAssyMatrix::Matrix_Add_Nodal(const uint& iMesh, const uint& iNodeID, const uint& jNodeID, double* NodalMatrix)
{
    mvMatrix[iMesh]->Matrix_Add_Nodal(iNodeID, jNodeID, NodalMatrix);
    return 1;
}

int CAssyMatrix::Matrix_Add_Elem(CAssyModel *pAssyModel, const uint& iMesh, const uint& iElem, double *ElemMatrix)
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
void CAssyMatrix::Matrix_Clear(const uint& iMesh)
{
    mvMatrix[iMesh]->Matrix_Clear();
}

// 対角項=val, mvD[inode]の非対角項=0
//
void CAssyMatrix::setValue(const uint& imesh, const uint& inode, const uint& idof, const double& value)
{
#ifdef ADVANCESOFT_DEBUG
    printf(" enter CAssyMatrix::setValue %d %d %d %e \n", imesh, inode, idof, value);
#endif

    mvMatrix[imesh]->setValue(inode, idof, value);

#ifdef ADVANCESOFT_DEBUG
    printf(" exit CAssyMatrix::setValue %d %d %d %e \n", imesh, inode, idof, value);
#endif
}

// 対角項=val
//
void CAssyMatrix::setValue_D(const uint& imesh, const uint& inode, const uint& idof, const double& value)
{
    mvMatrix[imesh]->setValue_D(inode, idof, value);
}

void CAssyMatrix::setZero_NonDiag(const uint& imesh, const uint& inode, const uint& idof)
{
    mvMatrix[imesh]->setZero_NonDiag(inode, idof);
}


// xに値を入れて、Ax=B でBを求める
void CAssyMatrix::multVector(const uint& imesh, CAssyVector* pX, CAssyVector* pB)
{
    // calc: B = A x
    mvMatrix[imesh]->multVector(pX->getVector(imesh), pB->getVector(imesh));
}

int CAssyMatrix::multVector(CAssyVector *pV, CAssyVector *pP, CAssyVector *pW) const
{
	pV->updateCommBoundary(); // for parallel processing

	if (mvMPCMatrix.size() == 0) {
		// p = (A + Ac) v
		// 1: calc p = A v
		for (int i = 0; i < mvMatrix.size(); i++) {
			mvMatrix[i]->multVector(pV->getVector(i), pP->getVector(i));
		}
		// 1a: calc p += Ac v
		if (mvContactMatrix.size() > 0) {
		}
	} else {
		// p = T (A + Ac) T' v
		// prepare w
		if (pW == 0) {
			pW = new CAssyVector(pP);
		}
		// 1: calc p = T'v
		for (int i = 0; i < mvMPCMatrix.size(); i++) {
			mvMPCMatrix[i]->multVector(pV, pP);
		}
		// 2: calc w = A p
		for (int i = 0; i < mvMatrix.size(); i++) {
			mvMatrix[i]->multVector(pP->getVector(i), pW->getVector(i));
		}
		// 2a: calc w += Ac p
		if (mvContactMatrix.size() > 0) {
			for (int i = 0; i < mvContactMatrix.size(); i++) {
//				mvContactMatrix[i]->multVectorAdd(pP, pW);
			}
		}
		// 3: calc p = T w
		for (int i = 0; i < mvMPCMatrix.size(); i++) {
			mvMPCMatrix[i]->transMultVector(pW, pP);
		}
	}

	pP->sumupCommBoundary(); // for parallel processing
	
	return 1;
}

int CAssyMatrix::multMPC(CAssyVector *pV, CAssyVector *pP) const
{
	if (mvMPCMatrix.size() == 0) {
	} else {
		// 1: calc p = T'v
		for (int i = 0; i < mvMPCMatrix.size(); i++) {
			mvMPCMatrix[i]->multVector(pV, pP);
		}
	}
	return 1;
}

int CAssyMatrix::residual(CAssyVector *pV, const CAssyVector *pF, CAssyVector *pR) const
{
	// 1: r = A v
	multVector(pV, pR);
	// 2: r = f - r
	for (int i = 0; i < pR->size(); i++) {
		(*pR)[i] = (*pF)[i] - (*pR)[i];
	}
	return 1;
}

int CAssyMatrix::setupSolver(int type)
{
	// TODO: implement setupSolver
	// type
	// iter
	// tolerance
	// others
	return 1;
}

int CAssyMatrix::setupPreconditioner(int type) const
{
	// TODO: implement setupPreconditioner
	// type
	// iter
	// level of fill-in???
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
	// TODO: implement setupSmoother
	// type
	// iter???
	for (CVMatIter im = mvMatrix.begin(); im != mvMatrix.end(); im++) {
		(*im)->setupSmoother(type);
	}
	return 1;
}

// pF:右辺ベクトル, pV:解ベクトル
//
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
		//   for each matrix
		int index = 0;
		for (CVMatConstIter im = mvMatrix.begin(); im != mvMatrix.end(); im++) {
			//     call preconditioner
			(*im)->precond(pR->getVector(index), z_resid.getVector(index));
			index++;
		}

		//   accumulate
		pZ->add(&z_resid);

		if (i == iter - 1) break;

		// z_resid = r - A z
		residual(pZ, pR, &z_resid);
	}
	return 1;
}

int CAssyMatrix::relax(const CAssyVector *pF, CAssyVector *pV, int iter) const
{
	CAssyVector v_resid(pF);
	pV->setZero();
	for (int i = 0; i < iter; i++) {
		//   for each matrix
		int index = 0;
		for (CVMatConstIter im = mvMatrix.begin(); im != mvMatrix.end(); im++) {
			//     call preconditioner
			(*im)->relax(pF->getVector(index), v_resid.getVector(index));
			index++;
		}

		//   accumulate
		pV->add(&v_resid);

		if (i == iter - 1) break;

		// v_resid = f - A v
		residual(pV, pF, &v_resid);
	}
	return 1;
}

int CAssyMatrix::MGCycle(const CAssyVector *pF, CAssyVector *pV, int iter, int alpha1, int alpha2) const
{
	relax(pF, pV, alpha1);
	if (mpCoarseMatrix != 0) {
		CAssyVector w(mpAssyModel, mnDOF);
		CAssyVector fc(mpCoarseMatrix->mpAssyModel, mnDOF);
		CAssyVector vc(mpCoarseMatrix->mpAssyModel, mnDOF);

		residual(pV, pF, &w);
		w.restrictTo(&fc);
		vc.setZero();

		for (int i = 0; i < iter; i++) {
                    mpCoarseMatrix->MGCycle(&vc, &fc, iter, alpha1, alpha2);
		}

		w.prolongateFrom(&vc);
		pV->add(&w);
	} else {
		// solve(pF, pV, iter_solve);
	}
	relax(pF, pV, alpha2);
	return 1;
}

int CAssyMatrix::MGInitialGuess(const CAssyVector *pF, CAssyVector *pV) const
{
	printf(" enter CAssyMatrix::MGInitialGuess %d \n", mMGLevel);
	int mAlpha; // TODO: TEMPORARY!!!

	//if (mpCoarseMatrix != 0) {
        if (mpCoarseMatrix != NULL) {
		CAssyVector w(mpAssyModel, mnDOF);
		CAssyVector fc(mpCoarseMatrix->mpAssyModel, mnDOF);
		CAssyVector vc(mpCoarseMatrix->mpAssyModel, mnDOF);

		residual(pV, pF, &w);
		w.restrictTo(&fc);
		vc.setZero();
		mpCoarseMatrix->MGInitialGuess(&fc, &vc);
                
                cout << "AssyMatrix::MGInitialGuess, mgLevel = " << mMGLevel << endl;

		w.prolongateFrom(&vc);
		pV->add(&w);
	} else {
		//   solve(pF, pV, iter_solve);
		precond(pF, pV, 1);
	}
        //////	relax(pF, pV, mAlpha);
	precond(pF, pV, 1);
	printf(" exit CAssyMatrix::MGInitialGuess %d \n", mMGLevel);
	return 1;
}

// デバッグ
// ------
// 2011.01.05
//
void CAssyMatrix::dump()
{
    uint i;
    for(i=0; i < mvMatrix.size(); i++){
        mvMatrix[i]->dump();
    };
    for(i=0; i < mvMPCMatrix.size(); i++){
        mvMPCMatrix[i]->dump();
    };
}

}//namespace pmw;
