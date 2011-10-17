/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/MatrixBCRS.cpp
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
#include "MatrixBCRS.h"
#include "Mesh.h"
#include "Node.h"
#include "Element.h"
#include "Vector.h"
#include <algorithm>
namespace pmw
{
CMatrixBCRS::CMatrixBCRS(CMesh *pMesh, const uiint& nDOF)
{
#ifdef ADVANCESOFT_DEBUG
	printf("enter CMatrixBCRS::CMatrixBCRS \n");
#endif
        mnDOF = nDOF;
	mnNode = pMesh->getNumOfNode();
        CIndexBucket *pBucket= pMesh->getBucket();
	mvIndexL.resize(mnNode+1); 
	mvIndexU.resize(mnNode+1); 
	mvIndexL[0] = 0;
	mvIndexU[0] = 0;
	for (uiint i_node = 0; i_node < mnNode; i_node++) {
            std::vector<uiint> v_item_l;
            std::vector<uiint> v_item_u;
            CNode *pNode= pMesh->getNodeIX(i_node);
            uiint i_node_id = pNode->getID();
            CElement *pElement;
            CAggregateElement *pAggElement= pMesh->getAggElem(i_node_id);
            uiint nNumOfElement= pAggElement->getNumOfElement();
            uiint i_elem;
            for(i_elem=0; i_elem < nNumOfElement; i_elem++){
                pElement= pAggElement->get(i_elem);
                uiint k_node, nNumOfNode= pElement->getNumOfNode();
                for (k_node = 0; k_node < nNumOfNode; k_node++) {
                    uiint k_node_id = pElement->getNode(k_node)->getID();
                    uiint k_index = pBucket->getIndexNode(k_node_id);
                    if(k_node_id < i_node_id){
                        v_item_l.push_back(k_index);
                    }else if(i_node_id < k_node_id){
                        v_item_u.push_back(k_index);
                    }
                };
            };
            std::sort(v_item_l.begin(), v_item_l.end());
            std::vector<uiint>::iterator new_end = std::unique(v_item_l.begin(), v_item_l.end());
            v_item_l.erase(new_end, v_item_l.end());
            mvIndexL[i_node + 1] = mvIndexL[i_node] + v_item_l.size();
            mvItemL.insert(mvItemL.end(), v_item_l.begin(), v_item_l.end());
            std::sort(v_item_u.begin(), v_item_u.end());
            new_end = std::unique(v_item_u.begin(), v_item_u.end());
            v_item_u.erase(new_end, v_item_u.end());
            mvIndexU[i_node + 1] = mvIndexU[i_node] + v_item_u.size();
            mvItemU.insert(mvItemU.end(), v_item_u.begin(), v_item_u.end());
	}
	mvD.resize(mnNode);
	for (uiint i = 0; i < mnNode; i++) {
		mvD[i].resize(mnDOF, mnDOF);
		for(uiint i1=0; i1<mnDOF; i1++) for(uiint i2=0; i2<mnDOF; i2++) mvD[i](i1,i2)=0.0;
	}
	mvALU.resize(mnNode);
	for (uiint i = 0; i < mnNode; i++) {
		mvALU[i].resize(mnDOF, mnDOF);
		for(uiint i1=0; i1<mnDOF; i1++) for(uiint i2=0; i2<mnDOF; i2++) mvALU[i](i1,i2)=0.0;
	}
	mINL = mvIndexL[mnNode];
	mvAL.resize(mINL);
	for (uiint i = 0; i < mINL; i++) {
		mvAL[i].resize(mnDOF, mnDOF);
		for(uiint i1=0; i1<mnDOF; i1++) for(uiint i2=0; i2<mnDOF; i2++) mvAL[i](i1,i2)=0.0;
	}
	mINU = mvIndexU[mnNode];
	mvAU.resize(mINU);
	for (uiint i = 0; i < mINU; i++) {
		mvAU[i].resize(mnDOF, mnDOF);
		for(uiint i1=0; i1<mnDOF; i1++) for(uiint i2=0; i2<mnDOF; i2++) mvAU[i](i1,i2)=0.0;
	}
#ifdef ADVANCESOFT_DEBUG	
	for (uiint i = 0; i < mnNode; i++) {
		printf(" %ld %ld %ld ; ", i, mvIndexL[i], mvIndexL[i+1]-1);
		for (uiint j = mvIndexL[i]; j < mvIndexL[i+1]; j++) {
			printf(" %ld",mvItemL[j]);
		}
		cout << endl;
	}
	for (uiint i = 0; i < mnNode; i++) {
		printf(" %ld %ld %ld ; ", i, mvIndexU[i], mvIndexU[i+1]-1);
		for (uiint j = mvIndexU[i]; j < mvIndexU[i+1]; j++) {
			printf(" %ld",mvItemU[j]);
		}
		cout << endl;
	}
	printf("exit CMatrixBCRS::CMatrixBCRS \n");
#endif
}
CMatrixBCRS::~CMatrixBCRS()
{
}
uiint CMatrixBCRS::Matrix_Add_Nodal(const uiint& inode, const uiint& jnode, const double* NodalMatrix)
{
    uiint nMatSize = mnDOF;
    uiint kL, kU;
    bool bL,bU;
    uiint irow = inode;
    uiint icol = jnode;
    if( icol == irow ) {
        for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
            mvD[icol](ii,jj) += NodalMatrix[nMatSize*ii+jj];
        }
    } else if( icol < irow ) {
        bL = false;
        for (uiint k = mvIndexL[irow]; k < mvIndexL[irow+1]; k++) {
            if( icol == mvItemL[k] ) {
                kL = k;
                bL = true;
            }
        }
        if( !bL ) printf("***** error in matrix index ***** %ld %ld %s \n",irow,icol,"-1");
        for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
            mvAL[kL](ii,jj) += NodalMatrix[nMatSize*ii+jj];
        }
    } else if( icol > irow ) {
        bU = false;
        for (uiint k = mvIndexU[irow]; k < mvIndexU[irow+1]; k++) {
            if( icol == mvItemU[k] ) {
                kU = k;
                bU = true;
            }
        }
        if( !bU ) printf("***** error in matrix index ***** %ld %ld %s \n",irow,icol,"-1");
        for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
            mvAU[kU](ii,jj) += NodalMatrix[nMatSize*ii+jj];
        }
    }
    return 1;
}
uiint CMatrixBCRS::Matrix_Add_Elem(CMesh *pMesh, const uiint& iElem, double *ElemMatrix)
{
#ifdef ADVANCESOFT_DEBUG
    printf(" enter CMatrixBCRS::Matrix_Add_Elem %ld %e \n", iElem, ElemMatrix[0]);
#endif
    CElement *pElement = pMesh->getElementIX(iElem);
    vector<CNode*> vNode= pElement->getNode();
    uiint nLocalNode = vNode.size();
    uiint nMatSize = nLocalNode * mnDOF;
    CIndexBucket *pBucket = pMesh->getBucket();
    for(uiint i=0; i< nLocalNode; i++) for(uiint j=0; j< nLocalNode; j++){
        uiint kL, kU;
        bool bL, bU;
        uiint iNodeID = vNode[i]->getID(); 
        uiint jNodeID = vNode[j]->getID();
        uiint irow = pBucket->getIndexNode(iNodeID);  
        uiint icol = pBucket->getIndexNode(jNodeID);  
        if( icol == irow ) {
            for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
                mvD[icol](ii,jj) += ElemMatrix[nMatSize*(ii+i*mnDOF) + jj+j*mnDOF];
            }
        } else if( icol < irow ) {
            bL = false;
            for (uiint k = mvIndexL[irow]; k < mvIndexL[irow+1]; k++) {
                if( icol == mvItemL[k] ) { 
                    kL = k;
                    bL = true;
                }
            }
            if( !bL ) printf("***** error in matrix index ***** %ld %ld %ld %ld %s \n",i,j,irow,icol,"-1");
            for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
                mvAL[kL](ii,jj) += ElemMatrix[nMatSize*(ii+i*mnDOF)+jj+j*mnDOF];
            }
        } else if( icol > irow ) {
            bU = false;
            for (uiint k = mvIndexU[irow]; k < mvIndexU[irow+1]; k++) {
                if( icol == mvItemU[k] ) { 
                    kU = k;
                    bU = true;
                }
            }
            if( !bU ) printf("***** error in matrix index ***** %ld %ld %ld %ld %s \n",i,j,irow,icol,"-1");
            for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
                mvAU[kU](ii,jj) += ElemMatrix[nMatSize*(ii+i*mnDOF) + jj+j*mnDOF];
            }
        }
    }
#ifdef ADVANCESOFT_DEBUG
    printf(" exit CMatrixBCRS::Matrix_Add_Elem \n");
#endif
    return 1;
}
void CMatrixBCRS::Matrix_Clear()
{
    for(uiint i=0; i< mnNode; i++){
        for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
            mvD[i](ii,jj) = 0.0;
        };
    };
    for(uiint i=0; i < mINL; i++){
        for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
            mvAL[i](ii,jj) = 0.0;
        };
    };
    for(uiint i=0; i < mINU; i++){
        for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
            mvAU[i](ii,jj) = 0.0;
        };
    };
}
void CMatrixBCRS::multVector(CVector *pV, CVector *pP) const
{
    for (uiint i = 0; i < mnNode; i++) {
        (*pP)[i] = prod(mvD[i], (*pV)[i]);
        for (uiint j = mvIndexL[i]; j < mvIndexL[i+1]; j++) {
            (*pP)[i] += prod(mvAL[j], (*pV)[mvItemL[j]]);
        };
        for (uiint j = mvIndexU[i]; j < mvIndexU[i+1]; j++) {
            (*pP)[i] += prod(mvAU[j], (*pV)[mvItemU[j]]);
        };
    };
}
void CMatrixBCRS::setValue_D(const uiint& inode, const uiint& idof, const double& value)
{
    mvD[inode](idof, idof) = value;
}
void CMatrixBCRS::setValue(const uiint& inode, const uiint& idof, const double& dDiag, CVector *pRHS, const double& dRHS)
{
    for(uiint j=0; j < mnDOF; j++){
        if(idof != j){
            mvD[inode](idof, j) = 0.0;
        }
    };
    for(uiint j=0; j < mnDOF; j++){
        if(idof != j){
            double diffVal = -mvD[inode](j, idof) * dRHS;
            pRHS->addValue(inode, j, diffVal);
            mvD[inode](j, idof) = 0.0;
        }
    };
    mvD[inode](idof,idof) = dDiag;
}
void CMatrixBCRS::setZero_NonDiag(const uiint& inode, const uiint& idof, CVector* pRHS, const double& dRHS)
{
    for(uiint j=mvIndexL[inode]; j < mvIndexL[inode+1]; j++){
        for(uiint jdof=0; jdof < mnDOF; jdof++){
            mvAL[j](idof,jdof) = 0.0;
        };
    };
    for(uiint j=mvIndexU[inode]; j < mvIndexU[inode+1]; j++){
        for(uiint jdof=0; jdof < mnDOF; jdof++){
            mvAU[j](idof,jdof) = 0.0;
        };
    };
    for(uiint jnode=0; jnode < mnNode; jnode++){
    for(uiint i=mvIndexL[jnode]; i < mvIndexL[jnode+1]; i++){
        uiint col=mvItemL[i];
        if(col==inode){
        for(uiint jdof=0; jdof < mnDOF; jdof++){
            double diffVal = -mvAL[i](jdof,idof) * dRHS;
            pRHS->addValue(jnode, jdof, diffVal);
            mvAL[i](jdof,idof) = 0.0;
        };
        }
    };
    };
    for(uiint jnode=0; jnode < mnNode; jnode++){
    for(uiint i=mvIndexU[jnode]; i < mvIndexU[jnode+1]; i++){
        uiint col=mvItemU[i];
        if(col==inode){
        for(uiint jdof=0; jdof < mnDOF; jdof++){
            double diffVal = -mvAU[i](jdof, idof) * dRHS;
            pRHS->addValue(jnode, jdof, diffVal);
            mvAU[i](jdof,idof) = 0.0;
        };
        }
    };
    };
}
uiint CMatrixBCRS::setupPreconditioner(iint type)
{
#ifdef ADVANCESOFT_DEBUG
  printf(" enter CMatrixBCRS::setupPreconditioner \n");
#endif
  ublas::matrix<double> pA(mnDOF,mnDOF), pB(mnDOF,mnDOF), pC(mnDOF,mnDOF);
  mPrecond = type;
  iint itype;
  itype = -99;
  if( mPrecond == 1 ) itype = 2;
  if( mPrecond == 2 ) itype = 2;
  if( mPrecond == 3 ) itype = 2;
  if( mPrecond == 4 ) itype = 1;
  if( itype == -99 ) printf("setup precondition; fatal error %ld\n", mPrecond);
  switch( itype ){
  case( 1 ):
    printf(" [TYPE:1] Preconditioner \n");
    for (int i = 0; i < mnNode; i++) {
      pA = mvD[i];
      for (uiint j = mvIndexL[i]; j < mvIndexL[i+1]; j++) {
	transpose(mvAL[j], &pB);
	pC  = prod(mvAL[j], mvALU[mvItemL[j]] );
	pA -= prod(pC, pB);
      };
      inverse(pA, &mvALU[i]);
    };
    break;
  case( 4 ):
    printf(" [TYPE:4] Preconditioner \n");
    for (uiint i = 0; i < mnNode; i++) mvALU[i] = mvD[i];
    for (uiint i = 0; i < mnNode; i++) {
      pA = mvALU[i];
      inverse(pA, &mvALU[i]);
      pA = mvALU[i];
      for (uiint j = mvIndexU[i]; j < mvIndexU[i+1]; j++) {
	pC  = prod(pA, mvAU[j]);
	transpose(mvAU[j], &pB);
	mvALU[mvItemU[j]] -= prod(pB, pC);
      };
    };
    break;
  case( 2 ):
    printf(" [TYPE:2] Preconditioner \n");
    for (uiint i = 0; i < mnNode; i++) {
      inverse(mvD[i], &mvALU[i]);
    };
    break;
  case( 3 ):
    printf(" [TYPE:3] Preconditioner \n");
    break;
  }
#ifdef ADVANCESOFT_DEBUG
  printf(" exit CMatrixBCRS::setupPreconditioner \n");
#endif
  return 1;
}
double CMatrixBCRS::inverse(ublas::matrix<double> pA, ublas::matrix<double> *pB)
{
	double det = determinant( pA );
	double Recip = 1.0 / det;
	(*pB)(0, 0) = Recip * ( pA(1, 1) * pA(2, 2) - pA(2, 1) * pA(1, 2) );
	(*pB)(0, 1) = Recip * (-pA(0, 1) * pA(2, 2) + pA(2, 1) * pA(0, 2) );
	(*pB)(0, 2) = Recip * ( pA(0, 1) * pA(1, 2) - pA(1, 1) * pA(0, 2) );
	(*pB)(1, 0) = Recip * (-pA(1, 0) * pA(2, 2) + pA(2, 0) * pA(1, 2) );
	(*pB)(1, 1) = Recip * ( pA(0, 0) * pA(2, 2) - pA(2, 0) * pA(0, 2) );
	(*pB)(1, 2) = Recip * (-pA(0, 0) * pA(1, 2) + pA(1, 0) * pA(0, 2) );
	(*pB)(2, 0) = Recip * ( pA(1, 0) * pA(2, 1) - pA(2, 0) * pA(1, 1) );
	(*pB)(2, 1) = Recip * (-pA(0, 0) * pA(2, 1) + pA(2, 0) * pA(0, 1) );
	(*pB)(2, 2) = Recip * ( pA(0, 0) * pA(1, 1) - pA(1, 0) * pA(0, 1) );
	return det;
}	
void CMatrixBCRS::transpose(ublas::matrix<double> pA, ublas::matrix<double> *pB)
{
	for (int i=0; i< mnDOF; i++) for (int j=0; j<mnDOF; j++) (*pB)(i, j) = pA(j, i);
}	
double CMatrixBCRS::determinant(ublas::matrix<double> pA)
{
	double det = pA(0, 0) * pA(1, 1) * pA(2, 2)
		+ pA(1, 0) * pA(2, 1) * pA(0, 2)
		+ pA(2, 0) * pA(0, 1) * pA(1, 2)
		- pA(2, 0) * pA(1, 1) * pA(0, 2)
		- pA(1, 0) * pA(0, 1) * pA(2, 2)
		- pA(0, 0) * pA(2, 1) * pA(1, 2);
	return det;
}	
uiint CMatrixBCRS::setupSmoother(iint type)
{
   	printf(" enter CMatrixBCRS::setupSmoother \n");
   	printf(" exit CMatrixBCRS::setupSmoother \n");
	return 1;
}
uiint CMatrixBCRS::precond(const CVector *pR, CVector *pZ) const
{
  CVector::ElemType WW(mnDOF);
  ublas::zero_vector<double> vzero(mnDOF);
  iint itype;
  itype = -99;
  if( mPrecond == 1 ) itype = 1;
  if( mPrecond == 2 ) itype = 2;
  if( mPrecond == 3 ) itype = 2;
  if( mPrecond == 4 ) itype = 1;
  if( itype == -99 ) printf("precond; fatal error %ld\n", mPrecond, itype);
  switch( itype ){
  case( 0 ):
    pZ->subst(pR);
    break;
  case( 1 ):
    pZ->subst(pR);
    for(uiint i=0; i< mnNode; i++) {
      for (uiint j = mvIndexL[i]; j < mvIndexL[i+1]; j++) {
	(*pZ)[i] -= prod(mvAL[j], (*pZ)[mvItemL[j]]);
      }
      (*pZ)[i] = prod(mvALU[i], (*pZ)[i]);
    }
    for(uiint i= mnNode-1; i >= 0 && i < mnNode; i--){
      WW(0)=0.0; WW(1)=0.0; WW(2)=0.0;
      for (uiint j = mvIndexU[i]; j < mvIndexU[i+1]; j++) {
	WW += prod(mvAU[j], (*pZ)[mvItemU[j]]);
      }
      (*pZ)[i] -= prod(mvALU[i], WW);
    }
    break;
  case( 2 ):
    for(uiint loop=0; loop < 10; loop++) {
      for(uiint i=0; i< mnNode; i++) {
	WW = (*pR)[i];
	for (uiint j = mvIndexL[i]; j < mvIndexL[i+1]; j++) {
	  WW -= prod(mvAL[j], (*pZ)[mvItemL[j]]);
	}
	for (uiint j = mvIndexU[i]; j < mvIndexU[i+1]; j++) {
	  WW -= prod(mvAU[j], (*pZ)[mvItemU[j]]);
	}
	WW = prod(mvALU[i], WW);
	double omega = 0.5;
	(*pZ)[i] = (*pZ)[i] + omega * ( WW - (*pZ)[i] );
      }
    }
    break;
  case( 3 ):
    for(uiint i=0; i< mnNode; i++) (*pZ)[i] = vzero;
    for(uiint i=0; i< 50; i++) relax(pR, pZ);
    break;
  }
    return 1;
}
uiint CMatrixBCRS::relax(const CVector *pR, CVector *pZ) const
{
  CVector::ElemType WW(mnDOF);
  for(uiint i=0; i< mnNode; i++) {
    WW = (*pR)[i];
    for (uiint j = mvIndexL[i]; j < mvIndexL[i+1]; j++) {
      WW -= prod(mvAL[j], (*pZ)[mvItemL[j]]);
    }
    for (uiint j = mvIndexU[i]; j < mvIndexU[i+1]; j++) {
      WW -= prod(mvAU[j], (*pZ)[mvItemU[j]]);
    }
    WW = prod(mvALU[i], WW);
    double omega = 0.5;
    (*pZ)[i] = (*pZ)[i] + omega * ( WW - (*pZ)[i] );
  }
  return 1;
}
void CMatrixBCRS::dump()
{
    cout << " ---- mvD ---- " << endl;
    for(uiint i=0; i< mnNode; i++){
        cout << "Node i:" << i << " ";
        for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
            cout << mvD[i](ii,jj) << "  ";
        };
        cout << endl;
    };
    cout << " ---- mvAL ---- " << endl;
    for(uiint i=0; i < mINL; i++){
        cout << "mINL:" << i << " ";
        for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
            cout << mvAL[i](ii,jj) << "  ";
        };
        cout << endl;
    };
    cout << " ---- mvAU ---- " << endl;
    for(uiint i=0; i < mINU; i++){
        cout << "mINU:" << i << " ";
        for(uiint ii=0; ii < mnDOF; ii++) for(uiint jj=0; jj < mnDOF; jj++) {
            cout << mvAU[i](ii,jj) << "  ";
        };
        cout << endl;
    };
}
}
