/*
 * MatrixBCRS.cpp
 *
 *  Created on: Oct 16, 2009
 *      Author: goto
 */

#include "MatrixBCRS.h"
#include "Mesh.h"
#include "Node.h"
#include "Element.h"
#include "Vector.h"
#include <algorithm>

namespace pmw
{

//template<int NDOF>
CMatrixBCRS::CMatrixBCRS(CMesh *pMesh, const uint& nDOF)
{
#ifdef ADVANCESOFT_DEBUG
	printf("enter CMatrixBCRS::CMatrixBCRS \n");
#endif
        mnDOF = nDOF;
	mnNode = pMesh->getNumOfNode();
	printf("initialize in CMatrixBCRS %d %d \n",mnDOF,mnNode);
	mvIndexL.resize(mnNode+1); // K.Matsubara
	mvIndexU.resize(mnNode+1); // K.Matsubara

	mvIndexL[0] = 0;
	mvIndexU[0] = 0;

	for (int i_node = 0; i_node < mnNode; i_node++) {
            std::vector<int> v_item_l;
            std::vector<int> v_item_u;

            CNode *pNode= pMesh->getNodeIX(i_node);
            uint i_node_id = pNode->getID();

            CElement *pElement;
            CAggregateElement *pAggElement= pMesh->getAggElem(i_node_id);

            uint numOfElement= pAggElement->getNumOfElement();
            uint i_elem;
            for(i_elem=0; i_elem < numOfElement; i_elem++){
                pElement= pAggElement->get(i_elem);

                uint numOfNode= pElement->getNumOfNode();
                uint k_node;
                for (k_node = 0; k_node < numOfNode; k_node++) {
                    uint k_node_id = pElement->getNode(k_node)->getID();
                    if(k_node_id < i_node_id){
                        v_item_l.push_back(k_node_id);
                    }else if(i_node_id < k_node_id){
                        v_item_u.push_back(k_node_id);
                    }
                };
            };

            std::sort(v_item_l.begin(), v_item_l.end());
            std::vector<int>::iterator new_end = std::unique(v_item_l.begin(), v_item_l.end());
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
	for (int i = 0; i < mnNode; i++) {
		mvD[i].resize(mnDOF, mnDOF);
		for(int i1=0; i1<mnDOF; i1++) for(int i2=0; i2<mnDOF; i2++) mvD[i](i1,i2)=0.0;
	}
	mvALU.resize(mnNode);
	for (int i = 0; i < mnNode; i++) {
		mvALU[i].resize(mnDOF, mnDOF);
		for(int i1=0; i1<mnDOF; i1++) for(int i2=0; i2<mnDOF; i2++) mvALU[i](i1,i2)=0.0;
	}
	mINL = mvIndexL[mnNode];
        cout << "MatrixBCRS::MatrixBCRS mINL=" << mINL << endl;

	mvAL.resize(mINL);
	for (int i = 0; i < mINL; i++) {
		mvAL[i].resize(mnDOF, mnDOF);
		for(int i1=0; i1<mnDOF; i1++) for(int i2=0; i2<mnDOF; i2++) mvAL[i](i1,i2)=0.0;
	}


	mINU = mvIndexU[mnNode];
        cout << "MatrixBCRS::MatrixBCRS mINU=" << mINU << endl;

	mvAU.resize(mINU);
	for (int i = 0; i < mINU; i++) {
		mvAU[i].resize(mnDOF, mnDOF);
		for(int i1=0; i1<mnDOF; i1++) for(int i2=0; i2<mnDOF; i2++) mvAU[i](i1,i2)=0.0;
	}

#ifdef ADVANCESOFT_DEBUG	
	for (int i = 0; i < mnNode; i++) {
		printf(" %d %d %d ; ", i, mvIndexL[i], mvIndexL[i+1]-1);
		for (int j = mvIndexL[i]; j < mvIndexL[i+1]; j++) {
			printf(" %d",mvItemL[j]);
		}
		cout << endl;
	}
	for (int i = 0; i < mnNode; i++) {
		printf(" %d %d %d ; ", i, mvIndexU[i], mvIndexU[i+1]-1);
		for (int j = mvIndexU[i]; j < mvIndexU[i+1]; j++) {
			printf(" %d",mvItemU[j]);
		}
		cout << endl;
	}
	printf("exit CMatrixBCRS::CMatrixBCRS \n");
#endif
}

//template<int NDOF>
CMatrixBCRS::~CMatrixBCRS()
{
    // TODO Auto-generated destructor stub
}

//
// add NodalStiff_Matrix 
//
int CMatrixBCRS::Matrix_Add_Nodal(const uint& iNodeID, const uint& jNodeID, const double* NodalMatrix)
{
    uint nMatSize = mnDOF;
    int kL, kU;
    uint irow = iNodeID;
    uint icol = jNodeID;
    if( icol == irow ) {//対角項
        
        for(int ii=0; ii < mnDOF; ii++) for(int jj=0; jj < mnDOF; jj++) {
            mvD[icol](ii,jj) += NodalMatrix[nMatSize*ii+jj];
        }
        
    } else if( icol < irow ) {// 下三角
        kL = -1;
        for (int k = mvIndexL[irow]; k < mvIndexL[irow+1]; k++) {
            if( icol == mvItemL[k] ) {
                kL = k;
            }
        }
        if( kL < 0 ) printf("***** error in matrix index ***** %d %d %d \n",irow,icol,kL);
        
        for(int ii=0; ii < mnDOF; ii++) for(int jj=0; jj < mnDOF; jj++) {
            mvAL[kL](ii,jj) += NodalMatrix[nMatSize*ii+jj];
        }

    } else if( icol > irow ) {// 上三角
        kU = -1;
        for (int k = mvIndexU[irow]; k < mvIndexU[irow+1]; k++) {
            if( icol == mvItemU[k] ) {
                kU = k;
            }
        }
        if( kU < 0 ) printf("***** error in matrix index ***** %d %d %d \n",irow,icol,kU);
        
        for(int ii=0; ii < mnDOF; ii++) for(int jj=0; jj < mnDOF; jj++) {
            mvAU[kU](ii,jj) += NodalMatrix[nMatSize*ii+jj];
        }
    }
    
    return 1;
}
//
// add ElementStiff_Matrix
//
int CMatrixBCRS::Matrix_Add_Elem(CMesh *pMesh, const uint& iElem, double *ElemMatrix)
{
#ifdef ADVANCESOFT_DEBUG
   	printf(" enter CMatrixBCRS::Matrix_Add_Elem %d %e \n", iElem, ElemMatrix[0]);
#endif

    CElement *pElement = pMesh->getElementIX(iElem);
    vector<CNode*> vNode= pElement->getNode();
    
    int nLocalNode = vNode.size();
    int nMatSize = nLocalNode * mnDOF;

    for(int i=0; i< nLocalNode; i++) for(int j=0; j< nLocalNode; j++){
		int kL, kU;
		uint irow = vNode[i]->getID();
		uint icol = vNode[j]->getID();
		if( icol == irow ) {
			for(int ii=0; ii < mnDOF; ii++) for(int jj=0; jj < mnDOF; jj++) {
			    mvD[icol](ii,jj) += ElemMatrix[nMatSize*(ii+i*mnDOF)+jj+j*mnDOF];
			}
		} else if( icol < irow ) {
			kL = -1;
			for (int k = mvIndexL[irow]; k < mvIndexL[irow+1]; k++) {
				if( icol == mvItemL[k] ) {
					kL = k;
				}
			}
			if( kL < 0 ) printf("***** error in matrix index ***** %d %d %d %d %d \n",i,j,irow,icol,kL);
			for(int ii=0; ii < mnDOF; ii++) for(int jj=0; jj < mnDOF; jj++) {
			    mvAL[kL](ii,jj) += ElemMatrix[nMatSize*(ii+i*mnDOF)+jj+j*mnDOF];
			}
		} else if( icol > irow ) {
			kU = -1;
			for (int k = mvIndexU[irow]; k < mvIndexU[irow+1]; k++) {
				if( icol == mvItemU[k] ) {
					kU = k;
				}
			}
			if( kU < 0 ) printf("***** error in matrix index ***** %d %d %d %d %d \n",i,j,irow,icol,kU);
			for(int ii=0; ii < mnDOF; ii++) for(int jj=0; jj < mnDOF; jj++) {
			    mvAU[kU](ii,jj) += ElemMatrix[nMatSize*(ii+i*mnDOF)+jj+j*mnDOF];
			}
		}
    }
#ifdef ADVANCESOFT_DEBUG
   	printf(" exit CMatrixBCRS::Matrix_Add_Elem \n");
#endif

    return 1;
}

// Matrix all 0 clear
//
void CMatrixBCRS::Matrix_Clear()
{
    for(uint i=0; i< mnNode; i++){
        for(int ii=0; ii < mnDOF; ii++) for(int jj=0; jj < mnDOF; jj++) {
            mvD[i](ii,jj) = 0.0;
        };
    };
    for(uint i=0; i < mINL; i++){
        for(int ii=0; ii < mnDOF; ii++) for(int jj=0; jj < mnDOF; jj++) {
            mvAL[i](ii,jj) = 0.0;
        };
    };
    for(uint i=0; i < mINU; i++){
        for(int ii=0; ii < mnDOF; ii++) for(int jj=0; jj < mnDOF; jj++) {
            mvAU[i](ii,jj) = 0.0;
        };
    };
}

void CMatrixBCRS::multVector(CVector *pV, CVector *pP) const
{
	for (int i = 0; i < mnNode; i++) {
		(*pP)[i] = prod(mvD[i], (*pV)[i]);
		for (int j = mvIndexL[i]; j < mvIndexL[i+1]; j++) {
			(*pP)[i] += prod(mvAL[j], (*pV)[mvItemL[j]]);
		}
		for (int j = mvIndexU[i]; j < mvIndexU[i+1]; j++) {
			(*pP)[i] += prod(mvAU[j], (*pV)[mvItemU[j]]);
		}
	}
}

void CMatrixBCRS::setValue(int inode, int idof, double value)
{
#ifdef ADVANCESOFT_DEBUG
   	printf(" enter CMatrixBCRS::setValue \n");
#endif
	mvD[inode](idof,idof) = value;
#ifdef ADVANCESOFT_DEBUG
   	printf(" exit CMatrixBCRS::setValue \n");
#endif
}

int CMatrixBCRS::setupPreconditioner(int type)
{
#ifdef ADVANCESOFT_DEBUG
  printf(" enter CMatrixBCRS::setupPreconditioner \n");
#endif
  ublas::matrix<double> pA(mnDOF,mnDOF), pB(mnDOF,mnDOF), pC(mnDOF,mnDOF);

  mPrecond = type;
  int itype;

  itype = -99;
  if( mPrecond == 1 ) itype = 2;
  if( mPrecond == 2 ) itype = 2;
  if( mPrecond == 3 ) itype = 2;
  if( mPrecond == 4 ) itype = 1;
  if( itype == -99 ) printf("setup precondition; fatal error %d\n", mPrecond);

  switch( itype ){
  case( 1 ):
    printf(" [TYPE:1] Preconditioner \n");
    // D - L * D**(-1) * LT
    for (int i = 0; i < mnNode; i++) {
      pA = mvD[i];
      for (int j = mvIndexL[i]; j < mvIndexL[i+1]; j++) {
	transpose(mvAL[j], &pB);
	pC  = prod(mvAL[j], mvALU[mvItemL[j]] );
	pA -= prod(pC, pB);
      };
      inverse(pA, &mvALU[i]);
    };
    break;
  case( 4 ):
    // D - UT * D**(-1) * U
    printf(" [TYPE:4] Preconditioner \n");
    for (int i = 0; i < mnNode; i++) mvALU[i] = mvD[i];
    for (int i = 0; i < mnNode; i++) {
      pA = mvALU[i];
      inverse(pA, &mvALU[i]);
      pA = mvALU[i];
      for (int j = mvIndexU[i]; j < mvIndexU[i+1]; j++) {
	pC  = prod(pA, mvAU[j]);
	transpose(mvAU[j], &pB);
	mvALU[mvItemU[j]] -= prod(pB, pC);
      };
    };
    break;
  case( 2 ):
    printf(" [TYPE:2] Preconditioner \n");
    for (int i = 0; i < mnNode; i++) {
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

int CMatrixBCRS::setupSmoother(int type)
{
	//
   	printf(" enter CMatrixBCRS::setupSmoother \n");
   	
   	printf(" exit CMatrixBCRS::setupSmoother \n");
	return 1;
}

int CMatrixBCRS::precond(const CVector *pR, CVector *pZ) const
{
  CVector::ElemType WW(mnDOF);
  ublas::zero_vector<double> vzero(mnDOF);

  int itype;
  itype = -99;
  if( mPrecond == 1 ) itype = 1;
  if( mPrecond == 2 ) itype = 2;
  if( mPrecond == 3 ) itype = 2;
  if( mPrecond == 4 ) itype = 1;
  if( itype == -99 ) printf("precond; fatal error %d\n", mPrecond, itype);

  switch( itype ){
  case( 0 ):
    pZ->subst(pR);
    break;
  case( 1 ):
    pZ->subst(pR);
    for(int i=0; i< mnNode; i++) {
      for (int j = mvIndexL[i]; j < mvIndexL[i+1]; j++) {
	(*pZ)[i] -= prod(mvAL[j], (*pZ)[mvItemL[j]]);
      }
      (*pZ)[i] = prod(mvALU[i], (*pZ)[i]);
    }
    for(int i = mnNode-1; i>-1 ; i--) {
      WW(0)=0.0; WW(1)=0.0; WW(2)=0.0;
      for (int j = mvIndexU[i]; j < mvIndexU[i+1]; j++) {
	WW += prod(mvAU[j], (*pZ)[mvItemU[j]]);
      }
      (*pZ)[i] -= prod(mvALU[i], WW);
    }
    break;
  case( 2 ):
    for(int loop=0; loop<10; loop++) {
      for(int i=0; i< mnNode; i++) {
	WW = (*pR)[i];
	for (int j = mvIndexL[i]; j < mvIndexL[i+1]; j++) {
	  WW -= prod(mvAL[j], (*pZ)[mvItemL[j]]);
	}
	for (int j = mvIndexU[i]; j < mvIndexU[i+1]; j++) {
	  WW -= prod(mvAU[j], (*pZ)[mvItemU[j]]);
	}
	WW = prod(mvALU[i], WW);
	double omega = 0.5;
	//					double omega = 1.0;
	(*pZ)[i] = (*pZ)[i] + omega * ( WW - (*pZ)[i] );
      }
    }
    break;
  case( 3 ):
    for(int i=0; i< mnNode; i++) (*pZ)[i] = vzero;
    for(int i=0; i<50; i++) relax(pR, pZ);
    break;
  }
	return 1;
}

int CMatrixBCRS::relax(const CVector *pR, CVector *pZ) const
{
  // call smoother (once)
  CVector::ElemType WW(mnDOF);
  for(int i=0; i< mnNode; i++) {
    WW = (*pR)[i];
    for (int j = mvIndexL[i]; j < mvIndexL[i+1]; j++) {
      WW -= prod(mvAL[j], (*pZ)[mvItemL[j]]);
    }
    for (int j = mvIndexU[i]; j < mvIndexU[i+1]; j++) {
      WW -= prod(mvAU[j], (*pZ)[mvItemU[j]]);
    }
    WW = prod(mvALU[i], WW);
    double omega = 0.5;
    (*pZ)[i] = (*pZ)[i] + omega * ( WW - (*pZ)[i] );
  }
  return 1;
}

// デバッグ
// ------
// 2011.01.05
//
void CMatrixBCRS::dump()
{
    cout << " ---- mvD ---- " << endl;
    for(uint i=0; i< mnNode; i++){
        cout << "Node i:" << i << " ";

        for(int ii=0; ii < mnDOF; ii++) for(int jj=0; jj < mnDOF; jj++) {
            cout << mvD[i](ii,jj) << "  ";
        };
        
        cout << endl;
    };

    cout << " ---- mvAL ---- " << endl;
    for(uint i=0; i < mINL; i++){
        cout << "mINL:" << i << " ";

        for(int ii=0; ii < mnDOF; ii++) for(int jj=0; jj < mnDOF; jj++) {
            cout << mvAL[i](ii,jj) << "  ";
        };
        
        cout << endl;
    };

    cout << " ---- mvAU ---- " << endl;
    for(uint i=0; i < mINU; i++){
        cout << "mINU:" << i << " ";

        for(int ii=0; ii < mnDOF; ii++) for(int jj=0; jj < mnDOF; jj++) {
            cout << mvAU[i](ii,jj) << "  ";
        };
        
        cout << endl;
    };
    
}

}//namespace pmw
