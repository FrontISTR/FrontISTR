/*
 * Vector.cpp
 *
 *  Created on: Jul 23, 2009
 *      Author: goto
 */

#include "Vector.h"
#include "Mesh.h"
#include "Node.h"

namespace pmw
{

CVector::CVector(/* const */ CMesh *pMesh)
{
  //   	printf(" enter CVector::CVector \n");

  mpMesh = pMesh;
  CNode *pNode = pMesh->getNode(0);
  mnDOF = pNode->numOfTotalParam();// '10.10.06 修正 k.Takeda
  mnNode = pMesh->getNumOfNode();
  // mnNodeInternal = pMesh->???;  // TODO check if this is necessary
  
  mvVector.resize(mnNode);
  for (int i = 0; i < mnNode; i++) {
    mvVector[i].resize(mnDOF);
  }
  
  //   	printf(" exit CVector::CVector \n");
}

CVector::CVector(const CVector *pVector)
{
	mnDOF = pVector->mnDOF;
	mnNode = pVector->mnNode;
	mnNodeInternal = pVector->mnNodeInternal;

	mvVector.resize(mnNode);
	for (int i = 0; i < mnNode; i++) {
		mvVector[i].resize(mnDOF);
	}
}

CVector::~CVector()
{
	// TODO Auto-generated destructor stub
}

size_t CVector::size() const
{
	return mnNode;
}

//int CVector::lenInternal() const
//{
//	//
//	return mnNodeInternal * mnDOF;
//}

const CVector::ElemType &CVector::operator[](size_t idx) const
{
	return mvVector[idx];
}

CVector::ElemType &CVector::operator[](size_t idx)
{
	return mvVector[idx];
}

void CVector::setZero()
{
	for (int i = 0; i < mnNode; i++) {
		mvVector[i].clear();
	}
}

void CVector::setValue(int inode, int idof, double value)
{
	mvVector[inode][idof] = value;
}

void CVector::addValue(int inode, int idof, double value)
{
	mvVector[inode][idof] += value;
}

double CVector::getValue(int inode, int idof)
{
	return( mvVector[inode][idof] );
}

void CVector::sumSV(double alpha, const CVector *pX, CVector *pY) const
{
	for (int i = 0; i < mnNode; i++) {
		pY->mvVector[i] = mvVector[i] + alpha * pX->mvVector[i];
	}
}

void CVector::addSV(double alpha, const CVector *pX)
{
	for (int i = 0; i < mnNode; i++) {
		mvVector[i] += alpha * pX->mvVector[i];
	}
}

void CVector::add(const CVector *pX)
{
	for (int i = 0; i < mnNode; i++) {
		mvVector[i] += pX->mvVector[i];
	}
}

void CVector::subst(const CVector *pX)
{
	for (int i = 0; i < mnNode; i++) {
		mvVector[i] = pX->mvVector[i];
	}
}

double CVector::norm2() const
{
	return innerProd(this);
}

double CVector::innerProd(const CVector *pX) const
{
	double sum = 0.0;
	for (int i = 0; i < mnNode; i++) {
		for (int j = 0; j < mnDOF; j++) {
			sum += mvVector[i](j) * pX->mvVector[i](j);
		}
	}
	return sum;
}

void CVector::updateCommBoundary()
{
	// TODO implement CVector::updateCommBoundary
}

int CVector::restrictTo(CVector *pV) const
{
	// TODO implement CVector::restrictTo
	
//	const CMesh *pMesh = getMesh();
	for( uint i=0; i< mpMesh->getNumOfNode(); i++) {
		CNode* node = mpMesh->getNodeIX(i);
		uint numP = node->getNumOfParentNode();
		if( numP == 0 ) {
			(*pV)[i] = mvVector[i];
		}
	}

	return 0;//2010.05.14
}

int CVector::prolongateFrom(const CVector *pV)
{
	// TODO implement CVector::prolongateFrom

	for( uint i=0; i< mpMesh->getNumOfNode(); i++) {
            CNode* node = mpMesh->getNodeIX(i);
            uint numP = node->getNumOfParentNode();
            if( numP == 0 ) {
                mvVector[i] = (*pV)[i];
            } else {
                mvVector[i](0) = 0.0;mvVector[i](1) = 0.0;mvVector[i](2) = 0.0;
                for(uint j=0; j<numP; j++) {
                    uint k = node->getParentNode(j)->getID();
                    mvVector[i] += (*pV)[k];
                }
                mvVector[i] = mvVector[i] / numP;
            }
	}

	return 0;//2010.05.14
}

}//namespace pmw
