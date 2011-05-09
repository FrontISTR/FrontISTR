/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Vector.cxx
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
  mpMesh = pMesh;
  CNode *pNode = pMesh->getNode(0);
  mnDOF = pNode->numOfVectorParam();
  mnNode = pMesh->getNumOfNode();
  mvVector.resize(mnNode);
  for (int i = 0; i < mnNode; i++) {
    mvVector[i].resize(mnDOF);
  }
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
}
size_t CVector::size() const
{
	return mnNode;
}
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
}
int CVector::restrictTo(CVector *pV) const
{
	for( uint i=0; i< mpMesh->getNumOfNode(); i++) {
		CNode* node = mpMesh->getNodeIX(i);
		uint numP = node->getNumOfParentNode();
		if( numP == 0 ) {
			(*pV)[i] = mvVector[i];
		}
	}
	return 0;
}
int CVector::prolongateFrom(const CVector *pV)
{
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
	return 0;
}
}
