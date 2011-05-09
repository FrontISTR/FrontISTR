/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   AssyVector.cxx
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
 * AssyVector.cpp
 *
 *  Created on: Oct 16, 2009
 *      Author: goto
 */
#include "AssyVector.h"
#include "AssyModel.h"
#include "CommNode.h"
#include "CommMesh2.h"
#include "HEC_MPI.h"
namespace pmw
{
typedef std::vector<CVector*> CVVec;
typedef CVVec::iterator CVVecIter;
typedef CVVec::const_iterator CVVecConstIter;
#define N_COMM_BUFF 1024*32
double BUFF[N_COMM_BUFF];
CAssyVector::CAssyVector(CAssyModel *pAssyModel)
: mpAssyModel( pAssyModel )
{
#ifdef ADVANCESOFT_DEBUG
#endif
	for (int i = 0; i < pAssyModel->getNumOfMesh(); i++) {
		CVector *pVector = new CVector( pAssyModel->getMesh(i) );
		mvVector.push_back(pVector);
	}
#ifdef ADVANCESOFT_DEBUG
#endif
}
CAssyVector::CAssyVector(const CAssyVector *pAssyVector)
{
	mpAssyModel = pAssyVector->mpAssyModel;
	for (int i = 0; i < pAssyVector->getNumOfVector(); i++) {
		CVector *pVector = new CVector( pAssyVector->getVector(i) );
		mvVector.push_back(pVector);
	}
}
CAssyVector::~CAssyVector()
{
}
size_t CAssyVector::size() const
{
	int sum = 0;
	for (CVVecConstIter icv = mvVector.begin(); icv != mvVector.end(); icv++) {
		sum += (*icv)->size();
	}
	return sum;
}
const CVector::ElemType &CAssyVector::operator[](size_t idx) const
{
	for (CVVecConstIter icv = mvVector.begin(); icv != mvVector.end(); icv++) {
		int len = (*icv)->size();
		if (idx < len) return (**icv)[idx];
		idx -= len;
	}
}
CVector::ElemType &CAssyVector::operator[](size_t idx)
{
	for (CVVecIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
		int len = (*iv)->size();
		if (idx < len) return (**iv)[idx];
		idx -= len;
	}
}
const double &CAssyVector::operator()(size_t meshID, size_t nodeID, size_t dof) const
{
	return (*mvVector[meshID])[nodeID](dof);
}
double &CAssyVector::operator()(size_t meshID, size_t nodeID, size_t dof)
{
	return (*mvVector[meshID])[nodeID](dof);
}
void CAssyVector::setZero()
{
	for (CVVecIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
		(*iv)->setZero();
	}
}
void CAssyVector::setValue(int imesh, int inode, int idof, double value)
{
	mvVector[imesh]->setValue(inode, idof, value);
}
double CAssyVector::getValue(int imesh, int inode, int idof)
{
	return( mvVector[imesh]->getValue(inode, idof) );
}
void CAssyVector::add(const CAssyVector *pV)
{
	CVVecConstIter icv = pV->mvVector.begin();
	for (CVVecIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
		(*iv)->add((*icv));
		icv++;
	}
}
void CAssyVector::subst(const CAssyVector *pV)
{
	CVVecConstIter icv = pV->mvVector.begin();
	for (CVVecIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
		(*iv)->subst((*icv));
		icv++;
	}
}
double CAssyVector::norm2() const
{
	double sum = 0;
	for (CVVecConstIter icv = mvVector.begin(); icv != mvVector.end(); icv++) {
		sum += (*icv)->norm2();
	}
	pmw::CHecMPI *pMPI= pmw::CHecMPI::Instance();
	double sum0 = sum;
	pMPI->Allreduce(&sum0, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return sum;
}
double CAssyVector::innerProd(const CAssyVector *pX) const
{
	double sum = 0;
	CVVecConstIter icX = pX->mvVector.begin();
	for (CVVecConstIter icv = mvVector.begin(); icv != mvVector.end(); icv++) {
		sum += (*icv)->innerProd((*icX));
		icX++; 
	}
	pmw::CHecMPI *pMPI= pmw::CHecMPI::Instance();
	double sum0 = sum;
	pMPI->Allreduce(&sum0, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return sum;
}
void CAssyVector::updateCommBoundary()
{
  pmw::CHecMPI *pMPI= pmw::CHecMPI::Instance();
  MPI_Status stat;
  uint imesh;
  int myrank = pMPI->getRank();
  imesh = 0;
  for (CVVecConstIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
    CMesh *pMesh= mpAssyModel->getMesh(imesh);
    uint numOfCommMesh2= pMesh->getCommMesh2Size();
    for(uint icomm=0; icomm< numOfCommMesh2; icomm++){
      CCommMesh2 *pCommMesh2 = pMesh->getCommMesh2IX(icomm);
      uint source	 = pCommMesh2->getRank();
      uint destination = pCommMesh2->getTrasmitRank();
      uint ic = 0;
      if( source < destination ) {
	int numOfCommNode= pCommMesh2->getCommVertNodeSize();
	if( numOfCommNode*3 > N_COMM_BUFF )
	  printf("Fatal Communication Error updateCommBoundary (%d > %d)"
		 ,numOfCommNode*3, N_COMM_BUFF);
	for(uint icnode=0; icnode< numOfCommNode; icnode++){
	  CCommNode *pCommNode= pCommMesh2->getCommVertNodeIX(icnode);
	  uint inode = pCommNode->getNode()->getID();
	  BUFF[ic++] = (*iv)->getValue(inode, 0);
	  BUFF[ic++] = (*iv)->getValue(inode, 1);
	  BUFF[ic++] = (*iv)->getValue(inode, 2);
	}
	pMPI->Send(BUFF, numOfCommNode*3, MPI_DOUBLE, destination, 100+myrank, MPI_COMM_WORLD);
      }
    }
    imesh++;
  }
  imesh = 0;
  for (CVVecConstIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
    CMesh *pMesh= mpAssyModel->getMesh(imesh);
    uint numOfCommMesh2= pMesh->getCommMesh2Size();
    for(uint icomm=0; icomm< numOfCommMesh2; icomm++){
      CCommMesh2 *pCommMesh2 = pMesh->getCommMesh2IX(icomm);
      uint source	 = pCommMesh2->getRank();
      uint destination = pCommMesh2->getTrasmitRank();
      uint ic = 0;
      if( source > destination ) {
	int numOfCommNode= pCommMesh2->getCommVertNodeSize();
	if( numOfCommNode*3 > N_COMM_BUFF )
	  printf("Fatal Communication Error updateCommBoundary (%d > %d)"
		 ,numOfCommNode*3, N_COMM_BUFF);
	pMPI->Recv(BUFF, numOfCommNode*3, MPI_DOUBLE, destination, 100+destination, MPI_COMM_WORLD, &stat);
	for(uint icnode=0; icnode< numOfCommNode; icnode++){
	  CCommNode *pCommNode= pCommMesh2->getCommVertNodeIX(icnode);
	  uint inode = pCommNode->getNode()->getID();
	  (*iv)->setValue(inode, 0, BUFF[ic++]);
	  (*iv)->setValue(inode, 1, BUFF[ic++]);
	  (*iv)->setValue(inode, 2, BUFF[ic++]);
	}
      }
    }
    imesh++;
  }
}
void CAssyVector::sumupCommBoundary()
{
  pmw::CHecMPI *pMPI= pmw::CHecMPI::Instance();
  MPI_Status stat;
  uint imesh;
  int myrank = pMPI->getRank();
  imesh = 0;
  for (CVVecConstIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
    CMesh *pMesh= mpAssyModel->getMesh(imesh);
    uint numOfCommMesh2= pMesh->getCommMesh2Size();
    for(uint icomm=0; icomm< numOfCommMesh2; icomm++){
      CCommMesh2 *pCommMesh2 = pMesh->getCommMesh2IX(icomm);
      uint source	 = pCommMesh2->getRank();
      uint destination = pCommMesh2->getTrasmitRank();
      uint ic = 0;
      if( source > destination ) {
	int numOfCommNode= pCommMesh2->getCommVertNodeSize();
	if( numOfCommNode*3 > N_COMM_BUFF )
	  printf("Fatal Communication Error updateCommBoundary (%d > %d)"
		 ,numOfCommNode*3, N_COMM_BUFF);
	for(uint icnode=0; icnode< numOfCommNode; icnode++){
	  CCommNode *pCommNode= pCommMesh2->getCommVertNodeIX(icnode);
	  uint inode = pCommNode->getNode()->getID();
	  BUFF[ic++] = (*iv)->getValue(inode, 0);
	  BUFF[ic++] = (*iv)->getValue(inode, 1);
	  BUFF[ic++] = (*iv)->getValue(inode, 2);
	  (*iv)->setValue(inode, 0, 0.0);
	  (*iv)->setValue(inode, 1, 0.0);
	  (*iv)->setValue(inode, 2, 0.0);
	}
	pMPI->Send(BUFF, numOfCommNode*3, MPI_DOUBLE,destination,101,MPI_COMM_WORLD);
      }
    }
    imesh++;
  }
  imesh = 0;
  for (CVVecConstIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
    CMesh *pMesh= mpAssyModel->getMesh(imesh);
    uint numOfCommMesh2= pMesh->getCommMesh2Size();
    for(uint icomm=0; icomm< numOfCommMesh2; icomm++){
      CCommMesh2 *pCommMesh2 = pMesh->getCommMesh2IX(icomm);
      uint source	 = pCommMesh2->getRank();
      uint destination = pCommMesh2->getTrasmitRank();
      uint ic = 0;
      if( source < destination ) {
	int numOfCommNode= pCommMesh2->getCommVertNodeSize();
	if( numOfCommNode*3 > N_COMM_BUFF )
	  printf("Fatal Communication Error updateCommBoundary (%d > %d)"
		 ,numOfCommNode*3, N_COMM_BUFF);
	pMPI->Recv(BUFF, numOfCommNode*3, MPI_DOUBLE,destination,101,MPI_COMM_WORLD,&stat);
	for(uint icnode=0; icnode< numOfCommNode; icnode++){
	  CCommNode *pCommNode= pCommMesh2->getCommVertNodeIX(icnode);
	  uint inode = pCommNode->getNode()->getID();
	  (*iv)->addValue(inode, 0, BUFF[ic++]);
	  (*iv)->addValue(inode, 1, BUFF[ic++]);
	  (*iv)->addValue(inode, 2, BUFF[ic++]);
	}
      }
    }
    imesh++;
  }
}
int CAssyVector::restrictTo(CAssyVector *pVc) const
{
	CVVecIter ivc = pVc->mvVector.begin();
	for (CVVecConstIter icv = mvVector.begin(); icv != mvVector.end(); icv++) {
		(*icv)->restrictTo((*ivc));
		ivc++;
	}
	return 1;
}
int CAssyVector::prolongateFrom(const CAssyVector *pcVc)
{
	CVVecConstIter icvc = pcVc->mvVector.begin();
	for (CVVecIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
		(*iv)->prolongateFrom((*icvc));
		icvc++;
	}
	return 1;
}
size_t CAssyVector::getNumOfVector() const
{
	return mvVector.size();
}
const CVector *CAssyVector::getVector(size_t index) const
{
	return mvVector[index];
}
CVector *CAssyVector::getVector(size_t index)
{
	return mvVector[index];
}
}
