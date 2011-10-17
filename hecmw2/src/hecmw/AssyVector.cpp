/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/AssyVector.cpp
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
CAssyVector::CAssyVector(CAssyModel *pAssyModel, const uiint& nDOF)
: mpAssyModel( pAssyModel )
{
#ifdef ADVANCESOFT_DEBUG
    printf(" enter CAssyVector::CAssyVector %d \n", pAssyModel->getNumOfMesh());
#endif
        mnDOF = nDOF;
	for (int i = 0; i < pAssyModel->getNumOfMesh(); i++) {
            CVector *pVector = new CVector( pAssyModel->getMesh(i), nDOF );
            mvVector.push_back(pVector);
	}
#ifdef ADVANCESOFT_DEBUG
    printf(" exit CAssyVector::CAssyVector \n");
#endif
}
CAssyVector::CAssyVector(const CAssyVector *pAssyVector)
{
	mpAssyModel = pAssyVector->mpAssyModel;
	for (int i = 0; i < pAssyVector->getNumOfVector(); i++) {
		CVector *pVector = new CVector( pAssyVector->getVector(i) );
		mvVector.push_back(pVector);
	}
        mnDOF = pAssyVector->mnDOF;
}
CAssyVector::~CAssyVector()
{
}
uiint& CAssyVector::getDOF()
{
    return mnDOF;
}
uiint CAssyVector::size() const
{
	uiint sum = 0;
	for (CVVecConstIter icv = mvVector.begin(); icv != mvVector.end(); icv++) {
		sum += (*icv)->size();
	}
	return sum;
}
const CVector::ElemType &CAssyVector::operator[](uiint idx) const
{
	for (CVVecConstIter icv = mvVector.begin(); icv != mvVector.end(); icv++) {
		uiint len = (*icv)->size();
		if (idx < len) return (**icv)[idx];
		idx -= len;
	}
}
CVector::ElemType &CAssyVector::operator[](uiint idx)
{
	for (CVVecIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
		uiint len = (*iv)->size();
		if (idx < len) return (**iv)[idx];
		idx -= len;
	}
}
const double &CAssyVector::operator()(uiint meshID, uiint nodeID, uiint dof) const
{
	return (*mvVector[meshID])[nodeID](dof);
}
double &CAssyVector::operator()(uiint meshID, uiint nodeID, uiint dof)
{
	return (*mvVector[meshID])[nodeID](dof);
}
void CAssyVector::Vector_Clear(const uiint& iMesh)
{
    mvVector[iMesh]->Vector_Clear();
}
void CAssyVector::setZero()
{
	for (CVVecIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
		(*iv)->setZero();
	}
}
void CAssyVector::setValue(uiint imesh, uiint inode, uiint idof, double value)
{
	mvVector[imesh]->setValue(inode, idof, value);
}
void CAssyVector::addValue(const uiint& imesh, const uiint& inode, const uiint& idof, const double& value)
{
    mvVector[imesh]->addValue(inode, idof, value);
}
double& CAssyVector::getValue(uiint imesh, uiint inode, uiint idof)
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
  uiint imesh;
  uiint myrank = pMPI->getRank();
  imesh = 0;
  for (CVVecConstIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
    CMesh *pMesh= mpAssyModel->getMesh(imesh);
    CIndexBucket *pBucket= pMesh->getBucket();
    uiint numOfCommMesh2= pMesh->getCommMesh2Size();
    for(uiint icomm=0; icomm< numOfCommMesh2; icomm++){
      CCommMesh2 *pCommMesh2 = pMesh->getCommMesh2IX(icomm);
      iint source      = pCommMesh2->getRank();
      iint destination = pCommMesh2->getTrasmitRank();
      uiint ic = 0;
      if( source < destination ) {
	uiint numOfCommNode= pCommMesh2->getCommNodeSize();
        if( numOfCommNode*mnDOF > N_COMM_BUFF )
	  printf("Fatal Communication Error updateCommBoundary (%ld > %d)"
		 ,numOfCommNode*mnDOF, N_COMM_BUFF);
	for(uiint icnode=0; icnode < numOfCommNode; icnode++){
	  CCommNode *pCommNode= pCommMesh2->getCommNodeIX(icnode);
          uiint nNodeID= pCommNode->getNode()->getID();
	  uiint inode = pBucket->getIndexNode(nNodeID);
          for(uiint idof=0; idof < mnDOF; idof++)
              BUFF[ic++] = (*iv)->getValue(inode, idof);
	}
	pMPI->Send(BUFF, numOfCommNode*mnDOF, MPI_DOUBLE, destination, 100+myrank, MPI_COMM_WORLD);
      }
    }
    imesh++;
  }
  imesh = 0;
  for (CVVecConstIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
    CMesh *pMesh= mpAssyModel->getMesh(imesh);
    CIndexBucket *pBucket= pMesh->getBucket();
    uiint numOfCommMesh2= pMesh->getCommMesh2Size();
    for(uiint icomm=0; icomm< numOfCommMesh2; icomm++){
      CCommMesh2 *pCommMesh2 = pMesh->getCommMesh2IX(icomm);
      iint source      = pCommMesh2->getRank();
      iint destination = pCommMesh2->getTrasmitRank();
      uiint ic = 0;
      if( source > destination ) {
	uiint numOfCommNode= pCommMesh2->getCommNodeSize();
	if( numOfCommNode*mnDOF > N_COMM_BUFF )
	  printf("Fatal Communication Error updateCommBoundary (%ld > %d)"
		 ,numOfCommNode*mnDOF, N_COMM_BUFF);
	pMPI->Recv(BUFF, numOfCommNode*mnDOF, MPI_DOUBLE, destination, 100+destination, MPI_COMM_WORLD, &stat);
	for(uiint icnode=0; icnode< numOfCommNode; icnode++){
	  CCommNode *pCommNode= pCommMesh2->getCommNodeIX(icnode);
          uiint nNodeID= pCommNode->getNode()->getID();
	  uiint inode = pBucket->getIndexNode(nNodeID);
          for(uiint idof=0; idof < mnDOF; idof++)
              (*iv)->setValue(inode, idof, BUFF[ic++]);
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
  uiint imesh;
  imesh = 0;
  for (CVVecConstIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
    CMesh *pMesh= mpAssyModel->getMesh(imesh);
    CIndexBucket *pBucket= pMesh->getBucket();
    uiint numOfCommMesh2= pMesh->getCommMesh2Size();
    for(uiint icomm=0; icomm< numOfCommMesh2; icomm++){
      CCommMesh2 *pCommMesh2 = pMesh->getCommMesh2IX(icomm);
      iint source      = pCommMesh2->getRank();
      iint destination = pCommMesh2->getTrasmitRank();
      uiint ic = 0;
      if( source > destination ) {
	uiint numOfCommNode= pCommMesh2->getCommNodeSize();
	if( numOfCommNode*mnDOF > N_COMM_BUFF )
	  printf("Fatal Communication Error updateCommBoundary (%ld > %d)"
		 ,numOfCommNode*mnDOF, N_COMM_BUFF);
	for(uiint icnode=0; icnode< numOfCommNode; icnode++){
	  CCommNode *pCommNode= pCommMesh2->getCommNodeIX(icnode);
          uiint nNodeID = pCommNode->getNode()->getID();
	  uiint inode = pBucket->getIndexNode(nNodeID);
          for(uiint idof=0; idof < mnDOF; idof++){
              BUFF[ic++] = (*iv)->getValue(inode, idof);
              (*iv)->setValue(inode, idof, 0.0);
          };
	};
	pMPI->Send(BUFF, numOfCommNode*mnDOF, MPI_DOUBLE,destination,101,MPI_COMM_WORLD);
      }
    }
    imesh++;
  }
  imesh = 0;
  for (CVVecConstIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
    CMesh *pMesh= mpAssyModel->getMesh(imesh);
    CIndexBucket *pBucket= pMesh->getBucket();
    uiint numOfCommMesh2= pMesh->getCommMesh2Size();
    for(uiint icomm=0; icomm< numOfCommMesh2; icomm++){
      CCommMesh2 *pCommMesh2 = pMesh->getCommMesh2IX(icomm);
      iint source      = pCommMesh2->getRank();
      iint destination = pCommMesh2->getTrasmitRank();
      uiint ic = 0;
      if( source < destination ) {
	uiint numOfCommNode= pCommMesh2->getCommNodeSize();
	if( numOfCommNode*mnDOF > N_COMM_BUFF )
	  printf("Fatal Communication Error updateCommBoundary (%ld > %d)"
		 ,numOfCommNode*mnDOF, N_COMM_BUFF);
	pMPI->Recv(BUFF, numOfCommNode*mnDOF, MPI_DOUBLE,destination,101,MPI_COMM_WORLD,&stat);
	for(uiint icnode=0; icnode< numOfCommNode; icnode++){
	  CCommNode *pCommNode= pCommMesh2->getCommNodeIX(icnode);
          uiint nNodeID= pCommNode->getNode()->getID();
	  uiint inode = pBucket->getIndexNode(nNodeID);
          for(uiint idof=0; idof < mnDOF; idof++)
              (*iv)->addValue(inode, idof, BUFF[ic++]);
	}
      }
    }
    imesh++;
  }
}
uiint CAssyVector::restrictTo(CAssyVector *pVc) const
{
	CVVecIter ivc = pVc->mvVector.begin();
	for (CVVecConstIter icv = mvVector.begin(); icv != mvVector.end(); icv++) {
		(*icv)->restrictTo((*ivc));
		ivc++;
	}
	return 1;
}
uiint CAssyVector::prolongateFrom(const CAssyVector *pcVc)
{
	CVVecConstIter icvc = pcVc->mvVector.begin();
	for (CVVecIter iv = mvVector.begin(); iv != mvVector.end(); iv++) {
		(*iv)->prolongateFrom((*icvc));
		icvc++;
	}
	return 1;
}
uiint CAssyVector::getNumOfVector() const
{
	return mvVector.size();
}
const CVector *CAssyVector::getVector(uiint index) const
{
	return mvVector[index];
}
CVector *CAssyVector::getVector(uiint index)
{
	return mvVector[index];
}
void CAssyVector::dump()
{
    uiint nNumOfParts = mvVector.size();
    uiint ipart;
    for(ipart=0; ipart < nNumOfParts; ipart++){
        cout << " ---- iParts : " << ipart << " Vector ---- " << endl;
        mvVector[ipart]->dump();
    };
    cout << endl;
}
}
