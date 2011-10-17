/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/Vector.cpp
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
#include "Vector.h"
#include "Mesh.h"
#include "Node.h"
namespace pmw
{
CVector::CVector(CMesh *pMesh, const uiint& nDOF)
{
#ifdef ADVANCESOFT_DEBUG
    printf(" enter CVector::CVector \n");
#endif
  mpMesh = pMesh;
  mnDOF = nDOF;
  mnNode = pMesh->getNumOfNode();
  mvVector.resize(mnNode);
  for (uiint i = 0; i < mnNode; i++) {
    mvVector[i].resize(mnDOF);
  }
#ifdef ADVANCESOFT_DEBUG
    printf(" exit CVector::CVector \n");
#endif
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
uiint CVector::size() const
{
	return mnNode;
}
bool CVector::isScopeNode(const uiint& idx) const
{
    uiint max = mvVector.size()-1;
    if(idx > max){
        return false;
    }else{
        return true;
    }
}
const CVector::ElemType& CVector::operator[](uiint idx) const
{
    return mvVector[idx];
}
CVector::ElemType& CVector::operator[](uiint idx)
{
    return mvVector[idx];
}
void CVector::Vector_Clear()
{
    uiint i,idof;
    for(i=0; i < mnNode; i++){
        for(idof=0; idof < mnDOF; idof++){
            mvVector[i][idof] = 0.0;
        };
    };
}
void CVector::setZero()
{
	for (uiint i = 0; i < mnNode; i++) {
		mvVector[i].clear();
	}
}
void CVector::setValue(uiint inode, uiint idof, double value)
{
	mvVector[inode][idof] = value;
}
void CVector::addValue(uiint inode, uiint idof, double value)
{
	mvVector[inode][idof] += value;
}
double& CVector::getValue(uiint inode, uiint idof)
{
	return( mvVector[inode][idof] );
}
void CVector::sumSV(double alpha, const CVector *pX, CVector *pY) const
{
	for (uiint i = 0; i < mnNode; i++) {
		pY->mvVector[i] = mvVector[i] + alpha * pX->mvVector[i];
	}
}
void CVector::addSV(double alpha, const CVector *pX)
{
	for (uiint i = 0; i < mnNode; i++) {
		mvVector[i] += alpha * pX->mvVector[i];
	}
}
void CVector::add(const CVector *pX)
{
	for (uiint i = 0; i < mnNode; i++) {
		mvVector[i] += pX->mvVector[i];
	}
}
void CVector::subst(const CVector *pX)
{
	for (uiint i = 0; i < mnNode; i++) {
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
	for (uiint i = 0; i < mnNode; i++) {
		for (uiint j = 0; j < mnDOF; j++) {
			sum += mvVector[i](j) * pX->mvVector[i](j);
		}
	}
	return sum;
}
void CVector::updateCommBoundary()
{
}
uiint CVector::restrictTo(CVector *pV) const
{
	for( uiint i=0; i< mpMesh->getNumOfNode(); i++) {
		CNode* node = mpMesh->getNodeIX(i);
		uiint numP = node->getNumOfParentNode();
		if( numP == 0 ) {
			(*pV)[i] = mvVector[i];
		}
	}
	return 0;
}
uiint CVector::prolongateFrom(const CVector *pV)
{
    vector<uiint> vQuadN;
    for( uiint i=0; i< mpMesh->getNumOfNode(); i++) {
        CNode* node = mpMesh->getNodeIX(i);
        uiint numP = node->getNumOfParentNode();
        if( numP == 0 ) {
            mvVector[i] = (*pV)[i];
        } else {
            for(uiint idof=0; idof < mnDOF; idof++) mvVector[i](idof) = 0.0;
            bool bQuad(false);
            if(numP == 2){
                for(uiint j=0; j < numP; j++){
                    uiint k = node->getParentNode(j)->getID();
                    if( !pV->isScopeNode(k) ) bQuad=true;
                };
            }
            if(bQuad) vQuadN.push_back(i);
            if( !bQuad ){
                for(uiint j=0; j < numP; j++) {
                    uiint k = node->getParentNode(j)->getID();
                    mvVector[i] += (*pV)[k];
                };
                mvVector[i] = mvVector[i] / numP;
            }
        }
    };
    return 0;
}
void CVector::dump()
{
    uiint i,j;
    for(i = 0; i < mnNode; i++){
        for(j = 0; j < mnDOF; j++){
                cout << mvVector[i](j) << " ";
        };
    };
    cout << endl;
}
}
