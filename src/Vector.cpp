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

CVector::CVector(CMesh *pMesh, const uint& nDOF)
{
#ifdef ADVANCESOFT_DEBUG
    printf(" enter CVector::CVector \n");
#endif

  mpMesh = pMesh;
  mnDOF = nDOF;
  mnNode = pMesh->getNumOfNode();
  // mnNodeInternal = pMesh->???;  // TODO check if this is necessary
  
  mvVector.resize(mnNode);
  for (int i = 0; i < mnNode; i++) {
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

//
// Vector範囲内に存在するか. 
//  *  prolongateFrom で利用 : 2次ノードの判定(2次ノードの親ノードはコースグリッドに存在しない)
//
bool CVector::isScopeNode(const uint& idx) const
{
    // 配列外を弾く
    //  *2次要素での prolongateFrom 2次ノードはコースグリッドに存在しない場合
    //
    uint max = mvVector.size()-1;
    if(idx > max){
        return false;
    }else{
        return true;
    }
}
const CVector::ElemType& CVector::operator[](size_t idx) const
{
    return mvVector[idx];
}

CVector::ElemType& CVector::operator[](size_t idx)
{
    return mvVector[idx];
}

// Matrix 0 clear
//
void CVector::Vector_Clear()
{
    uint i,idof;
    for(i=0; i < mnNode; i++){
        for(idof=0; idof < mnDOF; idof++){
            mvVector[i][idof] = 0.0;
        };
    };
}

void CVector::setZero()// use GMRES
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
    
    vector<uint> vQuadN;//2次ノードのインデックス番号 配列

    for( uint i=0; i< mpMesh->getNumOfNode(); i++) {
        CNode* node = mpMesh->getNodeIX(i);
        uint numP = node->getNumOfParentNode();

        if( numP == 0 ) {
            mvVector[i] = (*pV)[i];
        } else {
            //mvVector[i](0) = 0.0;mvVector[i](1) = 0.0;mvVector[i](2) = 0.0;
            for(uint idof=0; idof < mnDOF; idof++) mvVector[i](idof) = 0.0;// 2010.11.30

            //2次ノード判定
            bool bQuad(false);
            if(numP == 2){
                for(uint j=0; j < numP; j++){
                    uint k = node->getParentNode(j)->getID();
                    if( !pV->isScopeNode(k) ) bQuad=true;
                };
            }
            if(bQuad) vQuadN.push_back(i);//2次ノードのノード・インデックス番号を収集

            // 頂点のノードのプロロンゲート
            if( !bQuad ){
                for(uint j=0; j < numP; j++) {
                    uint k = node->getParentNode(j)->getID();

                    mvVector[i] += (*pV)[k];
                };
                mvVector[i] = mvVector[i] / numP;
            }

        }// if(numP==0):コースグリッド
    };

////    // 2次ノードの処理はあってもなくても変わらない:所詮推測値を入れているに過ぎないから
////    //
////    cout << "vQuadN.size == " << vQuadN.size() << endl;
////
////    // 2次ノードの値を同一レベルの頂点からセット
////    for(uint i=0; i < vQuadN.size(); i++){
////        uint idex = vQuadN[i];
////        CNode* pNode = mpMesh->getNodeIX(idex);
////        uint n1 = pNode->getParentNode(0)->getID();
////        uint n2 = pNode->getParentNode(1)->getID();
////
////        mvVector[idex] = (mvVector[n1] + mvVector[n2]) * 0.5;
////    }

    return 0;//2010.05.14
}

}//namespace pmw
