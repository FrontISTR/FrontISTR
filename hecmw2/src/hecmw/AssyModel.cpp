/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/AssyModel.cpp
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "AssyModel.h"
#include "AssyVector.h"
using namespace pmw;
CAssyModel::CAssyModel(void)
{
////    //--Testデータ begin--
////    mnNumOfGlobalComm= 4;
////
////    // rank pair
////    mvRankPair.resize(mnNumOfGlobalComm);
////    for(uiint i=0; i < mnNumOfGlobalComm; i++){
////        pair<uiint,uiint> rank_pair;
////        switch(i){
////            case(0): rank_pair.first= 1; rank_pair.second=2; break;
////            case(1): rank_pair.first= 0; rank_pair.second=1; break;
////            case(2): rank_pair.first= 1; rank_pair.second=2; break;
////            case(3): rank_pair.first= 0; rank_pair.second=1; break;
////        }
////        mvRankPair[i]= rank_pair;
////    };
////
////    // MeshID for each CommID
////    mvMeshID_CommID.resize(mnNumOfGlobalComm);
////    mvMeshID_CommID[0]= 0;//CommID:0 は Mesh:0についてる.
////    mvMeshID_CommID[1]= 0;//CommID:1 は Mesh:0についてる.
////    mvMeshID_CommID[2]= 1;//CommID:2 は Mesh:1についてる.
////    mvMeshID_CommID[3]= 1;//CommID:3 は Mesh:1についてる.
////    //--Testデータ end--
}
CAssyModel::~CAssyModel(void)
{
    for_each(mvMesh.begin(), mvMesh.end(), DeleteObject());
    for_each(mvContactMesh.begin(), mvContactMesh.end(), DeleteObject());
}
void CAssyModel::setNumGlobalCommMesh2(const uiint& nNumOfGlobalComm)
{
    mnNumOfGlobalComm= nNumOfGlobalComm;

    mvRankPair.resize(mnNumOfGlobalComm);
    mvMeshID_CommID.resize(mnNumOfGlobalComm);
}
void CAssyModel::setGlobalPairRank(const uiint& nCommID, const uiint& nFirstRank, const uiint& nSecondRank)
{
    mvRankPair[nCommID].first = nFirstRank;
    mvRankPair[nCommID].second= nSecondRank;
}
void CAssyModel::setMeshID_with_CommID(const uiint& nCommID, const uiint& nMeshID)
{
    mvMeshID_CommID[nCommID]= nMeshID;
}


CMesh* CAssyModel::getMesh_ID(const uiint& id)
{
    uiint index= moBucketMesh.getIndexMesh(id);
    return mvMesh[index];
}
void CAssyModel::addContactMesh(CContactMesh* pContactMesh, const uiint& id)
{
    mvContactMesh.push_back(pContactMesh);
    uiint index= mvContactMesh.size()-1;
    mmContactID2Index[id]= index;
}
//--
// 方程式Ax=b生成 : vvDOF[方程式番号][全域Mesh番号]= DOF : 全域Mesh番号=MeshID
//--
void CAssyModel::GeneLinearAlgebra(vvuint& vvDOF, CAssyModel *pCoarseAssyModel)
{
    //
    // vvDOF: 方程式数 - nGlobalNumOfMesh の2重配列
    //
    mNumOfEquation = vvDOF.size();
    mvAssyMatrix = new CAssyMatrix* [mNumOfEquation];
    mvRHSAssyVector = new CAssyVector* [mNumOfEquation];
    mvSolAssyVector = new CAssyVector* [mNumOfEquation];

    for(uiint i=0; i < mNumOfEquation; i++) {
        vuint vDOF = vvDOF[i];//---------- グローバルMesh数ぶんのDOF
        mvAssyMatrix[i] = new CAssyMatrix(this, vDOF, i);
        mvRHSAssyVector[i] = new CAssyVector(this, vDOF);
        mvSolAssyVector[i] = new CAssyVector(this, vDOF);
        //--
        // コースグリッド設定(Matrix,Vector) : VectorはMGCycle用途
        //--
        if(mMGLevel > 0) {
            mvAssyMatrix[i]->setCoarseMatrix(pCoarseAssyModel->getAssyMatrix(i));
            mvRHSAssyVector[i]->setCoarseVector(pCoarseAssyModel->getRHSAssyVector(i));     //-- '11.11.29
            mvSolAssyVector[i]->setCoarseVector(pCoarseAssyModel->getSolutionAssyVector(i));//-- '11.11.29
        } else {
            mvAssyMatrix[i]->setCoarseMatrix(NULL);
            mvRHSAssyVector[i]->setCoarseVector(NULL);//-- '11.11.29
            mvSolAssyVector[i]->setCoarseVector(NULL);//-- '11.11.29
        }
    };
}


CAssyMatrix* CAssyModel::getAssyMatrix(const uiint& index)
{
    return mvAssyMatrix[index];
}
CAssyVector* CAssyModel::getRHSAssyVector(const uiint& index)
{
    return mvRHSAssyVector[index];
}
CAssyVector* CAssyModel::getSolutionAssyVector(const uiint& index)
{
    return mvSolAssyVector[index];
}
//--
// IDのMeshは存在するか？
//--
bool CAssyModel::isSelfMesh(const uiint& id)
{
    bool bSelfMesh(false);
    uiint nNumOfMesh= mvMesh.size();

    iint nLow= 0;
    iint nHigh= nNumOfMesh-1;

    while(nLow <= nHigh) {
        uiint i= ( nLow + nHigh )/2;
        if( id == mvMesh[i]->getMeshID() ) {
            bSelfMesh=true;
            break;
        } else if( id < mvMesh[i]->getMeshID() ) {
            nHigh= i-1;
        } else {
            nLow = i+1;
        }
    };

    return bSelfMesh;
}





