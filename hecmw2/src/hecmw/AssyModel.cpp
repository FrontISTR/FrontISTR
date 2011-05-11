//
//  AssyModel.cpp
//
//			2009.04.20
//			2008.11.05
//			k.Takeda
#include "AssyModel.h"
#include "AssyVector.h"
using namespace pmw;

CAssyModel::CAssyModel(void)
{
//    mvContactMesh.resize(0);
//    mvMesh.resize(0);
}

CAssyModel::~CAssyModel(void)
{
    // Mesh, ContactMesh
    for_each(mvMesh.begin(), mvMesh.end(), DeleteObject());
    for_each(mvContactMesh.begin(), mvContactMesh.end(), DeleteObject());

    ////debug
    //cout << "~CAssyModel" << endl;
}

// Bucket経由によって,MeshID番号からMeshを取得
//
CMesh* CAssyModel::getMesh_ID(const uiint& id)
{
    uiint index= moBucketMesh.getIndexMesh(id);

    return mvMesh[index];
}

// ContactMeshの追加
// --
void CAssyModel::addContactMesh(CContactMesh* pContactMesh, const uiint& id)
{
    mvContactMesh.push_back(pContactMesh);

    uiint index= mvContactMesh.size()-1;
    mmContactID2Index[id]= index;
}


// Linear Algebra Equations
//
// Ax=b
//
void CAssyModel::GeneLinearAlgebra(const vuint& vDOF, CAssyModel *pCoarseAssyModel)//, vvuint& vPartsID)
{
    mNumOfEquation = vDOF.size();
    
    mvAssyMatrix = new CAssyMatrix* [mNumOfEquation];
    mvRHSAssyVector = new CAssyVector* [mNumOfEquation];
    mvSolAssyVector = new CAssyVector* [mNumOfEquation];
    
    uiint i;
    for(i=0; i < mNumOfEquation; i++){
        uiint nDOF = vDOF[i];
        mvAssyMatrix[i] = new CAssyMatrix(this, nDOF);
        
        if(mMGLevel > 0){
            mvAssyMatrix[i]->setCoarseMatrix(pCoarseAssyModel->getAssyMatrix(i));
        }else{
            mvAssyMatrix[i]->setCoarseMatrix(NULL);
        }
        
        mvRHSAssyVector[i] = new CAssyVector(this, nDOF);
        mvSolAssyVector[i] = new CAssyVector(this, nDOF);
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










