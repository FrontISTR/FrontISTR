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

    //debug
    cout << "~CAssyModel" << endl;
}

// Bucket経由によって,MeshID番号からMeshを取得
//
CMesh* CAssyModel::getMesh_ID(const uint& id)
{
    uint index= moBucketMesh.getIndexMesh(id);

    return mvMesh[index];
}

// ContactMeshの追加
// --
void CAssyModel::addContactMesh(CContactMesh* pContactMesh, const uint& id)
{
    mvContactMesh.push_back(pContactMesh);

    uint index= mvContactMesh.size()-1;
    mmContactID2Index[id]= index;
}


// Linear Algebra Equations
//
// Ax=b
//
void CAssyModel::GeneLinearAlgebra(const vuint& vDOF, CAssyModel *pCoarseAssyModel)//, vvuint& vPartsID)
{
    uint nNumOfEquation = vDOF.size();
    
    mvAssyMatrix = new CAssyMatrix* [nNumOfEquation];
    mvRHSAssyVector = new CAssyVector* [nNumOfEquation];
    mvSolAssyVector = new CAssyVector* [nNumOfEquation];
    
    uint i;
    for(i=0; i < nNumOfEquation; i++){
        uint nDOF = vDOF[i];
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
CAssyMatrix* CAssyModel::getAssyMatrix(const uint& index)
{
    return mvAssyMatrix[index];
}
CAssyVector* CAssyModel::getRHSAssyVector(const uint& index)
{
    return mvRHSAssyVector[index];
}
CAssyVector* CAssyModel::getSolutionAssyVector(const uint& index)
{
    return mvSolAssyVector[index];
}










