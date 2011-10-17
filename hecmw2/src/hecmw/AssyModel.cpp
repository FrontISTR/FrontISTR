/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/AssyModel.cpp
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
#include "AssyModel.h"
#include "AssyVector.h"
using namespace pmw;
CAssyModel::CAssyModel(void)
{
}
CAssyModel::~CAssyModel(void)
{
    for_each(mvMesh.begin(), mvMesh.end(), DeleteObject());
    for_each(mvContactMesh.begin(), mvContactMesh.end(), DeleteObject());
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
void CAssyModel::GeneLinearAlgebra(const vuint& vDOF, CAssyModel *pCoarseAssyModel)
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
