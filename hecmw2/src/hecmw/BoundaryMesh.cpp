/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/BoundaryMesh.cpp
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
#include <vector>
#include "BoundaryMesh.h"
using namespace pmw;
CBoundaryMesh::CBoundaryMesh()
{
    mnEdgeNodeCount = 0;
}
CBoundaryMesh::~CBoundaryMesh()
{
    if(mMaxMGLevel==mMGLevel) for_each(mvBNode.begin(), mvBNode.end(), DeleteObject());
}
void CBoundaryMesh::addDOF(const uiint& dof)
{
    mvDOF.push_back(dof);
    mmDOF2Index[dof] = mvDOF.size() - 1;
}
void CBoundaryMesh::setDOF(const uiint& index, const uiint& dof)
{
    mvDOF[index] = dof;
    mmDOF2Index[dof] = index;
}
void CBoundaryMesh::resizeDOF(const uiint& res_size)
{
    mvDOF.resize(res_size);
}
uiint& CBoundaryMesh::getDOF(const uiint& index)
{
    return mvDOF[index];
}
uiint& CBoundaryMesh::getDOF_Index(const uiint& dof)
{
    return mmDOF2Index[dof];
}
uiint CBoundaryMesh::getNumOfDOF()
{
    return mvDOF.size();
}
void CBoundaryMesh::resizeCGrid_BNodeValue(const uiint& maxLevel)
{
    uiint i, nNumOfBNode= mvBNode.size();
    for(i=0; i < nNumOfBNode; i++){
        mvBNode[i]->resizeValue(maxLevel+1);
    };
}
void CBoundaryMesh::setBNode(const uiint& index, CBoundaryNode* pBNode)
{
    uiint id;
    id= pBNode->getID();
    mmBNodeID2Index[id]= index;
    mvBNode[index]= pBNode;
}
void CBoundaryMesh::addBNode(CBoundaryNode* pBNode)
{
    mvBNode.push_back(pBNode);
    uiint id;
    id= pBNode->getID();
    mmBNodeID2Index[id]= mvBNode.size()-1;
}
void CBoundaryMesh::distValueBNode()
{
    switch(mnBndType){
        case(BoundaryType::Neumann):
            distNeumannValue();
            break;
        case(BoundaryType::Dirichlet):
            distDirichletValue();
            break;
    }
}
