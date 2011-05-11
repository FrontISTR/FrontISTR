//
//  BoundaryMesh.cpp
//
//
//              2010.05.06
//              k.Takeda
#include <vector>

#include "BoundaryMesh.h"
using namespace pmw;


CBoundaryMesh::CBoundaryMesh()
{
    ////mbSecondOrder = false;
    mnEdgeNodeCount = 0;
}
CBoundaryMesh::~CBoundaryMesh()
{
    if(mMaxMGLevel==mMGLevel) for_each(mvBNode.begin(), mvBNode.end(), DeleteObject());
}

//---
//DOFの追加・代入・配列領域
//---
//DOFの追加
void CBoundaryMesh::addDOF(const uiint& dof)
{
    mvDOF.push_back(dof);
    mmDOF2Index[dof] = mvDOF.size() - 1;
}
//DOFの代入
void CBoundaryMesh::setDOF(const uiint& index, const uiint& dof)
{
    mvDOF[index] = dof;
    mmDOF2Index[dof] = index;
}
//DOF配列の領域確保
void CBoundaryMesh::resizeDOF(const uiint& res_size)
{
    mvDOF.resize(res_size);
}
//---
//DOFの提供・DOFインデックス・DOF数
//---
//DOFの提供
uiint& CBoundaryMesh::getDOF(const uiint& index)
{
    return mvDOF[index];
}
//DOFのインデックスの提供
uiint& CBoundaryMesh::getDOF_Index(const uiint& dof)
{
    return mmDOF2Index[dof];
}
//DOF数
uiint CBoundaryMesh::getNumOfDOF()
{
    return mvDOF.size();
}

//CMWのRefineでコースグリッドのBNode境界条件の領域確保に使用 2011.04.22
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


// BNodeへの値の分配 => distNeumann, distDirichlet呼び出し
//
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



////// ディレクレ型境界値のBNodeへの再配分
//////
////void CBoundaryMesh::distDirichletValue()
////{
////    Utility::CLogger *pLogger= Utility::CLogger::Instance();
////
////    // ディレクレ条件が設定されていなければ,Errorを返す.
////    // --
////    if(mnBndType != BoundaryType::Dirichlet){
////        pLogger->Info(Utility::LoggerMode::Error, "BoundaryType Error,  CBoundaryMesh::distDirichletValue");
////        return;
////    }
////
////
////    if(mMGLevel==0){
////        distDirichletValue_at_CGrid();
////        if(mMaxMGLevel > 0) distDirichletValue_at_FGrid();
////    }else{
////        distDirichletValue_at_FGrid();
////    }
////}




