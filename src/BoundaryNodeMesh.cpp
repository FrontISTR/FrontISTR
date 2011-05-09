//
// BoundaryNodeMesh.cpp
//
//
//
//              2010.04.07
//              k.Takeda
#include "BoundaryNodeMesh.h"
using namespace pmw;

CBoundaryNodeMesh::CBoundaryNodeMesh()
{
    ;
}
CBoundaryNodeMesh::~CBoundaryNodeMesh()
{
    ;
}

// Dirichlet .or. Neumann
// --
void CBoundaryNodeMesh::setBndType(const uint& boundType)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();

    switch(boundType){
        case(BoundaryType::Dirichlet):
            mnBndType= boundType;
            break;
        case(BoundaryType::Neumann):
            mnBndType= boundType;
            break;
        default:
            pLogger->Info(Utility::LoggerMode::Error, "BoundaryType Error, CBoundaryNodeMesh::setType");
            break;
    }
}


void CBoundaryNodeMesh::resizeBNode(const uint& res_size)
{
    mvBNode.resize(res_size);
}

void CBoundaryNodeMesh::setBNode(const uint& index, CBoundarySBNode* pBNode)
{
    // BoundaryNodeID => Index
    //
    uint id= pBNode->getID();
    mmID2Index[id]= index;


    // NodeID => BoundaryNodeID
    //
    CNode *pNode= pBNode->getNode();
    uint node_id= pNode->getID();
    
    mmNodeID2BNodeID[node_id]= pBNode->getID();
}

void CBoundaryNodeMesh::addBNode(CBoundarySBNode* pBNode)
{
    mvBNode.push_back(pBNode);

    // BoundaryNodeID => Index
    //
    uint index= mvBNode.size()-1;
    uint id= pBNode->getID();

    mmID2Index[id]= index;
    
    
    // NodeID => BoundaryNodeID
    //
    CNode *pNode= pBNode->getNode();
    uint node_id= pNode->getID();
    
    mmNodeID2BNodeID[node_id]= pBNode->getID();
}









/*
// ----
// Neumann条件(等価節点力) #=> 面,辺に移行. 節点はRefineに関与しない.
// ----
//void CBoundary...Mesh::EquivalentNodeForce()
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();

    // ノイマン条件が設定されていなければ,Errorを返す.
    // --
    if(mnType != BoundaryType::Neumann){
        pLogger->Info(Utility::LoggerMode::Error, "BoundaryType Error,  CBoundaryNodeMesh::EquivalentNodeForce");
        return;
    }

    CShapeHexa *pHexaShape = CShapeHexa::Instance();
    CShapeTetra *pTetraShape= CShapeTetra::Instance();
    CShapePrism *pPrismShape= CShapePrism::Instance();
    CShapeQuad *pQuadShape = CShapeQuad::Instance();
    CShapeTriangle *pTriShape= CShapeTriangle::Instance();
    CShapeLine *pLineShape= CShapeLine::Instance();

    CNode *pNode;
    CBoundaryNode *pBndNode;
    CBoundaryNodeFace *pEntity;
    double entValue;  //エンティティ(面とか)の所有する境界値
    double integValue;//形状関数の積分値(正規化)
    double nodalValue;
    uint inode;
    uint ient, numOfEnt= mvBoundaryFace.size();

    for(ient=0; ient < numOfEnt; ient++){

        pEntity= mvBoundaryFace[ient];
        entValue= pEntity->getEntityValue();//エンティティの所有する値(面の値とか…)

        switch(pEntity->getType()){
            case(ElementType::Hexa):
                for(inode=0; inode < 8; inode++){
                    // 等価節点力計算
                    integValue= pHexaShape->getIntegralValue8(inode);
                    nodalValue= entValue * integValue;

                    // 境界Node取得
                    pNode= pEntity->getNode(inode);
                    pBndNode= pNode->getBoundaryNode();

                    // 境界値として節点"外力に加算"
                    pBndNode->addNeumannValue(mMGLevel, mnID, nodalValue);
                };
                break;
            case(ElementType::Hexa2):
                for(inode=0; inode < 20; inode++){
                    // 等価節点力計算
                    integValue= pHexaShape->getIntegralValue20(inode);
                    nodalValue= entValue * integValue;

                    // 境界Node取得
                    pNode= pEntity->getNode(inode);
                    pBndNode= pNode->getBoundaryNode();

                    // 境界値として節点"外力に加算"
                    pBndNode->addNeumannValue(mMGLevel, mnID, nodalValue);
                };
                break;
            case(ElementType::Tetra):
                for(inode=0; inode < 4; inode++){
                    // 等価節点力計算
                    integValue= pTetraShape->getIntegValue4(inode);
                    nodalValue= entValue * integValue;

                    // 境界Node取得
                    pNode= pEntity->getNode(inode);
                    pBndNode= pNode->getBoundaryNode();

                    // 境界値として節点"外力に加算"
                    pBndNode->addNeumannValue(mMGLevel, mnID, nodalValue);
                };
                break;
            case(ElementType::Tetra2):
                for(inode=0; inode < 10; inode++){
                    // 等価節点力計算
                    integValue= pTetraShape->getIntegValue10(inode);
                    nodalValue= entValue * integValue;

                    // 境界Node取得
                    pNode= pEntity->getNode(inode);
                    pBndNode= pNode->getBoundaryNode();

                    // 境界値として節点"外力に加算"
                    pBndNode->addNeumannValue(mMGLevel, mnID, nodalValue);
                };
                break;
            case(ElementType::Prism):
                for(inode=0; inode < 6; inode++){
                    // 等価節点力計算
                    integValue= pPrismShape->getIntegValue6(inode);
                    nodalValue= entValue * integValue;

                    // 境界Node取得
                    pNode= pEntity->getNode(inode);
                    pBndNode= pNode->getBoundaryNode();

                    // 境界値として節点"外力に加算"
                    pBndNode->addNeumannValue(mMGLevel, mnID, nodalValue);
                };
                break;
            case(ElementType::Prism2):
                for(inode=0; inode < 15; inode++){
                    // 等価節点力計算
                    integValue= pPrismShape->getIntegValue15(inode);
                    nodalValue= entValue * integValue;

                    // 境界Node取得
                    pNode= pEntity->getNode(inode);
                    pBndNode= pNode->getBoundaryNode();

                    // 境界値として節点"外力に加算"
                    pBndNode->addNeumannValue(mMGLevel, mnID, nodalValue);
                };
                break;
            case(ElementType::Quad):
                for(inode=0; inode < 4; inode++){
                    // 等価節点力計算
                    integValue= pQuadShape->getIntegValue4(inode);
                    nodalValue= entValue * integValue;

                    // 境界Node取得
                    pNode= pEntity->getNode(inode);
                    pBndNode= pNode->getBoundaryNode();

                    // 境界値として節点"外力に加算"
                    pBndNode->addNeumannValue(mMGLevel, mnID, nodalValue);
                };
                break;
            case(ElementType::Quad2):
                for(inode=0; inode < 8; inode++){
                    // 等価節点力計算
                    integValue= pQuadShape->getIntegValue8(inode);
                    nodalValue= entValue * integValue;

                    // 境界Node取得
                    pNode= pEntity->getNode(inode);
                    pBndNode= pNode->getBoundaryNode();

                    // 境界値として節点"外力に加算"
                    pBndNode->addNeumannValue(mMGLevel, mnID, nodalValue);
                };
                break;
            case(ElementType::Triangle):
                for(inode=0; inode < 3; inode++){
                    // 等価節点力計算
                    integValue= pTriShape->getIntegValue3(inode);
                    nodalValue= entValue * integValue;

                    // 境界Node取得
                    pNode= pEntity->getNode(inode);
                    pBndNode= pNode->getBoundaryNode();

                    // 境界値として節点"外力に加算"
                    pBndNode->addNeumannValue(mMGLevel, mnID, nodalValue);
                };
                break;
            case(ElementType::Triangle2):
                for(inode=0; inode < 6; inode++){
                    // 等価節点力計算
                    integValue= pTriShape->getIntegValue6(inode);
                    nodalValue= entValue * integValue;

                    // 境界Node取得
                    pNode= pEntity->getNode(inode);
                    pBndNode= pNode->getBoundaryNode();

                    // 境界値として節点"外力に加算"
                    pBndNode->addNeumannValue(mMGLevel, mnID, nodalValue);
                };
                break;
            case(ElementType::Beam):
                for(inode=0; inode < 2; inode++){
                    // 等価節点力計算
                    integValue= pLineShape->getIntegValue2(inode);
                    nodalValue= entValue * integValue;

                    // 境界Node取得
                    pNode= pEntity->getNode(inode);
                    pBndNode= pNode->getBoundaryNode();

                    // 境界値として節点"外力に加算"
                    pBndNode->addNeumannValue(mMGLevel, mnID, nodalValue);
                };
                break;
            case(ElementType::Beam2):
                for(inode=0; inode < 3; inode++){
                    // 等価節点力計算
                    integValue= pLineShape->getIntegValue3(inode);
                    nodalValue= entValue * integValue;

                    // 境界Node取得
                    pNode= pEntity->getNode(inode);
                    pBndNode= pNode->getBoundaryNode();

                    // 境界値として節点"外力に加算"
                    pBndNode->addNeumannValue(mMGLevel, mnID, nodalValue);
                };
                break;
            default:
                break;
        }
    };
}
*/



