//
//  BoundaryVolume.cpp
//
//
//
//                          2009.05.13
//                          k.Takeda
#include "EdgeTree.h"
#include "BoundaryParts.h"
#include "Element.h"
#include "BoundaryNode.h"
#include <vector>

#include "Logger.h"

#include "BoundaryVolume.h"
#include "ElementType.h"
#include "ElementProperty.h"
using namespace pmw;

CBoundaryVolume::CBoundaryVolume()
{
    ;
}

CBoundaryVolume::~CBoundaryVolume()
{
    ;
}


// 辺に隣接するVolumeのIDをセット
void CBoundaryVolume::setEdgeNeibVol(const uint& iedge, const uint& neibVolID)
{
    mvEdgeNeibVol[iedge].push_back(neibVolID);
}

// 面に隣接するVolumeのIDをセット
void CBoundaryVolume::setFaceNeibVol(const uint& iface, const uint& neibVolID)
{
    mvFaceNeibVol[iface]= neibVolID;
}

// 辺のBNodeをセット
void CBoundaryVolume::setEdgeBNode(const uint& iedge, CBoundaryNode* pBNode)
{
    mvEdgeBNode[iedge]= pBNode;
}

// 面のBNodeをセット
void CBoundaryVolume::setFaceBNode(const uint& iface, CBoundaryNode* pBNode)
{
    mvFaceBNode[iface]= pBNode;
}

// 体中心のBNodeをセット
void CBoundaryVolume::setVolBNode(CBoundaryNode* pBNode)
{
    mpVolBNode= pBNode;
}


// BNodeにNodeをセット
//// ----
//// *mpElementは,CBounaryPartsが所有
//// ----
//void CBoundaryVolume::setupNode()
//{
//    setupNode_Edge();
//    setupNode_Face();
//    setupNode_Volume();
//}


//// 辺のBNodeにNodeをセット
//// ----
//void CBoundaryVolume::setupNode_Edge()
//{
//    CNode *pEdgeNode, *pNode0, *pNode1;
//    uint elemEdgeID;
//
//    CBoundaryNode *pEdgeBNode, *pBNode0, *pBNode1;
//    uint iedge, *pLocalNum, numOfEdge= getNumOfEdge();
//    for(iedge=0; iedge < numOfEdge; iedge++){
//
//        pLocalNum= getLocalNode_Edge(iedge);// 辺に対応する両端のNodeIndex番号
//
//        pBNode0= mvBNode[pLocalNum[0]];
//        pBNode1= mvBNode[pLocalNum[1]];
//
//        pEdgeBNode= mvEdgeBNode[iedge];
//
//        pNode0= pBNode0->getNode();
//        pNode1= pBNode1->getNode();
//
//        elemEdgeID= mpElement->getEdgeIndex(pNode0, pNode1);
//
//        pEdgeNode= mpElement->getEdgeInterNode(elemEdgeID);
//
//        pEdgeBNode->setNode(pEdgeNode);
//    };
//}


//// 面のBNodeにNodeをセット
//// ----
//void CBoundaryVolume::setupNode_Face()
//{
//    CNode *pFaceNode;
//    CNode *pNode0, *pNode1, *pNode2;
//    uint elemFaceID;
//
//    CBoundaryNode *pFaceBNode;
//    CBoundaryNode *pBNode0, *pBNode1, *pBNode2;
//
//    uint *pLocalNum;
//    uint iface, numOfFace = getNumOfFace();
//    for(iface=0; iface < numOfFace; iface++){
//
//        pLocalNum= getLocalNode_Face(iface);// Faceを構成する頂点番号
//
//        // 頂点番号に対応するBNode => Node を取得
//        pBNode0= mvBNode[pLocalNum[0]]; pBNode1= mvBNode[pLocalNum[1]]; pBNode2= mvBNode[pLocalNum[2]];
//        pNode0= pBNode0->getNode();     pNode1= pBNode1->getNode();     pNode2= pBNode2->getNode();
//
//        // Nodeに対応するMesh-Elementの面番号を取得 => 面Nodeを取得
//        elemFaceID= mpElement->getFaceIndex(pNode0, pNode1, pNode2);
//        pFaceNode = mpElement->getFaceNode(elemFaceID);
//
//        pFaceBNode = mvFaceBNode[iface];
//        pFaceBNode->setNode(pFaceNode);// BNodeにNodeをセット
//    };
//}

//// 体積BNodeにNodeをセット
//void CBoundaryVolume::setupNode_Volume()
//{
//    CNode *pNode = mpElement->getVolumeNode();
//
//    mpVolBNode->setNode(pNode);
//}





// 要素の境界値を分配 => 子Volume
// ----
void CBoundaryVolume::distValue(CBoundaryVolume* pProgVol, const double& coef, const vuint& vDOF)
{
    uint idof, dof;
    double progValue;

    for(idof=0; idof < vDOF.size(); idof++){
        dof = vDOF[idof];
        progValue = mmValue[dof]*coef;

        //pProgVol->addDOF(dof);
        pProgVol->setBndValue(dof,progValue);
    };
}

// ----
// refineの子ルーチン
// ----
// BNodeセット: Hexa分割の子要素
// ----
uint CBoundaryVolume::dividHexa(const uint& iprog, CBoundaryVolume* pProgVol)
{
    switch(iprog){
        case(0):
            // 要素 0
            pProgVol->setBNode(0, mvBNode[0]);     pProgVol->setBNode(1, mvEdgeBNode[0]);
            pProgVol->setBNode(2, mvFaceBNode[0]); pProgVol->setBNode(3, mvEdgeBNode[3]);
            pProgVol->setBNode(4, mvEdgeBNode[8]); pProgVol->setBNode(5, mvFaceBNode[4]);
            pProgVol->setBNode(6, mpVolBNode);     pProgVol->setBNode(7, mvFaceBNode[3]);

            return 0;//子供がぶらさがる場所(BNodeの頂点番号)=＞Elementの子供も頂点にぶらさがるので番号を合わせてある.
            break;
        case(1):
            // 要素 1
            pProgVol->setBNode(0, mvEdgeBNode[0]); pProgVol->setBNode(1, mvBNode[1]);
            pProgVol->setBNode(2, mvEdgeBNode[1]); pProgVol->setBNode(3, mvFaceBNode[0]);
            pProgVol->setBNode(4, mvFaceBNode[4]); pProgVol->setBNode(5, mvEdgeBNode[9]);
            pProgVol->setBNode(6, mvFaceBNode[2]); pProgVol->setBNode(7, mpVolBNode);

            return 1;//子供がぶらさがる場所(BNodeの頂点番号)=＞Elementの子供も頂点にぶらさがるので番号を合わせてある.
            break;
        case(2):
            // 要素 2
            pProgVol->setBNode(0, mvEdgeBNode[8]); pProgVol->setBNode(1, mvFaceBNode[4]);
            pProgVol->setBNode(2, mpVolBNode);     pProgVol->setBNode(3, mvFaceBNode[3]);
            pProgVol->setBNode(4, mvBNode[4]);     pProgVol->setBNode(5, mvEdgeBNode[4]);
            pProgVol->setBNode(6, mvFaceBNode[1]); pProgVol->setBNode(7, mvEdgeBNode[7]);

            return 4;//子供がぶらさがる場所(BNodeの頂点番号)=＞Elementの子供も頂点にぶらさがるので番号を合わせてある.
            break;
        case(3):
            // 要素 3
            pProgVol->setBNode(0, mvFaceBNode[4]); pProgVol->setBNode(1, mvEdgeBNode[9]);
            pProgVol->setBNode(2, mvFaceBNode[2]); pProgVol->setBNode(3, mpVolBNode);
            pProgVol->setBNode(4, mvEdgeBNode[4]); pProgVol->setBNode(5, mvBNode[5]);
            pProgVol->setBNode(6, mvEdgeBNode[5]); pProgVol->setBNode(7, mvFaceBNode[1]);

            return 5;//子供がぶらさがる場所(BNodeの頂点番号)=＞Elementの子供も頂点にぶらさがるので番号を合わせてある.
            break;
        case(4):
            // 要素 4
            pProgVol->setBNode(0, mvEdgeBNode[3]); pProgVol->setBNode(1, mvFaceBNode[0]);
            pProgVol->setBNode(2, mvEdgeBNode[2]); pProgVol->setBNode(3, mvBNode[3]);
            pProgVol->setBNode(4, mvFaceBNode[3]); pProgVol->setBNode(5, mpVolBNode);
            pProgVol->setBNode(6, mvFaceBNode[5]); pProgVol->setBNode(7, mvEdgeBNode[11]);

            return 3;//子供がぶらさがる場所(BNodeの頂点番号)=＞Elementの子供も頂点にぶらさがるので番号を合わせてある.
            break;
        case(5):
            // 要素 5
            pProgVol->setBNode(0, mvFaceBNode[0]); pProgVol->setBNode(1, mvEdgeBNode[1]);
            pProgVol->setBNode(2, mvBNode[2]);     pProgVol->setBNode(3, mvEdgeBNode[2]);
            pProgVol->setBNode(4, mpVolBNode);     pProgVol->setBNode(5, mvFaceBNode[2]);
            pProgVol->setBNode(6, mvEdgeBNode[10]);pProgVol->setBNode(7, mvFaceBNode[5]);

            return 2;//子供がぶらさがる場所(BNodeの頂点番号)=＞Elementの子供も頂点にぶらさがるので番号を合わせてある.
            break;
        case(6):
            // 要素 6
            pProgVol->setBNode(0, mvFaceBNode[3]); pProgVol->setBNode(1, mpVolBNode);
            pProgVol->setBNode(2, mvFaceBNode[5]); pProgVol->setBNode(3, mvEdgeBNode[11]);
            pProgVol->setBNode(4, mvEdgeBNode[7]); pProgVol->setBNode(5, mvFaceBNode[1]);
            pProgVol->setBNode(6, mvEdgeBNode[6]); pProgVol->setBNode(7, mvBNode[7]);

            return 7;//子供がぶらさがる場所(BNodeの頂点番号)=＞Elementの子供も頂点にぶらさがるので番号を合わせてある.
            break;
        case(7):
            // 要素 7
            pProgVol->setBNode(0, mpVolBNode);      pProgVol->setBNode(1, mvFaceBNode[2]);
            pProgVol->setBNode(2, mvEdgeBNode[10]); pProgVol->setBNode(3, mvFaceBNode[5]);
            pProgVol->setBNode(4, mvFaceBNode[1]);  pProgVol->setBNode(5, mvEdgeBNode[5]);
            pProgVol->setBNode(6, mvBNode[6]);      pProgVol->setBNode(7, mvEdgeBNode[6]);

            return 6;//子供がぶらさがる場所(BNodeの頂点番号)=＞Elementの子供も頂点にぶらさがるので番号を合わせてある.
            break;
    }
}
// ----
// BNodeセット: Tetra分割の子要素
// ----
uint CBoundaryVolume::dividTetra(const uint& iprog, CBoundaryVolume* pProgVol)
{
    switch(iprog){
        case(0):
            // 要素 0
            pProgVol->setBNode(0, mvEdgeBNode[2]); pProgVol->setBNode(1, mvBNode[0]);
            pProgVol->setBNode(2, mvEdgeBNode[0]); pProgVol->setBNode(3, mvFaceBNode[0]);
            pProgVol->setBNode(4, mvFaceBNode[3]); pProgVol->setBNode(5, mvEdgeBNode[3]);
            pProgVol->setBNode(6, mvFaceBNode[1]); pProgVol->setBNode(7, mpVolBNode);

            return 0;//子供がぶらさがる場所(BNodeの頂点番号)=＞Elementの子供も頂点にぶらさがるので番号を合わせてある.
            break;
        case(1):
            // 要素 1
            pProgVol->setBNode(0, mvBNode[2]);     pProgVol->setBNode(1, mvEdgeBNode[2]);
            pProgVol->setBNode(2, mvFaceBNode[0]); pProgVol->setBNode(3, mvEdgeBNode[1]);
            pProgVol->setBNode(4, mvEdgeBNode[5]); pProgVol->setBNode(5, mvFaceBNode[3]);
            pProgVol->setBNode(6, mpVolBNode);     pProgVol->setBNode(7, mvFaceBNode[2]);

            return 2;//子供がぶらさがる場所(BNodeの頂点番号)=＞Elementの子供も頂点にぶらさがるので番号を合わせてある.
            break;
        case(2):
            // 要素 2
            pProgVol->setBNode(0, mvFaceBNode[0]); pProgVol->setBNode(1, mvEdgeBNode[0]);
            pProgVol->setBNode(2, mvBNode[1]);     pProgVol->setBNode(3, mvEdgeBNode[1]);
            pProgVol->setBNode(4, mpVolBNode);     pProgVol->setBNode(5, mvFaceBNode[1]);
            pProgVol->setBNode(6, mvEdgeBNode[4]); pProgVol->setBNode(7, mvFaceBNode[2]);

            return 1;//子供がぶらさがる場所(BNodeの頂点番号)=＞Elementの子供も頂点にぶらさがるので番号を合わせてある.
            break;
        case(3):
            // 要素 3
            pProgVol->setBNode(0, mvFaceBNode[3]); pProgVol->setBNode(1, mvEdgeBNode[3]);
            pProgVol->setBNode(2, mvFaceBNode[1]); pProgVol->setBNode(3, mpVolBNode);
            pProgVol->setBNode(4, mvEdgeBNode[5]); pProgVol->setBNode(5, mvBNode[3]);
            pProgVol->setBNode(6, mvEdgeBNode[4]); pProgVol->setBNode(7, mvFaceBNode[2]);

            return 3;//子供がぶらさがる場所(BNodeの頂点番号)=＞Elementの子供も頂点にぶらさがるので番号を合わせてある.
            break;
    }
}
// ----
// BNodeセット: Prism分割の子要素
// ----
uint CBoundaryVolume::dividPrism(const uint& iprog, CBoundaryVolume* pProgVol)
{
    switch(iprog){
        case(0):
            // 要素 0
            pProgVol->setBNode(0, mvBNode[2]);     pProgVol->setBNode(1, mvEdgeBNode[1]);
            pProgVol->setBNode(2, mvFaceBNode[0]); pProgVol->setBNode(3, mvEdgeBNode[2]);
            pProgVol->setBNode(4, mvEdgeBNode[5]); pProgVol->setBNode(5, mvFaceBNode[4]);
            pProgVol->setBNode(6, mpVolBNode);     pProgVol->setBNode(7, mvFaceBNode[3]);

            return 2;//子供がぶらさがる場所(BNodeの頂点番号)=＞Elementの子供も頂点にぶらさがるので番号を合わせてある.
            break;
        case(1):
            // 要素 1
            pProgVol->setBNode(0, mvEdgeBNode[1]); pProgVol->setBNode(1, mvBNode[0]);
            pProgVol->setBNode(2, mvEdgeBNode[0]); pProgVol->setBNode(3, mvFaceBNode[0]);
            pProgVol->setBNode(4, mvFaceBNode[4]); pProgVol->setBNode(5, mvEdgeBNode[3]);
            pProgVol->setBNode(6, mvFaceBNode[2]); pProgVol->setBNode(7, mpVolBNode);

            return 0;//子供がぶらさがる場所(BNodeの頂点番号)=＞Elementの子供も頂点にぶらさがるので番号を合わせてある.
            break;
        case(2):
            // 要素 2
            pProgVol->setBNode(0, mvFaceBNode[0]); pProgVol->setBNode(1, mvEdgeBNode[0]);
            pProgVol->setBNode(2, mvBNode[1]);     pProgVol->setBNode(3, mvEdgeBNode[2]);
            pProgVol->setBNode(4, mpVolBNode);     pProgVol->setBNode(5, mvFaceBNode[2]);
            pProgVol->setBNode(6, mvEdgeBNode[4]); pProgVol->setBNode(7, mvFaceBNode[3]);

            return 1;//子供がぶらさがる場所(BNodeの頂点番号)=＞Elementの子供も頂点にぶらさがるので番号を合わせてある.
            break;
        case(3):
            // 要素 3
            pProgVol->setBNode(0, mvEdgeBNode[5]); pProgVol->setBNode(1, mvFaceBNode[4]);
            pProgVol->setBNode(2, mpVolBNode);     pProgVol->setBNode(3, mvFaceBNode[3]);
            pProgVol->setBNode(4, mvBNode[5]);     pProgVol->setBNode(5, mvEdgeBNode[8]);
            pProgVol->setBNode(6, mvFaceBNode[1]); pProgVol->setBNode(7, mvEdgeBNode[7]);

            return 5;//子供がぶらさがる場所(BNodeの頂点番号)=＞Elementの子供も頂点にぶらさがるので番号を合わせてある.
            break;
        case(4):
            // 要素 4
            pProgVol->setBNode(0, mvFaceBNode[4]); pProgVol->setBNode(1, mvEdgeBNode[3]);
            pProgVol->setBNode(2, mvFaceBNode[2]); pProgVol->setBNode(3, mpVolBNode);
            pProgVol->setBNode(4, mvEdgeBNode[8]); pProgVol->setBNode(5, mvBNode[3]);
            pProgVol->setBNode(6, mvEdgeBNode[6]); pProgVol->setBNode(7, mvFaceBNode[1]);

            return 3;//子供がぶらさがる場所(BNodeの頂点番号)=＞Elementの子供も頂点にぶらさがるので番号を合わせてある.
            break;
        case(5):
            // 要素 5
            pProgVol->setBNode(0, mpVolBNode);     pProgVol->setBNode(1, mvFaceBNode[2]);
            pProgVol->setBNode(2, mvEdgeBNode[4]); pProgVol->setBNode(3, mvFaceBNode[3]);
            pProgVol->setBNode(4, mvFaceBNode[1]); pProgVol->setBNode(5, mvEdgeBNode[6]);
            pProgVol->setBNode(6, mvBNode[4]);     pProgVol->setBNode(7, mvEdgeBNode[7]);

            return 4;//子供がぶらさがる場所(BNodeの頂点番号)=＞Elementの子供も頂点にぶらさがるので番号を合わせてある.
            break;
    }
}

// 四面体の体積：スカラー三重積による計算
// ----
// # BoundaryVolumeを分割した四面体の体積計算:(1/6)*| u dot (v cross w) |
// ----
double CBoundaryVolume::tetraVolume(CNode* pNode0, CNode* pNode1, CNode* pNode2, CNode* pNode3)
{
    double vol(0.0);

    double u[3],v[3],w[3];

    u[0]= pNode1->getX() - pNode0->getX();
    u[1]= pNode1->getY() - pNode0->getY();
    u[2]= pNode1->getZ() - pNode0->getZ();

    v[0]= pNode2->getX() - pNode0->getX();
    v[1]= pNode2->getY() - pNode0->getY();
    v[2]= pNode2->getZ() - pNode0->getZ();

    w[0]= pNode3->getX() - pNode0->getX();
    w[1]= pNode3->getY() - pNode0->getY();
    w[2]= pNode3->getZ() - pNode0->getZ();

    // 外積
    // ---
    double cross[3];
    cross[0] = v[1]*w[2] - v[2]*w[1];
    cross[1] = v[2]*w[0] - v[0]*w[2];
    cross[2] = v[0]*w[1] - v[1]*w[0];

    // 内積
    // ---
    double scalar;
    scalar = u[0]*cross[0] + u[1]*cross[1] + u[2]*cross[2];
    
    if(scalar < 0.0) scalar *= -1.0;

    vol = (1.0/6.0)*scalar;
    
    return vol;
}





















