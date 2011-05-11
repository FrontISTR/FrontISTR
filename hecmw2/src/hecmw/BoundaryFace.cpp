//
//  BoundaryFace.cpp
//
//
//
//                      2009.05.18
//                      k.Takeda
#include <vector>

#include "EdgeTree.h"
#include "BoundaryNode.h"
#include "BoundaryParts.h"

#include "BoundaryFace.h"
#include "ElementType.h"
using namespace pmw;

CBoundaryFace::CBoundaryFace()
{
    mArea = 0.0;
    mpFaceBNode= NULL;
}
CBoundaryFace::~CBoundaryFace()
{
    delete []mvbMarkingEdge;
}

// 辺の数
//
uiint CBoundaryFace::getNumOfEdge()
{
    switch(mnShapeType){
        case(ElementType::Quad):case(ElementType::Quad2):
            return NumberOfEdge::Quad();
            break;
        case(ElementType::Triangle):case(ElementType::Triangle2):
            return NumberOfEdge::Triangle();
            break;
        default:
            return NumberOfEdge::Default();
            break;
    }
}

// 1.面の幾何形状のセット
// 2.EdgeNeibFaceの領域確保
// 3.markingの初期化
void CBoundaryFace::setBFaceShape(const uiint& elemType)
{ 
    Utility::CLogger* pLogger = Utility::CLogger::Instance();

    mnShapeType= elemType;
    uiint ivert;
    switch(mnShapeType){
        case(ElementType::Quad):case(ElementType::Quad2):
            mvEdgeNeibFace.resize(4);

            mvbMarkingEdge = new bool[4];
            for(ivert=0; ivert < 4; ivert++)
                mvbMarkingEdge[ivert] = false;//辺集合の為のマーキング初期化

            mvEdgeBNode.resize(4);

            break;
        case(ElementType::Triangle):case(ElementType::Triangle2):
            mvEdgeNeibFace.resize(3);

            mvbMarkingEdge = new bool[3];
            for(ivert=0; ivert < 3; ivert++)
                mvbMarkingEdge[ivert] = false;//辺集合の為のマーキング初期化

            mvEdgeBNode.resize(3);

            break;
        default:
            pLogger->Info(Utility::LoggerMode::Error, "invalid ElementType, CBoundaryFace::setBFaceShape");
            break;
    }

    // 1次-2次 Order
    switch(mnShapeType){
        case(ElementType::Quad):case(ElementType::Triangle):
            mnOrder = ElementOrder::First;
            break;
        case(ElementType::Quad2):case(ElementType::Triangle2):
            mnOrder = ElementOrder::Second;
            break;
    }
}

uiint CBoundaryFace::getNumOfVert()
{
    switch(mnShapeType){
        case(ElementType::Quad):case(ElementType::Quad2):
            return 4;
        case(ElementType::Triangle):case(ElementType::Triangle2):
            return 3;
        default:
            return 0;
    }
}

// 辺-マーキング
// ----
void CBoundaryFace::markingEdge(const uiint& iedge)
{
    mvbMarkingEdge[iedge]= true;
}

// 辺-隣接 面ID のセット(自身は含まない)
// ----
void CBoundaryFace::setEdgeNeibFace(const uiint& iedge, const uiint& neibFaceID)
{
    mvEdgeNeibFace[iedge]= neibFaceID;
}

// EdgeBNodeのセット
// ----
void CBoundaryFace::setEdgeBNode(const uiint& iedge, CBoundaryNode* pEdgeBNode)
{
    mvEdgeBNode[iedge]= pEdgeBNode;
}

// FaceBNodeのセット
// ----
void CBoundaryFace::setFaceBNode(CBoundaryNode* pFaceBNode)
{
    mpFaceBNode= pFaceBNode;
}


// 辺の両端のBNode
// ----
PairBNode CBoundaryFace::getPairBNode(const uiint& iedge)
{
    Utility::CLogger* pLogger = Utility::CLogger::Instance();

    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    uiint *pLocalNum;
    switch(mnShapeType){
        case(ElementType::Quad):case(ElementType::Quad2):
            pLocalNum= pEdgeTree->getQuadLocalNodeNum(iedge);
            break;

        case(ElementType::Triangle):case(ElementType::Triangle2):
            pLocalNum= pEdgeTree->getTriangleLocalNodeNum(iedge);
            break;

        default:
            pLogger->Info(Utility::LoggerMode::Error, "invalid ElementType, CBoundaryFace::getPairBNode");
            break;
    }

    uiint index1st = pLocalNum[0];
    uiint index2nd = pLocalNum[1];

    PairBNode pairBNode;
    pairBNode.first = mvBNode[index1st];
    pairBNode.second= mvBNode[index2nd];

    return pairBNode;
}

// 辺番号
// ----
uiint& CBoundaryFace::getEdgeID(PairBNode& pairBNode)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    
    uiint ivert= getVertIndex(pairBNode.first);
    uiint jvert= getVertIndex(pairBNode.second);
    
    switch(mnShapeType){
        case(ElementType::Quad):case(ElementType::Quad2):
            return pEdgeTree->getQuadEdgeIndex(ivert,jvert);
            break;
        case(ElementType::Triangle):case(ElementType::Triangle2):
            return pEdgeTree->getTriangleEdgeIndex(ivert, jvert);
            break;
        default:
            return pEdgeTree->getDisagTypeEdgeIndex(ivert, jvert);
            break;
    }
}


// ----
// EdgeBNodeにNodeをセット : *mpElementは,CBounaryPartsが所有
// ----
void CBoundaryFace::setupNode_Edge()
{
    uiint nNumOfEdge = getNumOfEdge();
    // 辺のBNodeにNodeをセット
    // ----
    CNode *pEdgeNode, *pNode0, *pNode1;
    uiint elemEdgeID;

    CBoundaryNode *pEdgeBNode, *pBNode0, *pBNode1;
    uiint iedge;

    for(iedge=0; iedge < nNumOfEdge; iedge++){
        pEdgeBNode= mvEdgeBNode[iedge];

        if(iedge != nNumOfEdge-1){
            pBNode0= mvBNode[iedge]; pBNode1= mvBNode[iedge+1];// 0-1, 1-2, 2-3 の辺
        }else{
            pBNode0= mvBNode[iedge]; pBNode1= mvBNode[0];// 3-0 の辺
        }
        pNode0= pBNode0->getNode(); pNode1= pBNode1->getNode();

        elemEdgeID= mpElement->getEdgeIndex(pNode0, pNode1);//要素の辺番号 取得
        pEdgeNode = mpElement->getEdgeInterNode(elemEdgeID);//要素辺のNode 取得

        pEdgeBNode->setNode(pEdgeNode);//辺のBNodeにNodeをセット
    };
}
// 
// ----
// FaceBNodeにNodeをセット  *mpElementは,CBounaryPartsが所有
// ----
void CBoundaryFace::setupNode_Face()
{
    // 面のBNodeにNodeをセット
    // ----
    CNode *pFaceNode;
    pFaceNode= mpElement->getFaceNode(mnElemFaceID);

    mpFaceBNode->setNode(pFaceNode);
}

//// ----
//void CBoundaryFace::setupNode()
//{
//    Utility::CLogger* pLogger = Utility::CLogger::Instance();
//
//    uint numOfEdge;
//    switch(mnShapeType){
//        case(ElementType::Quad):case(ElementType::Quad2):
//            numOfEdge= NumberOfEdge::Quad();
//            break;
//        case(ElementType::Triangle):case(ElementType::Triangle2):
//            numOfEdge= NumberOfEdge::Triangle();
//            break;
//        default:
//            pLogger->Info(Utility::LoggerMode::Error, "invalid ElementType, CBoundaryFace::setupNode");
//            numOfEdge= NumberOfEdge::Default();
//            break;
//    }
//
//    // 辺のBNodeにNodeをセット
//    // ----
//    CNode *pEdgeNode, *pNode0, *pNode1;
//    uint elemEdgeID;
//
//    CBoundaryNode *pEdgeBNode, *pBNode0, *pBNode1;
//    uint iedge;
//
//    for(iedge=0; iedge < numOfEdge; iedge++){
//        pEdgeBNode= mvEdgeBNode[iedge];
//
//        if(iedge != numOfEdge-1){
//            pBNode0= mvBNode[iedge]; pBNode1= mvBNode[iedge+1];// 0-1, 1-2, 2-3 の辺
//        }else{
//            pBNode0= mvBNode[iedge]; pBNode1= mvBNode[0];// 3-0 の辺
//        }
//        pNode0= pBNode0->getNode(); pNode1= pBNode1->getNode();
//
//        elemEdgeID= mpElement->getEdgeIndex(pNode0, pNode1);//要素の辺番号 取得
//        pEdgeNode = mpElement->getEdgeInterNode(elemEdgeID);//要素辺のNode 取得
//
//        pEdgeBNode->setNode(pEdgeNode);//辺のBNodeにNodeをセット
//    };
//
//
//    // 面のBNodeにNodeをセット
//    // ----
//    CNode *pFaceNode;
//    pFaceNode= mpElement->getFaceNode(mnElemFaceID);
//
//    mpFaceBNode->setNode(pFaceNode);
//}


// ----
// 再分割 : countID==新FaceのID
// ----
void CBoundaryFace::refine(uiint& countID, const vuint& vDOF)
{
    Utility::CLogger* pLogger = Utility::CLogger::Instance();

    uiint numOfProg;
    uiint iprog;
    CBoundaryFace *pProgBFace;
    double progArea;
    double progValue;
    double coef; //面積比
    uiint idof, dof;//dofインデックスと,DOF番号
    
    CElement  *pProgElem;// 子要素
    uiint     nProgEntity;// 子要素Entity番号
    uiint         ipos;   // 親要素の頂点番号(子要素の位置)
    CNode      *pNode;   // 親要素の頂点のNode
    vector<CNode*>  vNode;// progBFaceの頂点のNode
    CBoundaryNode *pBNode;// progBFaceのBNode
    uiint           ibnode;// progBFaceのindex

    switch(mnShapeType){
        //----
        //Quad
        //----
        case(ElementType::Quad):case(ElementType::Quad2):
            numOfProg= 4;
            mvProgBFace.reserve(numOfProg);

            for(iprog=0; iprog < numOfProg; iprog++){

                pProgBFace = new CBoundaryFace;// <<<<<<<<<<<<<<<< new

                pProgBFace->setID(countID);
                countID++;
                if(mnOrder==ElementOrder::First ) pProgBFace->setBFaceShape(ElementType::Quad);
                if(mnOrder==ElementOrder::Second) pProgBFace->setBFaceShape(ElementType::Quad2);
                pProgBFace->resizeBNode(4);

                // BNodeのセット
                switch(iprog){
                    case(0):
                        pProgBFace->setBNode(0, mvBNode[0]);
                        pProgBFace->setBNode(1, mvEdgeBNode[0]);
                        pProgBFace->setBNode(2, mpFaceBNode);
                        pProgBFace->setBNode(3, mvEdgeBNode[3]);
                        
                        progArea= pProgBFace->calcArea();

                        break;
                    case(1):
                        pProgBFace->setBNode(0, mvEdgeBNode[0]);
                        pProgBFace->setBNode(1, mvBNode[1]);
                        pProgBFace->setBNode(2, mvEdgeBNode[1]);
                        pProgBFace->setBNode(3, mpFaceBNode);

                        progArea= pProgBFace->calcArea();

                        break;
                    case(2):
                        pProgBFace->setBNode(0, mvEdgeBNode[1]);
                        pProgBFace->setBNode(1, mvBNode[2]);
                        pProgBFace->setBNode(2, mvEdgeBNode[2]);
                        pProgBFace->setBNode(3, mpFaceBNode);

                        progArea= pProgBFace->calcArea();

                        break;
                    case(3):
                        pProgBFace->setBNode(0, mvEdgeBNode[2]);
                        pProgBFace->setBNode(1, mvBNode[3]);
                        pProgBFace->setBNode(2, mvEdgeBNode[3]);
                        pProgBFace->setBNode(3, mpFaceBNode);

                        progArea= pProgBFace->calcArea();

                        break;
                }
                // ProgBFaceコンテナ
                mvProgBFace.push_back(pProgBFace);

                //----
                // 面積比による境界値再配分
                //----
                coef = progArea / mArea;

                for(idof=0; idof < vDOF.size(); idof++){
                    dof = vDOF[idof];
                    progValue = mmValue[dof]*coef;

                    //pProgBFace->addDOF(dof);
                    pProgBFace->setBndValue(dof,progValue);
                };

                //----
                // progElem, 要素ID のセット
                //----
                pNode = mvBNode[iprog]->getNode();

                ipos  = mpElement->getLocalVertNum(pNode->getID());
                pProgElem = mpElement->getProgElem(ipos);
                
                pProgBFace->setElement(pProgElem);
                pProgBFace->setElementID(pProgElem->getID());
                
                //----
                // progElemのFace番号(ElemFaceID:エンティティ番号)のセット
                //----
                vNode.clear();
                for(ibnode=0; ibnode < 3; ibnode++){
                    pBNode = pProgBFace->getBNode(ibnode);
                    pNode = pBNode->getNode();
                    
                    vNode.push_back(pNode);
                };
                nProgEntity = pProgElem->getFaceIndex(vNode[0], vNode[1], vNode[2]);
                pProgBFace->setElementFaceID(nProgEntity);
            };
            
            break;
        //----
        //Triangle
        //----
        case(ElementType::Triangle):case(ElementType::Triangle2):
            numOfProg= 3;
            mvProgBFace.reserve(numOfProg);

            for(iprog=0; iprog < numOfProg; iprog++){

                pProgBFace = new CBoundaryFace;

                pProgBFace->setID(countID);
                countID++;
                //親が三角形であっても,子はQuad形状
                if(mnOrder==ElementOrder::First ) pProgBFace->setBFaceShape(ElementType::Quad);
                if(mnOrder==ElementOrder::Second) pProgBFace->setBFaceShape(ElementType::Quad2);
                pProgBFace->resizeBNode(4);

                // BNodeのセット
                switch(iprog){
                    case(0):
                        pProgBFace->setBNode(0, mvBNode[0]);
                        pProgBFace->setBNode(1, mvEdgeBNode[0]);
                        pProgBFace->setBNode(2, mpFaceBNode);
                        pProgBFace->setBNode(3, mvEdgeBNode[2]);

                        progArea= pProgBFace->calcArea();
                        break;
                    case(1):
                        pProgBFace->setBNode(0, mvEdgeBNode[0]);
                        pProgBFace->setBNode(1, mvBNode[1]);
                        pProgBFace->setBNode(2, mvEdgeBNode[1]);
                        pProgBFace->setBNode(3, mpFaceBNode);

                        progArea= pProgBFace->calcArea();
                        break;
                    case(2):
                        pProgBFace->setBNode(0, mvEdgeBNode[1]);
                        pProgBFace->setBNode(1, mvBNode[2]);
                        pProgBFace->setBNode(2, mvEdgeBNode[2]);
                        pProgBFace->setBNode(3, mpFaceBNode);

                        progArea= pProgBFace->calcArea();
                        break;
                }
                // ProgBFaceコンテナ
                mvProgBFace.push_back(pProgBFace);

                //----
                // 面積比による境界値再配分
                //----
                coef = progArea / mArea;

                for(idof=0; idof < vDOF.size(); idof++){
                    dof = vDOF[idof];
                    progValue = mmValue[dof]*coef;

                    //pProgBFace->addDOF(dof);
                    pProgBFace->setBndValue(dof,progValue);
                };

                //----
                // progElem, 要素ID のセット
                //----
                pNode = mvBNode[iprog]->getNode();
                ipos  = mpElement->getLocalVertNum(pNode->getID());
                pProgElem = mpElement->getProgElem(ipos);

                pProgBFace->setElement(pProgElem);
                pProgBFace->setElementID(pProgElem->getID());

                //----
                // progElemのFace番号(ElemFaceID:エンティティ番号)のセット
                //----
                vNode.clear();
                for(ibnode=0; ibnode < 3; ibnode++){
                    pBNode = pProgBFace->getBNode(ibnode);
                    pNode = pBNode->getNode();

                    vNode.push_back(pNode);
                };
                nProgEntity = pProgElem->getFaceIndex(vNode[0], vNode[1], vNode[2]);
                pProgBFace->setElementFaceID(nProgEntity);
            };
            break;
            
        default:
            pLogger->Info(Utility::LoggerMode::Error, "invalid ElementType, CBoundaryFace::refine");
            break;
    }
    
}


// 外積による三角形面積
double CBoundaryFace::triArea(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    // 外積: 絶対値==三角形面積の2倍
    // ----
    // x = y1 z2 - z1 y2
    // y = z1 x2 - x1 z2
    // z = x1 y2 - y1 x2
    double x, y, z, x1, y1, z1, x2, y2, z2;

    // vector 0
    x1 = pNode1->getX() - pNode0->getX();
    y1 = pNode1->getY() - pNode0->getY();
    z1 = pNode1->getZ() - pNode0->getZ();

    // vector 1
    x2 = pNode2->getX() - pNode0->getX();
    y2 = pNode2->getY() - pNode0->getY();
    z2 = pNode2->getZ() - pNode0->getZ();

    x = y1*z2 - z1*y2;
    y = z1*x2 - x1*z2;
    z = x1*y2 - y1*x2;

    return sqrt(x*x + y*y + z*z)*0.5;
}

// Faceの面積計算(境界値再配分の為)
// ----
double& CBoundaryFace::calcArea()
{
    Utility::CLogger* pLogger = Utility::CLogger::Instance();
    
    CNode *pNode0, *pNode1, *pNode2;

    switch(mnShapeType){
        case(ElementType::Quad):case(ElementType::Quad2):
            
            mArea = 0.0;
            
            pNode0 = mvBNode[0]->getNode();
            pNode1 = mvBNode[1]->getNode();
            pNode2 = mvBNode[2]->getNode();

            mArea += triArea(pNode0, pNode1, pNode2);
            
            pNode0 = mvBNode[2]->getNode();
            pNode1 = mvBNode[3]->getNode();
            pNode2 = mvBNode[0]->getNode();

            mArea += triArea(pNode0, pNode1, pNode2);
            
            break;
            
        case(ElementType::Triangle):case(ElementType::Triangle2):

            pNode0 = mvBNode[0]->getNode();
            pNode1 = mvBNode[1]->getNode();
            pNode2 = mvBNode[2]->getNode();

            mArea = triArea(pNode0, pNode1, pNode2);

            break;
        default:
            pLogger->Info(Utility::LoggerMode::Error, "invalid ElementType, CBoundaryFace::Area");
            break;
    }
    
    return mArea;
}

// 上位GridのBNodeへディレクレ値を配分
//
void CBoundaryFace::distDirichletVal(const uiint& dof, const uiint& mgLevel, const uiint& nMaxMGLevel)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();

    CEdgeTree *pEdgeTree = CEdgeTree::Instance();
    uiint *pnEdgeVert;
    uiint iedge, numOfEdge;
    double dAveVal;

    // 辺のDirichlet値
    //
    switch(mnShapeType){
        case(ElementType::Quad):
        case(ElementType::Quad2):
            numOfEdge=4;
            for(iedge=0; iedge < numOfEdge; iedge++){
                pnEdgeVert = pEdgeTree->getQuadLocalNodeNum(iedge);

                dAveVal=0.0;
                dAveVal += mvBNode[pnEdgeVert[0]]->getValue(dof, mgLevel);
                dAveVal += mvBNode[pnEdgeVert[1]]->getValue(dof, mgLevel);

                dAveVal *= 0.5;

                if(mnShapeType==ElementType::Quad){
                    if(mgLevel!=nMaxMGLevel)
                      mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, dAveVal);//上位グリッドのディレクレ値
                }
                if(mnShapeType==ElementType::Quad2){
                    mvEdgeBNode[iedge]->setValue(dof, mgLevel,   dAveVal);//カレント・グリッドのディレクレ値
                    if(mgLevel!=nMaxMGLevel)
                      mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, dAveVal);//上位グリッドのディレクレ値
                }
            };
            break;

        case(ElementType::Triangle):
        case(ElementType::Triangle2):
            numOfEdge=3;
            for(iedge=0; iedge < numOfEdge; iedge++){
                pnEdgeVert = pEdgeTree->getTriangleLocalNodeNum(iedge);//QuadからTriangleに変更'10.10.07

                dAveVal=0.0;
                dAveVal += mvBNode[pnEdgeVert[0]]->getValue(dof, mgLevel);
                dAveVal += mvBNode[pnEdgeVert[1]]->getValue(dof, mgLevel);

                dAveVal *= 0.5;

                if(mnShapeType==ElementType::Triangle){
                    if(mgLevel!=nMaxMGLevel)
                      mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, dAveVal);//上位グリッドのディレクレ値
                }

                if(mnShapeType==ElementType::Triangle2){
                    mvEdgeBNode[iedge]->setValue(dof, mgLevel,   dAveVal);//カレント・グリッドのディレクレ値
                    if(mgLevel!=nMaxMGLevel)
                      mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, dAveVal);//上位グリッドのディレクレ値
                }
            };
            break;
            
        default:
            pLogger->Info(Utility::LoggerMode::Error, "BoundaryFace::distDirichletVal, invalid ElementType");
            break;
    }

    //面のDirichlet値 & ここで頂点になっているBNodeの上位GridへのDirichlet値
    //
    uiint ivert, numOfVert=getNumOfVert();
    double dVal;
    dAveVal=0.0;
    for(ivert=0; ivert < numOfVert; ivert++){
        dVal = mvBNode[ivert]->getValue(dof, mgLevel);
        dAveVal += dVal;

        if(mgLevel!=nMaxMGLevel)
          mvBNode[ivert]->setValue(dof, mgLevel+1, dVal);//上位Gridへ下位の値をそのまま渡す:ディレクレ値(頂点)
    };
    //面中心
    dAveVal /= (double)numOfVert;
    if(mgLevel!=nMaxMGLevel)
      mpFaceBNode->setValue(dof, mgLevel+1, dAveVal);//上位Gridへディレクレ値(面中心)
    
////    if(mgLevel > 0){
////        dAveVal /= (double)numOfVert;
////
////        if(mgLevel!=nMaxMGLevel)
////          mpFaceBNode->setValue(dof, mgLevel+1, dAveVal);//上位Gridへディレクレ値(面中心)
////    }
////    if(mgLevel==0){
////        if(mgLevel!=nMaxMGLevel)
////          mpFaceBNode->setValue(dof, mgLevel+1, mmValue[dof]);//Level==0の場合は、要素境界値をそのまま渡す(BNode自体は上位Grid)
////    }
}


// 2次要素の場合、辺BNodeをmvBNodeに移設.
//
void CBoundaryFace::replaceEdgeBNode()
{
    // 四辺形2次
    if(mnShapeType==ElementType::Quad2){
        mvBNode.resize(NumberOfNode::Quad2());

        uiint nNumOfVert=NumberOfVertex::Quad();
        uiint nNumOfEdge=NumberOfEdge::Quad();
        uiint iedge;

        for(iedge=0; iedge < nNumOfEdge; iedge++){
            CBoundaryNode *pBNode = mvEdgeBNode[iedge];
            mvBNode[nNumOfVert + iedge] = pBNode;
        };
    }
    // 三角形2次
    if(mnShapeType==ElementType::Triangle2){
        mvBNode.resize(NumberOfNode::Triangle2());

        uiint nNumOfVert=NumberOfVertex::Triangle();
        uiint nNumOfEdge=NumberOfEdge::Triangle();
        uiint iedge;

        for(iedge=0; iedge < nNumOfEdge; iedge++){
            CBoundaryNode *pBNode = mvEdgeBNode[iedge];
            mvBNode[nNumOfVert + iedge] = pBNode;
        };
    }
}

// Refine時の辺BNodeの解放
//
void CBoundaryFace::deleteProgData()
{
    // 2次要素以外の場合、辺ノード解放
    if(mnOrder==ElementOrder::First || mnOrder==ElementOrder::Zero){
        vector<CBoundaryNode*>().swap(mvEdgeBNode);
    }
}




