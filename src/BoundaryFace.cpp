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
}
CBoundaryFace::~CBoundaryFace()
{
//    //debug
//    cout << "~CBoundaryFace" << endl;
}

// 辺の数
//
uint CBoundaryFace::getNumOfEdge()
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
void CBoundaryFace::setBFaceShape(const uint& elemType)
{ 
    mnShapeType= elemType;
    uint ivert;
    switch(mnShapeType){
        case(ElementType::Quad):case(ElementType::Quad2):
            mvEdgeNeibFace.resize(4);

            mvbMarkingEdge.reserve(4);
            for(ivert=0; ivert < 4; ivert++)
                mvbMarkingEdge.push_back(false);//辺集合の為のマーキング初期化

            mvEdgeBNode.resize(4);

            break;
        case(ElementType::Triangle):case(ElementType::Triangle2):
            mvEdgeNeibFace.resize(3);

            mvbMarkingEdge.reserve(3);
            for(ivert=0; ivert < 3; ivert++)
                mvbMarkingEdge.push_back(false);//辺集合の為のマーキング初期化

            mvEdgeBNode.resize(3);

            break;
        default:
            mpLogger->Info(Utility::LoggerMode::Error, "invalid ElementType, CBoundaryFace::setBFaceShape");
            break;
    }
}

// 辺-マーキング
// ----
void CBoundaryFace::markingEdge(const uint& iedge)
{
    mvbMarkingEdge[iedge]= true;
}

// 辺-隣接 面ID のセット(自身は含まない)
// ----
void CBoundaryFace::setEdgeNeibFace(const uint& iedge, const uint& neibFaceID)
{
    mvEdgeNeibFace[iedge]= neibFaceID;
}

// EdgeBNodeのセット
// ----
void CBoundaryFace::setEdgeBNode(const uint& iedge, CBoundaryNode* pEdgeBNode)
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
PairBNode& CBoundaryFace::getPairBNode(const uint& iedge)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    uint *pLocalNum;
    switch(mnShapeType){
        case(ElementType::Quad):case(ElementType::Quad2):
            pLocalNum= pEdgeTree->getQuadLocalNodeNum(iedge);
            break;

        case(ElementType::Triangle):case(ElementType::Triangle2):
            pLocalNum= pEdgeTree->getTriangleLocalNodeNum(iedge);
            break;

        default:
            mpLogger->Info(Utility::LoggerMode::Error, "invalid ElementType, CBoundaryFace::getPairBNode");
            break;
    }

    uint index1st = pLocalNum[0];
    uint index2nd = pLocalNum[1];
    
    mPairBNode.first = mvBNode[index1st];
    mPairBNode.second= mvBNode[index2nd];

    return mPairBNode;
}

// 辺番号
// ----
uint& CBoundaryFace::getEdgeID(PairBNode& pairBNode)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    
    uint ivert= getVertIndex(pairBNode.first);
    uint jvert= getVertIndex(pairBNode.second);

    ////debug
    ////-----
    ////vert index チェック
    //uint index;
    //for(index=0; index < mvBNode.size(); index++){
    //    if(mvBNode[index]->getID()==pairBNode.first->getID()){
    //        if(ivert!=index) cout << "CBoundaryFace::getEdgeID, getVertIndex Error" << endl;
    //    }
    //    if(mvBNode[index]->getID()==pairBNode.second->getID()){
    //        if(jvert!=index) cout << "CBoundaryFace::getEdgeID, getVertIndex Error" << endl;
    //    }
    //};
    
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

// BNodeにNodeをセット
// ----
// *mpElementは,CBounaryPartsが所有
// ----
void CBoundaryFace::setupNode()
{
    uint numOfEdge;
    switch(mnShapeType){
        case(ElementType::Quad):case(ElementType::Quad2):
            numOfEdge= NumberOfEdge::Quad();
            break;
        case(ElementType::Triangle):case(ElementType::Triangle2):
            numOfEdge= NumberOfEdge::Triangle();
            break;
        default:
            mpLogger->Info(Utility::LoggerMode::Error, "invalid ElementType, CBoundaryFace::setupNode");
            numOfEdge= NumberOfEdge::Default();
            break;
    }

    // 辺のBNodeにNodeをセット
    // ----
    CNode *pEdgeNode, *pNode0, *pNode1;
    uint elemEdgeID;

    CBoundaryNode *pEdgeBNode, *pBNode0, *pBNode1;
    uint iedge;

    for(iedge=0; iedge < numOfEdge; iedge++){
        pEdgeBNode= mvEdgeBNode[iedge];

        if(iedge != numOfEdge-1){
            pBNode0= mvBNode[iedge]; pBNode1= mvBNode[iedge+1];// 0-1, 1-2, 2-3 の辺
        }else{
            pBNode0= mvBNode[iedge]; pBNode1= mvBNode[0];// 3-0 の辺
        }
        pNode0= pBNode0->getNode(); pNode1= pBNode1->getNode();

        elemEdgeID= mpElement->getEdgeIndex(pNode0, pNode1);//要素の辺番号 取得
        pEdgeNode = mpElement->getEdgeInterNode(elemEdgeID);//要素辺のNode 取得

        pEdgeBNode->setNode(pEdgeNode);//辺のBNodeにNodeをセット
    };


    // 面のBNodeにNodeをセット
    // ----
    CNode *pFaceNode;
    pFaceNode= mpElement->getFaceNode(mnElemFaceID);

    mpFaceBNode->setNode(pFaceNode);
}


// ----
// 再分割 : countID==新FaceのID
// ----
void CBoundaryFace::refine(uint& countID, const vuint& vDOF)
{
    uint numOfProg;
    uint iprog;
    CBoundaryFace *pProgBFace;
    double progArea;
    double progValue;
    double coef; //面積比
    uint idof, dof;//dofインデックスと,DOF番号
    
    CElement  *pProgElem;// 子要素
    uint     nProgEntity;// 子要素Entity番号
    uint         ipos;   // 親要素の頂点番号(子要素の位置)
    CNode      *pNode;   // 親要素の頂点のNode
    vector<CNode*>  vNode;// progBFaceの頂点のNode
    CBoundaryNode *pBNode;// progBFaceのBNode
    uint           ibnode;// progBFaceのindex

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
                pProgBFace->setBFaceShape(ElementType::Quad);
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
                pProgBFace->setBFaceShape(ElementType::Quad);//親が三角形であっても,子はQuad
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
            mpLogger->Info(Utility::LoggerMode::Error, "invalid ElementType, CBoundaryFace::refine");
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
            mpLogger->Info(Utility::LoggerMode::Error, "invalid ElementType, CBoundaryFace::Area");
            break;
    }
    
    return mArea;
}

// 上位GridのBNodeへディレクレ値を配分
//
void CBoundaryFace::distDirichletVal(const uint& dof, const uint& mgLevel)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();

    CEdgeTree *pEdgeTree = CEdgeTree::Instance();
    uint *pnEdgeVert;
    uint iedge, numOfEdge;
    double dAveVal;

    // 辺のDirichlet値
    //
    switch(mnShapeType){
        case(ElementType::Quad):
        case(ElementType::Quad2):// !!!! ここのGridでの辺BNodeのディレクレ値はセットしない !!!
            numOfEdge=4;
            for(iedge=0; iedge < numOfEdge; iedge++){
                pnEdgeVert = pEdgeTree->getQuadLocalNodeNum(iedge);

                dAveVal=0.0;
                dAveVal += mvBNode[pnEdgeVert[0]]->getValue(dof, mgLevel);
                dAveVal += mvBNode[pnEdgeVert[1]]->getValue(dof, mgLevel);

                dAveVal *= 0.5;

                mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, dAveVal);//上位グリッドのディレクレ値
            };
            break;
        case(ElementType::Triangle):
        case(ElementType::Triangle2):// !!!! ここのGridでの辺BNodeのディレクレ値はセットしない !!!
            numOfEdge=3;
            for(iedge=0; iedge < numOfEdge; iedge++){
                pnEdgeVert = pEdgeTree->getQuadLocalNodeNum(iedge);

                dAveVal=0.0;
                dAveVal += mvBNode[pnEdgeVert[0]]->getValue(dof, mgLevel);
                dAveVal += mvBNode[pnEdgeVert[1]]->getValue(dof, mgLevel);

                dAveVal *= 0.5;

                mvEdgeBNode[iedge]->setValue(dof, mgLevel+1, dAveVal);//上位グリッドのディレクレ値
            };
            break;
        default:
            pLogger->Info(Utility::LoggerMode::Error, "BoundaryFace::distDirichletVal, invalid ElementType");
            break;
    }
    
    //面のDirichlet値 & ここで頂点になっているBNodeの上位GridへのDirichlet値
    //
    uint ivert, numOfVert=mvBNode.size();
    double dVal;
    dAveVal=0.0;
    for(ivert=0; ivert < numOfVert; ivert++){
        dVal = mvBNode[ivert]->getValue(dof, mgLevel);
        dAveVal += dVal;

        mvBNode[ivert]->setValue(dof, mgLevel+1, dVal);//上位Gridへ下位の値をそのまま渡す:ディレクレ値(頂点)
    };
    if(mgLevel > 0){
        dAveVal /= (double)numOfVert;
        mpFaceBNode->setValue(dof, mgLevel+1, dAveVal);//上位Gridへディレクレ値(面中心)
    }
    if(mgLevel==0){
        mpFaceBNode->setValue(dof, mgLevel+1, mmValue[dof]);//Level==0の場合は、要素境界値をそのまま渡す(BNode自体は上位Grid)
    }
}








