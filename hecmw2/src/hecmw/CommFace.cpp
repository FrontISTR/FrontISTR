//
//  CommFace.cpp
//
//
//
//              2010.03.02
//              k.Takeda
#include <vector>
#include "Element.h"

#include "CommFace.h"
#include "ElementType.h"
using namespace pmw;

CCommFace::CCommFace()
{
    ;
}
CCommFace::~CCommFace()
{
    ;
}
// 初期化
void CCommFace::initialize(const uiint& numOfVert, const uiint& numOfEdge, const uiint& nOrder)
{
    // CommNodeサイズ、FaceTypeの設定
    switch(numOfVert){
        case(4):
            if(ElementOrder::First== nOrder){ mFaceType= ElementType::Quad;  mvCommNode.resize(4);}
            if(ElementOrder::Second==nOrder){ mFaceType= ElementType::Quad2; mvCommNode.resize(8);}
            break;
        case(3):
            if(ElementOrder::First== nOrder){ mFaceType= ElementType::Triangle;  mvCommNode.resize(3);}
            if(ElementOrder::Second==nOrder){ mFaceType= ElementType::Triangle2; mvCommNode.resize(6);}
            break;
        case(2):
            if(ElementOrder::First== nOrder){ mFaceType= ElementType::Beam;  mvCommNode.resize(2);}
            if(ElementOrder::Second==nOrder){ mFaceType= ElementType::Beam2; mvCommNode.resize(3);}
            break;
        default:
            break;
    }
    mvEdgeCommNode.resize(numOfEdge);//辺ノード
    
    mvEdgeCommFace.resize(numOfEdge);//辺に隣接する”面”

    mNumOfEdge= numOfEdge;

    //辺マーキングの初期化
    mvbEdgeMarking = new bool[numOfEdge];
    uiint iedge;
    for(iedge=0; iedge< numOfEdge; iedge++){
        mvbEdgeMarking[iedge]= false;
    };
}

uiint CCommFace::getOrder()
{
    switch(mFaceType){
        case(ElementType::Quad):
        case(ElementType::Triangle):
        case(ElementType::Beam):
            return ElementOrder::First;

        case(ElementType::Quad2):
        case(ElementType::Triangle2):
        case(ElementType::Beam2):
            return ElementOrder::Second;
    }
}
uiint CCommFace::getNumOfVert()
{
    switch(mFaceType){
        case(ElementType::Quad):case(ElementType::Quad2):
            return 4;
        case(ElementType::Triangle):case(ElementType::Triangle2):
            return 3;
        case(ElementType::Beam):case(ElementType::Beam2):
            return 2;
    }
}


uiint& CCommFace::getEdgeIndex(PairCommNode& pairCommNode)
{
    // Node0, Node1のインデックス番号から辺番号を特定
    uiint first_id=  pairCommNode.first->getID();
    uiint second_id= pairCommNode.second->getID();
    uiint self_id;
    uiint localNum0, localNum1;
    uiint numOfVert;
    uiint inode;

    ////debug
    //cout << "CommFace::getEdgeIndex(pairCommNode), numOfNode= " << numOfNode << endl;
    if(mNumOfEdge==3) numOfVert = 3;
    if(mNumOfEdge==4) numOfVert = 4;
    if(mNumOfEdge==1) numOfVert = 2;

    for(inode=0; inode< numOfVert; inode++){
        self_id=  mvCommNode[inode]->getID();

        ////debug
        //cout << "CCommFace::getEdgeIndex(pairCommNode), self_id= " << self_id << endl;

        if(self_id==first_id)  localNum0= inode;
        if(self_id==second_id) localNum1= inode;

        //if(mvCommNode[inode]->getID() == pairCommNode.first->getID() ) localNum0=inode;
        //if(mvCommNode[inode]->getID() == pairCommNode.second->getID()) localNum1=inode;
    };

    ////debug
    //cout << "CCommFace::getEdgeIndex(pairCommNode), first_id= " << first_id << ", second_id= " << second_id << endl;
    //cout << "CCommFace::getEdgeIndex(pairCommNode), localNum0= " << localNum0 << ", localNum1= " << localNum1 << endl;
    
    
    CEdgeTree* pEdgeTree= CEdgeTree::Instance();// 局所番号の組み合わせで辺番号が特定
    
    switch(mFaceType){
        case(ElementType::Quad):case(ElementType::Quad2):
            return pEdgeTree->getQuadEdgeIndex(localNum0, localNum1);

        case(ElementType::Triangle):case(ElementType::Triangle2):
            return pEdgeTree->getTriangleEdgeIndex(localNum0, localNum1);

        case(ElementType::Beam):case(ElementType::Beam2):
            return pEdgeTree->getBeamEdgeIndex(localNum0, localNum1);
            
        default:
            return pEdgeTree->getDisagTypeEdgeIndex(localNum0, localNum1);//999が返る.
    }
}

PairCommNode CCommFace::getEdgePairCommNode(const uiint& iedge)
{
    Utility::CLogger *pLogger;

    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    uiint* vertNum;
    switch(mFaceType){
        case(ElementType::Quad):case(ElementType::Quad2):
            vertNum= pEdgeTree->getQuadLocalNodeNum(iedge);
            break;
        case(ElementType::Triangle):case(ElementType::Triangle2):
            vertNum= pEdgeTree->getTriangleLocalNodeNum(iedge);
            break;
        case(ElementType::Beam):case(ElementType::Beam2):
            vertNum= pEdgeTree->getBeamLocalNodeNum(iedge);
            break;
        default:
            //Errorメッセージ出力
            pLogger= Utility::CLogger::Instance();
            pLogger->Info(Utility::LoggerMode::Error,"CommMesh2 Face, ShapeType Error @CommFace::getEdgePairCommNode");
            break;
    }

    PairCommNode pairCommNode;

    pairCommNode.first = mvCommNode[vertNum[0]];
    pairCommNode.second= mvCommNode[vertNum[1]];

    return pairCommNode;
}


// ----
// Edge集合処理のFaceへの処理
// ----
// 辺に隣接するFaceのセット タイプ1
void CCommFace::setEdgeCommFace(CCommFace* pNeibFace, const uiint& iedge)
{
    mvEdgeCommFace[iedge]= pNeibFace;
}
// 辺に隣接するFaceのセット タイプ2
void CCommFace::setEdgeCommFace(CCommFace* pNeibFace, PairCommNode& pairCommNode)
{
////    //debug
////    cout << "CommFace::setEdgeCommFace(pFace, pairCommNode),  pair.first=  " << pairCommNode.first->getID()
////                                                        << ",  pair.second= " << pairCommNode.second->getID() << endl;

    uiint iedge= getEdgeIndex(pairCommNode);
    
    mvEdgeCommFace[iedge]= pNeibFace;
}

// 辺ノードのセット 1
void CCommFace::setEdgeCommNode(CCommNode* pEdgeCommNode, const uiint& iedge)
{
    mvEdgeCommNode[iedge]= pEdgeCommNode;

////    //debug
////    if(pEdgeCommNode->getID()==9){
////        cout << "EdgeCommNode ID==9, iedge=" << iedge << ", CommFace ID= " << mID << ",  mgLevel= " << mMGLevel << endl;
////    }
////    cout << "CommFace::setEdgeCommNode, EdgeCommNode ID= " << pEdgeCommNode->getID()
////         << ", iedge= " << iedge
////         << ", level= " << mMGLevel
////         << ", faceID= " << mID << endl;
}
// 辺ノードのセット 2
void CCommFace::setEdgeCommNode(CCommNode* pEdgeCommNode, PairCommNode& pairCommNode)
{
    uiint iedge= getEdgeIndex(pairCommNode);
    
    mvEdgeCommNode[iedge]= pEdgeCommNode;
}

// 辺ノード設置済みマーキング 1
void CCommFace::markingEdgeNode(const uiint& iedge)
{
    mvbEdgeMarking[iedge]= true;;
}
// 辺ノード設置済みマーキング 2
void CCommFace::markingEdgeNode(PairCommNode& pairCommNode)
{
    uiint iedge= getEdgeIndex(pairCommNode);

    mvbEdgeMarking[iedge]= true;
}


// ----
// Refine
// ----
vector<CCommFace*>& CCommFace::refine(CElement *pElement)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();

    if(!pElement->isCommMesh2()){
        pLogger->Info(Utility::LoggerMode::Error, "invalid Element, at CCommFace::refine");

        cout << "CommFace ID = " << mID << ",  elemID= " << mElementID << ",  pElementID= " << pElement->getID() << endl;
    }

    CElement *pProgElem;//再分割された子供要素
    CNode    *pNode;//CommNodeが所有する,頂点のノード
    uiint nodeID;
    uiint progElemID;
    uiint progEntityID;
    
    uiint numOfVert;
    uiint numOfEdge;
    uiint nOrder;// 1次 || 2次 要素
    CCommFace* pCommFace;

    switch(mFaceType){
        case(ElementType::Quad):case(ElementType::Quad2):
            numOfVert= 4;
            numOfEdge= 4;
            nOrder = pElement->getOrder();

            
            // 頂点"0"にぶら下がる新CommFace
            pCommFace= new CCommFace;//// <<<<<<<<<<<< 新CommFace生成
            pCommFace->initialize(numOfVert, numOfEdge, nOrder);
            pCommFace->setCommNode(0, mvCommNode[0]);
            pCommFace->setCommNode(1, mvEdgeCommNode[0]);
            pCommFace->setCommNode(2, mpFaceCommNode);
            pCommFace->setCommNode(3, mvEdgeCommNode[3]);
            pCommFace->setMGLevel(mMGLevel+1);


            // CommNode[0]に対応するNodeにぶら下がるprogElement
            pNode= mvCommNode[0]->getNode();
            nodeID= pNode->getID();
            pProgElem= pElement->getProgElem_NodeID(nodeID);/////<<<<<<<<<< 頂点0に生成された子要素
            
            progElemID= pProgElem->getID();
            progEntityID= pProgElem->getCommEntityID();

            pCommFace->setElementID(progElemID);
            pCommFace->setElementFaceID(progEntityID);
            mvProgCommFace.push_back(pCommFace);



            // 頂点"1"にぶら下がる新CommFace
            pCommFace= new CCommFace;//// <<<<<<<<<<<< 新CommFace生成
            pCommFace->initialize(numOfVert, numOfEdge, nOrder);
            pCommFace->setCommNode(0, mvEdgeCommNode[0]);
            pCommFace->setCommNode(1, mvCommNode[1]);
            pCommFace->setCommNode(2, mvEdgeCommNode[1]);
            pCommFace->setCommNode(3, mpFaceCommNode);
            pCommFace->setMGLevel(mMGLevel+1);

            // CommNode[1]に対応するNodeにぶら下がるprogElement
            pNode= mvCommNode[1]->getNode();
            nodeID= pNode->getID();
            pProgElem= pElement->getProgElem_NodeID(nodeID);/////<<<<<<<<<< 頂点1に生成された子要素

            progElemID= pProgElem->getID();
            progEntityID= pProgElem->getCommEntityID();

            pCommFace->setElementID(progElemID);
            pCommFace->setElementFaceID(progEntityID);
            mvProgCommFace.push_back(pCommFace);


            // 頂点"2"にぶら下がる新CommFace
            pCommFace= new CCommFace;//// <<<<<<<<<<<< 新CommFace生成
            pCommFace->initialize(numOfVert, numOfEdge, nOrder);
            pCommFace->setCommNode(0, mpFaceCommNode);
            pCommFace->setCommNode(1, mvEdgeCommNode[1]);
            pCommFace->setCommNode(2, mvCommNode[2]);
            pCommFace->setCommNode(3, mvEdgeCommNode[2]);
            pCommFace->setMGLevel(mMGLevel+1);

            // CommNode[2]に対応するNodeにぶら下がるprogElement
            pNode= mvCommNode[2]->getNode();
            nodeID= pNode->getID();
            pProgElem= pElement->getProgElem_NodeID(nodeID);/////<<<<<<<<<< 頂点2に生成された子要素

            progElemID= pProgElem->getID();
            progEntityID= pProgElem->getCommEntityID();

            pCommFace->setElementID(progElemID);
            pCommFace->setElementFaceID(progEntityID);
            mvProgCommFace.push_back(pCommFace);


            // 頂点"3"にぶら下がる新CommFace
            pCommFace= new CCommFace;//// <<<<<<<<<<<< 新CommFace生成
            pCommFace->initialize(numOfVert, numOfEdge, nOrder);
            pCommFace->setCommNode(0, mvEdgeCommNode[3]);
            pCommFace->setCommNode(1, mpFaceCommNode);
            pCommFace->setCommNode(2, mvEdgeCommNode[2]);
            pCommFace->setCommNode(3, mvCommNode[3]);
            pCommFace->setMGLevel(mMGLevel+1);

            // CommNode[3]に対応するNodeにぶら下がるElement
            pNode= mvCommNode[3]->getNode();
            nodeID= pNode->getID();
            pProgElem= pElement->getProgElem_NodeID(nodeID);/////<<<<<<<<<< 頂点3に生成された子要素

            progElemID= pProgElem->getID();
            progEntityID= pProgElem->getCommEntityID();

            pCommFace->setElementID(progElemID);
            pCommFace->setElementFaceID(progEntityID);
            mvProgCommFace.push_back(pCommFace);

            break;

        case(ElementType::Triangle):case(ElementType::Triangle2):
            //三角形 => 四辺形を三個生成
            numOfVert= 4;
            numOfEdge= 4;
            nOrder= pElement->getOrder();

            // 頂点"0"にぶら下がる新CommFace
            pCommFace= new CCommFace;//// <<<<<<<<<<<< 新CommFace生成
            pCommFace->initialize(numOfVert, numOfEdge, nOrder);
            pCommFace->setCommNode(0, mvCommNode[0]);
            pCommFace->setCommNode(1, mvEdgeCommNode[0]);
            pCommFace->setCommNode(2, mpFaceCommNode);
            pCommFace->setCommNode(3, mvEdgeCommNode[2]);
            pCommFace->setMGLevel(mMGLevel+1);

            // CommNode[0]に対応するNodeにぶら下がるprogElement
            pNode= mvCommNode[0]->getNode();
            nodeID= pNode->getID();
            pProgElem= pElement->getProgElem_NodeID(nodeID);/////<<<<<<<<<< 頂点0に生成された子要素

            progElemID= pProgElem->getID();
            progEntityID= pProgElem->getCommEntityID();

            pCommFace->setElementID(progElemID);
            pCommFace->setElementFaceID(progEntityID);
            mvProgCommFace.push_back(pCommFace);



            // 頂点"1"にぶら下がる新CommFace
            pCommFace= new CCommFace;//// <<<<<<<<<<<< 新CommFace生成
            pCommFace->initialize(numOfVert, numOfEdge, nOrder);
            pCommFace->setCommNode(0, mvEdgeCommNode[0]);
            pCommFace->setCommNode(1, mvCommNode[1]);
            pCommFace->setCommNode(2, mvEdgeCommNode[1]);
            pCommFace->setCommNode(3, mpFaceCommNode);
            pCommFace->setMGLevel(mMGLevel+1);

            // CommNode[1]に対応するNodeにぶら下がるprogElement
            pNode= mvCommNode[1]->getNode();
            nodeID= pNode->getID();
            pProgElem= pElement->getProgElem_NodeID(nodeID);/////<<<<<<<<<< 頂点1に生成された子要素

            progElemID= pProgElem->getID();
            progEntityID= pProgElem->getCommEntityID();

            pCommFace->setElementID(progElemID);
            pCommFace->setElementFaceID(progEntityID);
            mvProgCommFace.push_back(pCommFace);



            // 頂点"2"にぶら下がる新CommFace
            pCommFace= new CCommFace;//// <<<<<<<<<<<< 新CommFace生成
            pCommFace->initialize(numOfVert, numOfEdge, nOrder);
            pCommFace->setCommNode(0, mvEdgeCommNode[1]);
            pCommFace->setCommNode(1, mvCommNode[2]);
            pCommFace->setCommNode(2, mvEdgeCommNode[2]);
            pCommFace->setCommNode(3, mpFaceCommNode);
            pCommFace->setMGLevel(mMGLevel+1);

            // CommNode[2]に対応するNodeにぶら下がるprogElement
            pNode= mvCommNode[2]->getNode();
            nodeID= pNode->getID();
            pProgElem= pElement->getProgElem_NodeID(nodeID);/////<<<<<<<<<< 頂点2に生成された子要素

            progElemID= pProgElem->getID();
            progEntityID= pProgElem->getCommEntityID();

            pCommFace->setElementID(progElemID);
            pCommFace->setElementFaceID(progEntityID);
            mvProgCommFace.push_back(pCommFace);
            
            break;

        case(ElementType::Beam):case(ElementType::Beam2):
            numOfVert= 2;
            numOfEdge= 1;
            nOrder= pElement->getOrder();

            // 頂点"0"にぶら下がる新CommFace
            pCommFace= new CCommFace;//// <<<<<<<<<<<< 新CommFace生成
            pCommFace->initialize(numOfVert, numOfEdge, nOrder);
            pCommFace->setCommNode(0, mvCommNode[0]);
            pCommFace->setCommNode(1, mvEdgeCommNode[0]);
            pCommFace->setMGLevel(mMGLevel+1);

            // CommNode[0]に対応するNodeにぶら下がるprogElement
            pNode= mvCommNode[0]->getNode();
            nodeID= pNode->getID();
            pProgElem= pElement->getProgElem_NodeID(nodeID);/////<<<<<<<<<< 頂点0に生成された子要素

            progElemID= pProgElem->getID();
            progEntityID= pProgElem->getCommEntityID();
            pCommFace->setElementID(progElemID);
            pCommFace->setElementFaceID(progEntityID);
            mvProgCommFace.push_back(pCommFace);


            // 頂点"1"にぶら下がる新CommFace
            pCommFace= new CCommFace;//// <<<<<<<<<<<< 新CommFace生成
            pCommFace->initialize(numOfVert, numOfEdge, nOrder);
            pCommFace->setCommNode(0, mvEdgeCommNode[0]);
            pCommFace->setCommNode(1, mvCommNode[1]);
            pCommFace->setMGLevel(mMGLevel+1);

            // CommNode[1]に対応するNodeにぶら下がるprogElement
            pNode= mvCommNode[1]->getNode();
            nodeID= pNode->getID();
            pProgElem= pElement->getProgElem_NodeID(nodeID);/////<<<<<<<<<< 頂点1に生成された子要素

            progElemID= pProgElem->getID();
            progEntityID= pProgElem->getCommEntityID();
            pCommFace->setElementID(progElemID);
            pCommFace->setElementFaceID(progEntityID);
            mvProgCommFace.push_back(pCommFace);
            
            break;
            
        default:
            break;
    }

    return mvProgCommFace;
}

//
// 2次要素面のときの辺ノード => mvCommNodeへ移設
//
void CCommFace::replaceEdgeCommNode()
{
    uiint nNumOfVert;
    switch(mFaceType){
        case(ElementType::Quad):
            return;////////////// quit
        case(ElementType::Triangle):
            return;////////////// quit
        case(ElementType::Beam):
            return;////////////// quit

        case(ElementType::Quad2):
            nNumOfVert=4;
            break;

        case(ElementType::Triangle2):
            nNumOfVert=3;
            break;

        case(ElementType::Beam2):
            nNumOfVert=2;
            break;
    }

    uiint iedge;
    for(iedge=0; iedge < mNumOfEdge; iedge++){
        mvCommNode[nNumOfVert + iedge] = mvEdgeCommNode[iedge];
    };
}

//
// Refine後のvector解放
// 
void CCommFace::deleteProgData()
{
    vector<CCommNode*>().swap(mvEdgeCommNode);// 辺に存在するNode (Refine用途)
    
    vector<CCommFace*>().swap(mvEdgeCommFace);// 辺に隣接するFace
    
    
    delete []mvbEdgeMarking;// 辺ノード生成マーキング
}







