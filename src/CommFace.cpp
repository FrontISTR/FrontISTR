
#include <vector>

#include "Element.h"

//
//  CommFace.cpp
//
//
//
//              2010.03.02
//              k.Takeda
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
void CCommFace::initialize(const uint& numOfVert, const uint& numOfEdge)
{
    mvCommNode.resize(numOfVert);    //頂点ノード
    mvEdgeCommNode.resize(numOfEdge);//辺ノード
    
    mvEdgeCommFace.resize(numOfEdge);//辺に隣接する”面”

    mNumOfEdge= numOfEdge;

    //辺マーキングの初期化
    mvbEdgeMarking = new bool[numOfEdge];
    uint iedge;
    for(iedge=0; iedge< numOfEdge; iedge++){
        mvbEdgeMarking[iedge]= false;
    };


    //FaceTypeの設定
    switch(numOfVert){
        case(4):
            mFaceType= ElementType::Quad;
            break;
        case(3):
            mFaceType= ElementType::Triangle;
            break;
        case(2):
            mFaceType= ElementType::Beam;
            break;
        default:
            break;
    }
}


uint& CCommFace::getEdgeIndex(PairCommNode& pairCommNode)
{
    // Node0, Node1のインデックス番号から辺番号を特定
    uint first_id=  pairCommNode.first->getID();
    uint second_id= pairCommNode.second->getID();
    uint self_id;
    uint localNum0, localNum1;
    uint numOfNode= mvCommNode.size();
    uint inode;

    ////debug
    //cout << "CommFace::getEdgeIndex(pairCommNode), numOfNode= " << numOfNode << endl;

    for(inode=0; inode< numOfNode; inode++){
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

PairCommNode CCommFace::getEdgePairCommNode(const uint& iedge)
{
    Utility::CLogger *pLogger;

    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    uint* vertNum;
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
void CCommFace::setEdgeCommFace(CCommFace* pNeibFace, const uint& iedge)
{
    mvEdgeCommFace[iedge]= pNeibFace;
}
// 辺に隣接するFaceのセット タイプ2
void CCommFace::setEdgeCommFace(CCommFace* pNeibFace, PairCommNode& pairCommNode)
{
//    //debug
//    cout << "CCommFace::setEdgeCommFace(pFace, pairCommNode),  pair.first=  " << pairCommNode.first->getID()
//                                                        << ",  pair.second= " << pairCommNode.second->getID() << endl;

    uint iedge= getEdgeIndex(pairCommNode);
    
    mvEdgeCommFace[iedge]= pNeibFace;
}

// 辺ノードのセット 1
void CCommFace::setEdgeCommNode(CCommNode* pEdgeCommNode, const uint& iedge)
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
    uint iedge= getEdgeIndex(pairCommNode);
    
    mvEdgeCommNode[iedge]= pEdgeCommNode;
}

// 辺ノード設置済みマーキング 1
void CCommFace::markingEdgeNode(const uint& iedge)
{
    mvbEdgeMarking[iedge]= true;;
}
// 辺ノード設置済みマーキング 2
void CCommFace::markingEdgeNode(PairCommNode& pairCommNode)
{
    uint iedge= getEdgeIndex(pairCommNode);

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
    uint nodeID;
    uint progElemID;
    uint progEntityID;
    
    uint numOfVert= mvCommNode.size();
    uint numOfEdge;
    CCommFace* pCommFace;
    switch(numOfVert){
        case(4):
            numOfEdge= 4;

            // 頂点"0"にぶら下がる新CommFace
            pCommFace= new CCommFace;//// <<<<<<<<<<<< 新CommFace生成
            pCommFace->initialize(numOfVert, numOfEdge);
            pCommFace->setVertCommNode(0, mvCommNode[0]);
            pCommFace->setVertCommNode(1, mvEdgeCommNode[0]);
            pCommFace->setVertCommNode(2, mpFaceCommNode);
            pCommFace->setVertCommNode(3, mvEdgeCommNode[3]);
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
            pCommFace->initialize(numOfVert, numOfEdge);
            pCommFace->setVertCommNode(0, mvEdgeCommNode[0]);
            pCommFace->setVertCommNode(1, mvCommNode[1]);
            pCommFace->setVertCommNode(2, mvEdgeCommNode[1]);
            pCommFace->setVertCommNode(3, mpFaceCommNode);
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
            pCommFace->initialize(numOfVert, numOfEdge);
            pCommFace->setVertCommNode(0, mpFaceCommNode);
            pCommFace->setVertCommNode(1, mvEdgeCommNode[1]);
            pCommFace->setVertCommNode(2, mvCommNode[2]);
            pCommFace->setVertCommNode(3, mvEdgeCommNode[2]);
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
            pCommFace->initialize(numOfVert, numOfEdge);
            pCommFace->setVertCommNode(0, mvEdgeCommNode[3]);
            pCommFace->setVertCommNode(1, mpFaceCommNode);
            pCommFace->setVertCommNode(2, mvEdgeCommNode[2]);
            pCommFace->setVertCommNode(3, mvCommNode[3]);
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
        case(3):
            numOfEdge= 3;

            // 頂点"0"にぶら下がる新CommFace
            pCommFace= new CCommFace;//// <<<<<<<<<<<< 新CommFace生成
            pCommFace->initialize(numOfVert, numOfEdge);
            pCommFace->setVertCommNode(0, mvCommNode[0]);
            pCommFace->setVertCommNode(1, mvEdgeCommNode[0]);
            pCommFace->setVertCommNode(2, mpFaceCommNode);
            pCommFace->setVertCommNode(3, mvEdgeCommNode[2]);
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
            pCommFace->initialize(numOfVert, numOfEdge);
            pCommFace->setVertCommNode(0, mvEdgeCommNode[0]);
            pCommFace->setVertCommNode(1, mvCommNode[1]);
            pCommFace->setVertCommNode(2, mvEdgeCommNode[1]);
            pCommFace->setVertCommNode(3, mpFaceCommNode);
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
            pCommFace->initialize(numOfVert, numOfEdge);
            pCommFace->setVertCommNode(0, mvEdgeCommNode[1]);
            pCommFace->setVertCommNode(1, mvCommNode[2]);
            pCommFace->setVertCommNode(2, mvEdgeCommNode[2]);
            pCommFace->setVertCommNode(3, mpFaceCommNode);
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
        case(2):
            numOfEdge= 1;

            // 頂点"0"にぶら下がる新CommFace
            pCommFace= new CCommFace;//// <<<<<<<<<<<< 新CommFace生成
            pCommFace->initialize(numOfVert, numOfEdge);
            pCommFace->setVertCommNode(0, mvCommNode[0]);
            pCommFace->setVertCommNode(1, mvEdgeCommNode[0]);
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
            pCommFace->initialize(numOfVert, numOfEdge);
            pCommFace->setVertCommNode(0, mvEdgeCommNode[0]);
            pCommFace->setVertCommNode(1, mvCommNode[1]);
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



// Refine後のvector解放
// ----
void CCommFace::deleteProgData()
{
    vector<CCommNode*>().swap(mvEdgeCommNode);// 辺に存在するNode (Refine用途)
    
    vector<CCommFace*>().swap(mvEdgeCommFace);// 辺に隣接するFace
    
    
    delete []mvbEdgeMarking;// 辺ノード生成マーキング
}







