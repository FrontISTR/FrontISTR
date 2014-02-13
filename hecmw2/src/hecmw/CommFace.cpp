/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/CommFace.cpp
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include <vector>
#include "Element.h"
#include "CommFace.h"
#include "ElementType.h"
#include "HEC_MPI.h"
using namespace pmw;
CCommFace::CCommFace()
{
    ;
}
CCommFace::~CCommFace()
{
    ;
}
void CCommFace::initialize(const uiint& numOfVert, const uiint& numOfEdge, const uiint& nOrder)
{
    switch(numOfVert) {
    case(4):
        if(ElementOrder::First== nOrder) {
            mFaceType= ElementType::Quad;
            mvCommNode.resize(4);
        }
        if(ElementOrder::Second==nOrder) {
            mFaceType= ElementType::Quad2;
            mvCommNode.resize(8);
        }
        break;
    case(3):
        if(ElementOrder::First== nOrder) {
            mFaceType= ElementType::Triangle;
            mvCommNode.resize(3);
        }
        if(ElementOrder::Second==nOrder) {
            mFaceType= ElementType::Triangle2;
            mvCommNode.resize(6);
        }
        break;
    case(2):
        if(ElementOrder::First== nOrder) {
            mFaceType= ElementType::Beam;
            mvCommNode.resize(2);
        }
        if(ElementOrder::Second==nOrder) {
            mFaceType= ElementType::Beam2;
            mvCommNode.resize(3);
        }
        break;
    case(1):
        mFaceType= ElementType::Point;//2012.04.27
        mvCommNode.resize(1);
        break;
    default:
        break;
    }
    //辺・面; 点の場合:0
    mvEdgeCommNode.resize(numOfEdge);
    mvEdgeCommFace.resize(numOfEdge);
    mNumOfEdge= numOfEdge;

    mvbEdgeMarking = new bool[numOfEdge];
    uiint iedge;
    for(iedge=0; iedge< numOfEdge; iedge++) {
        mvbEdgeMarking[iedge]= false;
    };
}
uiint CCommFace::getOrder()
{
    switch(mFaceType) {
        //1次
    case(ElementType::Quad):
    case(ElementType::Triangle):
    case(ElementType::Beam):
        return ElementOrder::First;
        //2次
    case(ElementType::Quad2):
    case(ElementType::Triangle2):
    case(ElementType::Beam2):
        return ElementOrder::Second;
        //0次:Point
    case(ElementType::Point):
        return ElementOrder::Zero;
    }
}
uiint CCommFace::getNumOfVert()
{
    switch(mFaceType) {
    case(ElementType::Quad):
    case(ElementType::Quad2):
        return 4;
    case(ElementType::Triangle):
    case(ElementType::Triangle2):
        return 3;
    case(ElementType::Beam):
    case(ElementType::Beam2):
        return 2;
    case(ElementType::Point):
        return 1;
    }
}
uiint& CCommFace::getEdgeIndex(PairCommNode& pairCommNode)
{
    CEdgeTree* pEdgeTree= CEdgeTree::Instance();

    //--
    //点
    //--
    if(mFaceType==ElementType::Point) return pEdgeTree->getDummyEdgeIndex();//辺は無い

    //--
    //辺・面
    //--
    uiint first_id=  pairCommNode.first->getID();
    uiint second_id= pairCommNode.second->getID();
    uiint self_id;
    uiint localNum0, localNum1;
    uiint numOfVert;
    uiint inode;
    if(mNumOfEdge==3) numOfVert = 3;
    if(mNumOfEdge==4) numOfVert = 4;
    if(mNumOfEdge==1) numOfVert = 2;

    for(inode=0; inode< numOfVert; inode++) {
        self_id=  mvCommNode[inode]->getID();
        if(self_id==first_id)  localNum0= inode;
        if(self_id==second_id) localNum1= inode;
    };

    switch(mFaceType) {
    case(ElementType::Quad):
    case(ElementType::Quad2):
        return pEdgeTree->getQuadEdgeIndex(localNum0, localNum1);
    case(ElementType::Triangle):
    case(ElementType::Triangle2):
        return pEdgeTree->getTriangleEdgeIndex(localNum0, localNum1);
    case(ElementType::Beam):
    case(ElementType::Beam2):
        return pEdgeTree->getBeamEdgeIndex(localNum0, localNum1);
    default:
        return pEdgeTree->getDisagTypeEdgeIndex(localNum0, localNum1);
    }
}
PairCommNode CCommFace::getEdgePairCommNode(const uiint& iedge)
{
    Utility::CLogger *pLogger;
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    PairCommNode pairCommNode;
    //--
    //点
    //--
    if(mFaceType==ElementType::Point) {
        pairCommNode.first = mvCommNode[0];
        pairCommNode.second= mvCommNode[0];
        return pairCommNode;
    }

    //--
    //辺・面
    //--
    uiint* vertNum;
    switch(mFaceType) {
    case(ElementType::Quad):
    case(ElementType::Quad2):
        vertNum= pEdgeTree->getQuadLocalNodeNum(iedge);
        break;
    case(ElementType::Triangle):
    case(ElementType::Triangle2):
        vertNum= pEdgeTree->getTriangleLocalNodeNum(iedge);
        break;
    case(ElementType::Beam):
    case(ElementType::Beam2):
        vertNum= pEdgeTree->getBeamLocalNodeNum(iedge);
        break;
    default:
        pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"CommMesh2 Face, ShapeType Error @CommFace::getEdgePairCommNode");
        break;
    }

    pairCommNode.first = mvCommNode[vertNum[0]];
    pairCommNode.second= mvCommNode[vertNum[1]];
    return pairCommNode;
}
void CCommFace::setEdgeCommFace(CCommFace* pNeibFace, const uiint& iedge)
{
    //Point
    if(mFaceType==ElementType::Point) {
        //Error
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "CommFace::setEdgeCommFace(CCommFace* pNeibFace, const uiint& iedge), this is Point");
        return;
    }

    //Edge, Face
    mvEdgeCommFace[iedge]= pNeibFace;
}
void CCommFace::setEdgeCommFace(CCommFace* pNeibFace, PairCommNode& pairCommNode)
{
    //Point
    if(mFaceType==ElementType::Point) {
        //Error
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "CommFace::setEdgeCommFace(CCommFace* pNeibFace, PairCommNode& pairCommNode), this is Point");
        return;
    }
    //Edge, Face
    uiint iedge= getEdgeIndex(pairCommNode);
    mvEdgeCommFace[iedge]= pNeibFace;
}
void CCommFace::setEdgeCommNode(CCommNode* pEdgeCommNode, const uiint& iedge)
{
    //Point
    if(mFaceType==ElementType::Point) {
        //Error
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "CommFace::setEdgeCommNode(CCommNode* pEdgeCommNode, const uiint& iedge), this is Point");
        return;
    }
    //Edge, Face
    mvEdgeCommNode[iedge]= pEdgeCommNode;
}
void CCommFace::setEdgeCommNode(CCommNode* pEdgeCommNode, PairCommNode& pairCommNode)
{
    //Point
    if(mFaceType==ElementType::Point) {
        //Error
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "CommFace::setEdgeCommNode(CCommNode* pEdgeCommNode, PairCommNode& pairCommNode), this is Point");
        return;
    }
    //Edge, Face
    uiint iedge= getEdgeIndex(pairCommNode);
    mvEdgeCommNode[iedge]= pEdgeCommNode;
}
void CCommFace::markingEdgeNode(const uiint& iedge)
{
    //Point
    if(mFaceType==ElementType::Point) {
        //Error
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "CommFace::markingEdgeNode(const uiint& iedge), this is Point");
        return;
    }
    //Edge, Face
    mvbEdgeMarking[iedge]= true;;
}
void CCommFace::markingEdgeNode(PairCommNode& pairCommNode)
{
    //Point
    if(mFaceType==ElementType::Point) {
        //Error
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "CommFace::markingEdgeNode(PairCommNode& pairCommNode), this is Point");
        return;
    }
    //Edge, Face
    uiint iedge= getEdgeIndex(pairCommNode);
    mvbEdgeMarking[iedge]= true;
}
//--
// 自身をRefine
//--
vector<CCommFace*>& CCommFace::refine(CElement *pElement)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();

    CHecMPI *pMPI=CHecMPI::Instance();
    ////cout << "CommFace::refine -------------- entry   rank:" << pMPI->getRank() << endl;

    //----------------------
    // Edge,Face : Refine対象
    //----------------------
    if(!pElement->isCommMesh2()) {
        pLogger->Info(Utility::LoggerMode::Error, "invalid Element, at CCommFace::refine");
        string sTypeName;
        switch(mFaceType) {
        case(ElementType::Quad):
            sTypeName="Quad";
            break;
        case(ElementType::Quad2):
            sTypeName="Quad2";
            break;
        case(ElementType::Triangle):
            sTypeName="Triangle";
            break;
        case(ElementType::Triangle2):
            sTypeName="Triangle2";
            break;
        case(ElementType::Beam):
            sTypeName="Beam";
            break;
        case(ElementType::Beam2):
            sTypeName="Beam2";
            break;
        case(ElementType::Point):
            sTypeName="Point";
            break;
        }
        cout << "CommFace ID:" << mID << " Type:" << sTypeName
             << ",  elemID:" << mElementID
             << ",  pElementID:" << pElement->getID() << " rank:" << pMPI->getRank() << endl;
    }
    CElement *pProgElem;
    CNode    *pNode;
    uiint nodeID;
    uiint progElemID;
    uiint progEntityID;
    uiint numOfVert;
    uiint numOfEdge;
    uiint nOrder;
    CCommFace* pCommFace;
    switch(mFaceType) {
        //--
        // 四辺形
        //--
    case(ElementType::Quad):
    case(ElementType::Quad2):
        numOfVert= 4;
        numOfEdge= 4;
        nOrder = pElement->getOrder();

        pCommFace= new CCommFace;
        pCommFace->initialize(numOfVert, numOfEdge, nOrder);
        pCommFace->setCommNode(0, mvCommNode[0]);
        pCommFace->setCommNode(1, mvEdgeCommNode[0]);
        pCommFace->setCommNode(2, mpFaceCommNode);
        pCommFace->setCommNode(3, mvEdgeCommNode[3]);
        pCommFace->setMGLevel(mMGLevel+1);
        pNode= mvCommNode[0]->getNode();
        nodeID= pNode->getID();
        pProgElem= pElement->getProgElem_NodeID(nodeID);
        progElemID= pProgElem->getID();
        pCommFace->setElementID(progElemID);
        mvProgCommFace.push_back(pCommFace);

        pCommFace= new CCommFace;
        pCommFace->initialize(numOfVert, numOfEdge, nOrder);
        pCommFace->setCommNode(0, mvEdgeCommNode[0]);
        pCommFace->setCommNode(1, mvCommNode[1]);
        pCommFace->setCommNode(2, mvEdgeCommNode[1]);
        pCommFace->setCommNode(3, mpFaceCommNode);
        pCommFace->setMGLevel(mMGLevel+1);
        pNode= mvCommNode[1]->getNode();
        nodeID= pNode->getID();
        pProgElem= pElement->getProgElem_NodeID(nodeID);
        progElemID= pProgElem->getID();
        pCommFace->setElementID(progElemID);
        mvProgCommFace.push_back(pCommFace);

        pCommFace= new CCommFace;
        pCommFace->initialize(numOfVert, numOfEdge, nOrder);
        pCommFace->setCommNode(0, mpFaceCommNode);
        pCommFace->setCommNode(1, mvEdgeCommNode[1]);
        pCommFace->setCommNode(2, mvCommNode[2]);
        pCommFace->setCommNode(3, mvEdgeCommNode[2]);
        pCommFace->setMGLevel(mMGLevel+1);
        pNode= mvCommNode[2]->getNode();
        nodeID= pNode->getID();
        pProgElem= pElement->getProgElem_NodeID(nodeID);
        progElemID= pProgElem->getID();
        pCommFace->setElementID(progElemID);
        mvProgCommFace.push_back(pCommFace);

        pCommFace= new CCommFace;
        pCommFace->initialize(numOfVert, numOfEdge, nOrder);
        pCommFace->setCommNode(0, mvEdgeCommNode[3]);
        pCommFace->setCommNode(1, mpFaceCommNode);
        pCommFace->setCommNode(2, mvEdgeCommNode[2]);
        pCommFace->setCommNode(3, mvCommNode[3]);
        pCommFace->setMGLevel(mMGLevel+1);
        pNode= mvCommNode[3]->getNode();
        nodeID= pNode->getID();
        pProgElem= pElement->getProgElem_NodeID(nodeID);
        progElemID= pProgElem->getID();
        pCommFace->setElementID(progElemID);
        mvProgCommFace.push_back(pCommFace);
        break;
        //--
        // 三角形 => refine後:四辺形
        //--
    case(ElementType::Triangle):
    case(ElementType::Triangle2):
        //--
        // EdgeCommNodeには,まだNodeはセットされていない.
        //--
        numOfVert= 4;
        numOfEdge= 4;
        nOrder= pElement->getOrder();

        pCommFace= new CCommFace;//-------------------------------生成
        pCommFace->initialize(numOfVert, numOfEdge, nOrder);
        pCommFace->setCommNode(0, mvCommNode[0]);
        pCommFace->setCommNode(1, mvEdgeCommNode[0]);
        pCommFace->setCommNode(2, mpFaceCommNode);
        pCommFace->setCommNode(3, mvEdgeCommNode[2]);
        pCommFace->setMGLevel(mMGLevel+1);
        pNode= mvCommNode[0]->getNode();
        nodeID= pNode->getID();
        pProgElem= pElement->getProgElem_NodeID(nodeID);
        progElemID= pProgElem->getID();
        pCommFace->setElementID(progElemID);
        mvProgCommFace.push_back(pCommFace);

        pCommFace= new CCommFace;//-------------------------------生成
        pCommFace->initialize(numOfVert, numOfEdge, nOrder);
        pCommFace->setCommNode(0, mvEdgeCommNode[0]);
        pCommFace->setCommNode(1, mvCommNode[1]);
        pCommFace->setCommNode(2, mvEdgeCommNode[1]);
        pCommFace->setCommNode(3, mpFaceCommNode);
        pCommFace->setMGLevel(mMGLevel+1);
        pNode= mvCommNode[1]->getNode();
        nodeID= pNode->getID();
        pProgElem= pElement->getProgElem_NodeID(nodeID);
        progElemID= pProgElem->getID();
        pCommFace->setElementID(progElemID);
        mvProgCommFace.push_back(pCommFace);

        pCommFace= new CCommFace;//-------------------------------生成
        pCommFace->initialize(numOfVert, numOfEdge, nOrder);
        pCommFace->setCommNode(0, mvEdgeCommNode[1]);
        pCommFace->setCommNode(1, mvCommNode[2]);
        pCommFace->setCommNode(2, mvEdgeCommNode[2]);
        pCommFace->setCommNode(3, mpFaceCommNode);
        pCommFace->setMGLevel(mMGLevel+1);
        pNode= mvCommNode[2]->getNode();
        nodeID= pNode->getID();
        pProgElem= pElement->getProgElem_NodeID(nodeID);
        progElemID= pProgElem->getID();
        pCommFace->setElementID(progElemID);
        mvProgCommFace.push_back(pCommFace);
        break;
        //--
        // 辺
        //--
    case(ElementType::Beam):
    case(ElementType::Beam2):
        numOfVert= 2;
        numOfEdge= 1;
        nOrder= pElement->getOrder();

        pCommFace= new CCommFace;//-----------------------------生成
        pCommFace->initialize(numOfVert, numOfEdge, nOrder);

        pCommFace->setCommNode(0, mvCommNode[0]);
        pCommFace->setCommNode(1, mvEdgeCommNode[0]);

        pCommFace->setMGLevel(mMGLevel+1);
        pNode= mvCommNode[0]->getNode();
        nodeID= pNode->getID();
        pProgElem= pElement->getProgElem_NodeID(nodeID);
        progElemID= pProgElem->getID();
        pCommFace->setElementID(progElemID);
        mvProgCommFace.push_back(pCommFace);

        pCommFace= new CCommFace;//----------------------------生成
        pCommFace->initialize(numOfVert, numOfEdge, nOrder);

        pCommFace->setCommNode(0, mvEdgeCommNode[0]);
        pCommFace->setCommNode(1, mvCommNode[1]);

        pCommFace->setMGLevel(mMGLevel+1);
        pNode= mvCommNode[1]->getNode();
        nodeID= pNode->getID();
        pProgElem= pElement->getProgElem_NodeID(nodeID);
        progElemID= pProgElem->getID();
        pCommFace->setElementID(progElemID);
        mvProgCommFace.push_back(pCommFace);
        break;
        //--
        // 点 '12.04.27
        //--
    case(ElementType::Point):
        numOfVert= 1;
        numOfEdge= 0;
        nOrder= ElementOrder::Zero;

        pCommFace= new CCommFace;//----------- Level値,要素が異なるのでLevel毎に新規生成
        pCommFace->initialize(numOfVert, numOfEdge, nOrder);
        pCommFace->setCommNode(0, mvCommNode[0]);
        pCommFace->setMGLevel(mMGLevel+1);

        pNode= mvCommNode[0]->getNode();
        nodeID= pNode->getID();
        pProgElem= pElement->getProgElem_NodeID(nodeID);//--- 頂点にぶら下がるElement(Refine先)
        progElemID= pProgElem->getID();
        pCommFace->setElementID(progElemID);
        mvProgCommFace.push_back(pCommFace);
        break;
    default:
        break;
    }

    ////cout << "CommFace::refine -------------- exit   rank:" << pMPI->getRank() << endl;

    return mvProgCommFace;
}
//
// Factoryの面中心ノード取得で使用
//
uiint& CCommFace::getElementFaceID(CElement *pElem)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();

    if(mFaceType==ElementType::Beam ||  mFaceType==ElementType::Beam2 || mFaceType==ElementType::Point) {
        //error message
        pLogger->Info(Utility::LoggerMode::Error, " CommFace::getElementFaceID ---- invalid FaceType");

        return pLogger->getUDummyValue();
    }

    CNode *pNode0= mvCommNode[0]->getNode();
    CNode *pNode1= mvCommNode[1]->getNode();
    CNode *pNode2= mvCommNode[2]->getNode();

    return pElem->getFaceIndex(pNode0, pNode1, pNode2);
}



void CCommFace::replaceEdgeCommNode()
{
    //--
    // 点 :何もしない
    //--
    if(mFaceType==ElementType::Point) return;

    //--
    // 辺・面
    //--
    uiint nNumOfVert;
    switch(mFaceType) {
    case(ElementType::Quad):
        return;
    case(ElementType::Triangle):
        return;
    case(ElementType::Beam):
        return;
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
    for(iedge=0; iedge < mNumOfEdge; iedge++) {
        mvCommNode[nNumOfVert + iedge] = mvEdgeCommNode[iedge];
    };
}
void CCommFace::deleteProgData()
{
    vector<CCommNode*>().swap(mvEdgeCommNode);
    vector<CCommFace*>().swap(mvEdgeCommFace);
    delete []mvbEdgeMarking;
}
