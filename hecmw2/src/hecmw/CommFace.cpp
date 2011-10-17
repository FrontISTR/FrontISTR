/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/CommFace.cpp
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
void CCommFace::initialize(const uiint& numOfVert, const uiint& numOfEdge, const uiint& nOrder)
{
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
    mvEdgeCommNode.resize(numOfEdge);
    mvEdgeCommFace.resize(numOfEdge);
    mNumOfEdge= numOfEdge;
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
    uiint first_id=  pairCommNode.first->getID();
    uiint second_id= pairCommNode.second->getID();
    uiint self_id;
    uiint localNum0, localNum1;
    uiint numOfVert;
    uiint inode;
    if(mNumOfEdge==3) numOfVert = 3;
    if(mNumOfEdge==4) numOfVert = 4;
    if(mNumOfEdge==1) numOfVert = 2;
    for(inode=0; inode< numOfVert; inode++){
        self_id=  mvCommNode[inode]->getID();
        if(self_id==first_id)  localNum0= inode;
        if(self_id==second_id) localNum1= inode;
    };
    CEdgeTree* pEdgeTree= CEdgeTree::Instance();
    switch(mFaceType){
        case(ElementType::Quad):case(ElementType::Quad2):
            return pEdgeTree->getQuadEdgeIndex(localNum0, localNum1);
        case(ElementType::Triangle):case(ElementType::Triangle2):
            return pEdgeTree->getTriangleEdgeIndex(localNum0, localNum1);
        case(ElementType::Beam):case(ElementType::Beam2):
            return pEdgeTree->getBeamEdgeIndex(localNum0, localNum1);
        default:
            return pEdgeTree->getDisagTypeEdgeIndex(localNum0, localNum1);
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
            pLogger= Utility::CLogger::Instance();
            pLogger->Info(Utility::LoggerMode::Error,"CommMesh2 Face, ShapeType Error @CommFace::getEdgePairCommNode");
            break;
    }
    PairCommNode pairCommNode;
    pairCommNode.first = mvCommNode[vertNum[0]];
    pairCommNode.second= mvCommNode[vertNum[1]];
    return pairCommNode;
}
void CCommFace::setEdgeCommFace(CCommFace* pNeibFace, const uiint& iedge)
{
    mvEdgeCommFace[iedge]= pNeibFace;
}
void CCommFace::setEdgeCommFace(CCommFace* pNeibFace, PairCommNode& pairCommNode)
{
    uiint iedge= getEdgeIndex(pairCommNode);
    mvEdgeCommFace[iedge]= pNeibFace;
}
void CCommFace::setEdgeCommNode(CCommNode* pEdgeCommNode, const uiint& iedge)
{
    mvEdgeCommNode[iedge]= pEdgeCommNode;
}
void CCommFace::setEdgeCommNode(CCommNode* pEdgeCommNode, PairCommNode& pairCommNode)
{
    uiint iedge= getEdgeIndex(pairCommNode);
    mvEdgeCommNode[iedge]= pEdgeCommNode;
}
void CCommFace::markingEdgeNode(const uiint& iedge)
{
    mvbEdgeMarking[iedge]= true;;
}
void CCommFace::markingEdgeNode(PairCommNode& pairCommNode)
{
    uiint iedge= getEdgeIndex(pairCommNode);
    mvbEdgeMarking[iedge]= true;
}
vector<CCommFace*>& CCommFace::refine(CElement *pElement)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    if(!pElement->isCommMesh2()){
        pLogger->Info(Utility::LoggerMode::Error, "invalid Element, at CCommFace::refine");
        cout << "CommFace ID = " << mID << ",  elemID= " << mElementID << ",  pElementID= " << pElement->getID() << endl;
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
    switch(mFaceType){
        case(ElementType::Quad):case(ElementType::Quad2):
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
            progEntityID= pProgElem->getCommEntityID();
            pCommFace->setElementID(progElemID);
            pCommFace->setElementFaceID(progEntityID);
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
            progEntityID= pProgElem->getCommEntityID();
            pCommFace->setElementID(progElemID);
            pCommFace->setElementFaceID(progEntityID);
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
            progEntityID= pProgElem->getCommEntityID();
            pCommFace->setElementID(progElemID);
            pCommFace->setElementFaceID(progEntityID);
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
            progEntityID= pProgElem->getCommEntityID();
            pCommFace->setElementID(progElemID);
            pCommFace->setElementFaceID(progEntityID);
            mvProgCommFace.push_back(pCommFace);
            break;
        case(ElementType::Triangle):case(ElementType::Triangle2):
            numOfVert= 4;
            numOfEdge= 4;
            nOrder= pElement->getOrder();
            pCommFace= new CCommFace;
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
            progEntityID= pProgElem->getCommEntityID();
            pCommFace->setElementID(progElemID);
            pCommFace->setElementFaceID(progEntityID);
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
            progEntityID= pProgElem->getCommEntityID();
            pCommFace->setElementID(progElemID);
            pCommFace->setElementFaceID(progEntityID);
            mvProgCommFace.push_back(pCommFace);
            pCommFace= new CCommFace;
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
            progEntityID= pProgElem->getCommEntityID();
            pCommFace->setElementID(progElemID);
            pCommFace->setElementFaceID(progEntityID);
            mvProgCommFace.push_back(pCommFace);
            break;
        case(ElementType::Beam):case(ElementType::Beam2):
            numOfVert= 2;
            numOfEdge= 1;
            nOrder= pElement->getOrder();
            pCommFace= new CCommFace;
            pCommFace->initialize(numOfVert, numOfEdge, nOrder);
            pCommFace->setCommNode(0, mvCommNode[0]);
            pCommFace->setCommNode(1, mvEdgeCommNode[0]);
            pCommFace->setMGLevel(mMGLevel+1);
            pNode= mvCommNode[0]->getNode();
            nodeID= pNode->getID();
            pProgElem= pElement->getProgElem_NodeID(nodeID);
            progElemID= pProgElem->getID();
            progEntityID= pProgElem->getCommEntityID();
            pCommFace->setElementID(progElemID);
            pCommFace->setElementFaceID(progEntityID);
            mvProgCommFace.push_back(pCommFace);
            pCommFace= new CCommFace;
            pCommFace->initialize(numOfVert, numOfEdge, nOrder);
            pCommFace->setCommNode(0, mvEdgeCommNode[0]);
            pCommFace->setCommNode(1, mvCommNode[1]);
            pCommFace->setMGLevel(mMGLevel+1);
            pNode= mvCommNode[1]->getNode();
            nodeID= pNode->getID();
            pProgElem= pElement->getProgElem_NodeID(nodeID);
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
void CCommFace::replaceEdgeCommNode()
{
    uiint nNumOfVert;
    switch(mFaceType){
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
    for(iedge=0; iedge < mNumOfEdge; iedge++){
        mvCommNode[nNumOfVert + iedge] = mvEdgeCommNode[iedge];
    };
}
void CCommFace::deleteProgData()
{
    vector<CCommNode*>().swap(mvEdgeCommNode);
    vector<CCommFace*>().swap(mvEdgeCommFace);
    delete []mvbEdgeMarking;
}
