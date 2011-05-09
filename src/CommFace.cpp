/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   CommFace.cxx
|
|                     Written by T.Takeda,    2010/06/01
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
using namespace pmw;
CCommFace::CCommFace()
{
    ;
}
CCommFace::~CCommFace()
{
    ;
}
void CCommFace::initialize(const uint& numOfVert, const uint& numOfEdge)
{
    mvCommNode.resize(numOfVert);    
    mvEdgeCommNode.resize(numOfEdge);
    mvEdgeCommFace.resize(numOfEdge);
    mNumOfEdge= numOfEdge;
    mvbEdgeMarking.resize(numOfEdge);
    uint iedge;
    for(iedge=0; iedge< numOfEdge; iedge++){
        mvbEdgeMarking[iedge]= false;
    };
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
    uint localNum0, localNum1;
    uint numOfNode= mvCommNode.size();
    uint inode;
    for(inode=0; inode< numOfNode; inode++){
        if(mvCommNode[inode]->getID() == pairCommNode.first->getID() ) localNum0=inode;
        if(mvCommNode[inode]->getID() == pairCommNode.second->getID()) localNum1=inode;
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
PairCommNode& CCommFace::getEdgePairCommNode(const uint& iedge)
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
            pLogger= Utility::CLogger::Instance();
            pLogger->Info(Utility::LoggerMode::Error,"CommMesh2 Face, ShapeType Error @CommFace::getEdgePairCommNode");
            break;
    }
    mPairCommNode.first = mvCommNode[vertNum[0]];
    mPairCommNode.second= mvCommNode[vertNum[1]];
    return mPairCommNode;
}
void CCommFace::setEdgeCommFace(CCommFace* pNeibFace, const uint& iedge)
{
    mvEdgeCommFace[iedge]= pNeibFace;
}
void CCommFace::setEdgeCommFace(CCommFace* pNeibFace, PairCommNode& pairCommNode)
{
    uint iedge= getEdgeIndex(pairCommNode);
    mvEdgeCommFace[iedge]= pNeibFace;
}
void CCommFace::setEdgeCommNode(CCommNode* pEdgeCommNode, const uint& iedge)
{
    mvEdgeCommNode[iedge]= pEdgeCommNode;
}
void CCommFace::setEdgeCommNode(CCommNode* pEdgeCommNode, PairCommNode& pairCommNode)
{
    uint iedge= getEdgeIndex(pairCommNode);
    mvEdgeCommNode[iedge]= pEdgeCommNode;
}
void CCommFace::markingEdgeNode(const uint& iedge)
{
    mvbEdgeMarking[iedge]= true;;
}
void CCommFace::markingEdgeNode(PairCommNode& pairCommNode)
{
    uint iedge= getEdgeIndex(pairCommNode);
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
    uint nodeID;
    uint progElemID;
    uint progEntityID;
    uint numOfVert= mvCommNode.size();
    uint numOfEdge;
    CCommFace* pCommFace;
    switch(numOfVert){
        case(4):
            numOfEdge= 4;
            pCommFace= new CCommFace;
            pCommFace->initialize(numOfVert, numOfEdge);
            pCommFace->setVertCommNode(0, mvCommNode[0]);
            pCommFace->setVertCommNode(1, mvEdgeCommNode[0]);
            pCommFace->setVertCommNode(2, mpFaceCommNode);
            pCommFace->setVertCommNode(3, mvEdgeCommNode[3]);
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
            pCommFace->initialize(numOfVert, numOfEdge);
            pCommFace->setVertCommNode(0, mvEdgeCommNode[0]);
            pCommFace->setVertCommNode(1, mvCommNode[1]);
            pCommFace->setVertCommNode(2, mvEdgeCommNode[1]);
            pCommFace->setVertCommNode(3, mpFaceCommNode);
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
            pCommFace->initialize(numOfVert, numOfEdge);
            pCommFace->setVertCommNode(0, mpFaceCommNode);
            pCommFace->setVertCommNode(1, mvEdgeCommNode[1]);
            pCommFace->setVertCommNode(2, mvCommNode[2]);
            pCommFace->setVertCommNode(3, mvEdgeCommNode[2]);
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
            pCommFace->initialize(numOfVert, numOfEdge);
            pCommFace->setVertCommNode(0, mvEdgeCommNode[3]);
            pCommFace->setVertCommNode(1, mpFaceCommNode);
            pCommFace->setVertCommNode(2, mvEdgeCommNode[2]);
            pCommFace->setVertCommNode(3, mvCommNode[3]);
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
        case(3):
            numOfEdge= 3;
            pCommFace= new CCommFace;
            pCommFace->initialize(numOfVert, numOfEdge);
            pCommFace->setVertCommNode(0, mvCommNode[0]);
            pCommFace->setVertCommNode(1, mvEdgeCommNode[0]);
            pCommFace->setVertCommNode(2, mpFaceCommNode);
            pCommFace->setVertCommNode(3, mvEdgeCommNode[2]);
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
            pCommFace->initialize(numOfVert, numOfEdge);
            pCommFace->setVertCommNode(0, mvEdgeCommNode[0]);
            pCommFace->setVertCommNode(1, mvCommNode[1]);
            pCommFace->setVertCommNode(2, mvEdgeCommNode[1]);
            pCommFace->setVertCommNode(3, mpFaceCommNode);
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
            pCommFace->initialize(numOfVert, numOfEdge);
            pCommFace->setVertCommNode(0, mvEdgeCommNode[1]);
            pCommFace->setVertCommNode(1, mvCommNode[2]);
            pCommFace->setVertCommNode(2, mvEdgeCommNode[2]);
            pCommFace->setVertCommNode(3, mpFaceCommNode);
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
        case(2):
            numOfEdge= 1;
            pCommFace= new CCommFace;
            pCommFace->initialize(numOfVert, numOfEdge);
            pCommFace->setVertCommNode(0, mvCommNode[0]);
            pCommFace->setVertCommNode(1, mvEdgeCommNode[0]);
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
            pCommFace->initialize(numOfVert, numOfEdge);
            pCommFace->setVertCommNode(0, mvEdgeCommNode[0]);
            pCommFace->setVertCommNode(1, mvCommNode[1]);
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
