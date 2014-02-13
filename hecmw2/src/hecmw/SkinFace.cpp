/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/SkinFace.cpp
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
#include "ElementType.h"
#include "Element.h"
#include <vector>
#include "SkinFace.h"
using namespace pmw;
CSkinFace::CSkinFace()
{
    mbSelfDom=false;
    mvNormalVector.resize(3);
}
CSkinFace::~CSkinFace()
{
}
void CSkinFace::setShapeType(const uiint& shapeType)
{
    Utility::CLogger *pLogger;
    switch(shapeType) {
    case(ElementType::Quad):
        mShapeType= ElementType::Quad;
        mNumOfEdge= NumberOfEdge::Quad();
        mnOrder= ElementOrder::First;
        break;
    case(ElementType::Quad2):
        mShapeType= ElementType::Quad2;
        mNumOfEdge= NumberOfEdge::Quad();
        mnOrder= ElementOrder::Second;
        break;
    case(ElementType::Triangle):
        mShapeType= ElementType::Triangle;
        mNumOfEdge= NumberOfEdge::Triangle();
        mnOrder= ElementOrder::First;
        break;
    case(ElementType::Triangle2):
        mShapeType= ElementType::Triangle2;
        mNumOfEdge= NumberOfEdge::Triangle();
        mnOrder= ElementOrder::Second;
        break;
    case(ElementType::Beam):
        mShapeType= ElementType::Beam;
        mNumOfEdge= NumberOfEdge::Beam();
        mnOrder= ElementOrder::First;
        break;
    case(ElementType::Beam2):
        mShapeType= ElementType::Beam2;
        mNumOfEdge= NumberOfEdge::Beam();
        mnOrder= ElementOrder::Second;
        break;
    case(ElementType::Point):
        mShapeType= ElementType::Point;
        mNumOfEdge= 0;
        mnOrder= ElementOrder::Zero;
        break;
    default:
        pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"MPC Face, ShapeType Error @CSkinFace::setShapeType");
        break;
    }
    mvEdgeNode.resize(mNumOfEdge);
    mvbEdgeMarking = new bool[mNumOfEdge];
    mvEdgeFace.resize(mNumOfEdge);
    uiint iedge;
    for(iedge=0; iedge< mNumOfEdge; iedge++) {
        mvbEdgeMarking[iedge]=false;
    };
}
uiint CSkinFace::getNumOfVert()
{
    switch(mShapeType) {
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
void CSkinFace::setEdgeFace(CSkinFace* pFace, const uiint& iedge)
{
    mvEdgeFace[iedge]= pFace;
}
void CSkinFace::setEdgeConNode(CContactNode* pEdgeConNode, const uiint& iedge)
{
    mvEdgeNode[iedge]= pEdgeConNode;
    mvbEdgeMarking[iedge]= true;
}
uiint& CSkinFace::getEdgeIndex(PairConNode& pairConNode)
{
    uiint localNum0, localNum1;
    uiint nNumOfVert = getNumOfVert();
    uiint icnode;
    for(icnode=0; icnode< nNumOfVert; icnode++) {
        if(mvConNode[icnode]->getID() == pairConNode.first->getID())  localNum0=icnode;
        if(mvConNode[icnode]->getID() == pairConNode.second->getID()) localNum1=icnode;
    };
    CEdgeTree* pEdgeTree= CEdgeTree::Instance();
    switch(mShapeType) {
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
void CSkinFace::setEdgeFace(CSkinFace* pFace, PairConNode& pairConNode)
{
    uiint edgeNum= getEdgeIndex(pairConNode);
    mvEdgeFace[edgeNum]= pFace;
}
void CSkinFace::setEdgeConNode(CContactNode* pEdgeConNode, PairConNode& pairConNode)
{
    uiint edgeNum= getEdgeIndex(pairConNode);
    mvEdgeNode[edgeNum]= pEdgeConNode;
}
bool CSkinFace::isEdgeNodeMarking(PairConNode& pairConNode)
{
    uiint edgeNum= getEdgeIndex(pairConNode);
    return mvbEdgeMarking[edgeNum];
}
void CSkinFace::markingEdgeNode(PairConNode& pairConNode)
{
    uiint edgeNum= getEdgeIndex(pairConNode);
    mvbEdgeMarking[edgeNum]= true;
}
PairConNode CSkinFace::getEdgePairNode(const uiint& iedge)
{
    PairConNode pairConNode;
    Utility::CLogger *pLogger;
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    uiint* vertNum;
    switch(mShapeType) {
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
        pLogger->Info(Utility::LoggerMode::Error,"MPC Face, ShapeType Error @CSkinFace::getEdgePairNode");
        break;
    }
    pairConNode.first = mvConNode[vertNum[0]];
    pairConNode.second= mvConNode[vertNum[1]];
    return pairConNode;
}
CContactNode* CSkinFace::getEdgeConNode(PairConNode& pairConNode)
{
    uiint edgeNum= getEdgeIndex(pairConNode);
    return mvEdgeNode[edgeNum];
}
CSkinFace* CSkinFace::generateFace()
{
    mpOtherFace = new CSkinFace;
    return mpOtherFace;
}
void CSkinFace::setupNodeID_progFace(CElement* pElem, const uiint& numOfVert)
{
    CNode *pFaceNode;
    pFaceNode= pElem->getFaceNode(mElementFaceID);
    mpFaceNode->setNodeID(pFaceNode->getID());
    setupEdgeNodeID(pElem, numOfVert);
}
void CSkinFace::setupEdgeNodeID(CElement* pElem, const uiint& numOfVert)
{
    uiint ivert,nvert;
    uiint vnodeID[2];
    uiint elemEdgeIndex;
    CNode *pEdgeNode;
    for(ivert=0; ivert< numOfVert; ivert++) {
        if(ivert==numOfVert-1) {
            nvert=0;
        } else {
            nvert=ivert+1;
        }
        vnodeID[0]= mvConNode[ivert]->getNodeID();
        vnodeID[1]= mvConNode[nvert]->getNodeID();
        elemEdgeIndex= pElem->getEdgeIndex(vnodeID[0], vnodeID[1]);
        pEdgeNode= pElem->getEdgeInterNode(elemEdgeIndex);
        mvEdgeNode[ivert]->setNodeID(pEdgeNode->getID());
    };
}
void CSkinFace::setupNodeID_2nd_LastLevel(CElement* pElem)
{
    uiint numOfVert = this->getNumOfVert();
    setupEdgeNodeID(pElem, numOfVert);
}
void CSkinFace::refine(CElement *pElem, uiint& faceID)
{
    CElement *pProgElem;
    uiint nodeID;
    uiint localNum;
    uiint numOfVert;
    CNode *pEdgeNode;
    uiint iprog;
    CSkinFace* pFace;
    switch(mShapeType) {
    case(ElementType::Quad):
    case(ElementType::Quad2):
        mvProgFace.reserve(4);
        for(iprog=0; iprog< 4; iprog++) {
            pFace = generateFace();//---------------------------------- 生成
            pFace->setRank(mRank);
            pFace->setLevel(mLevel+1);
            if(mnOrder==ElementOrder::First)  pFace->resizeNode(4);
            if(mnOrder==ElementOrder::Second) pFace->resizeNode(8);
            pFace->setShapeType(mShapeType);
            pFace->setID(faceID);
            faceID++;
            if(mbSelfDom) {
                pFace->markingSelf();
                pFace->setMeshID(mMeshID);
                nodeID= mvConNode[iprog]->getNodeID();
                pProgElem= pElem->getProgElem_NodeID(nodeID);
                pFace->setElementID(pProgElem->getID());

                switch(pElem->getType()) {
                case(ElementType::Hexa):
                case(ElementType::Quad):
                case(ElementType::Hexa2):
                case(ElementType::Quad2):
                    pFace->setFaceID(mElementFaceID);
                    break;
                case(ElementType::Prism):
                case(ElementType::Prism2):
                    localNum= pElem->getLocalVertNum(nodeID);
                    switch(localNum) {
                    case(0):
                        if(mElementFaceID==2) pFace->setFaceID(2);
                        if(mElementFaceID==4) pFace->setFaceID(4);
                        break;
                    case(1):
                        if(mElementFaceID==2) pFace->setFaceID(2);
                        if(mElementFaceID==3) pFace->setFaceID(5);
                        break;
                    case(2):
                        if(mElementFaceID==3) pFace->setFaceID(3);
                        if(mElementFaceID==4) pFace->setFaceID(4);
                        break;
                    case(3):
                        if(mElementFaceID==2) pFace->setFaceID(2);
                        if(mElementFaceID==4) pFace->setFaceID(4);
                        break;
                    case(4):
                        if(mElementFaceID==2) pFace->setFaceID(2);
                        if(mElementFaceID==3) pFace->setFaceID(5);
                        break;
                    case(5):
                        if(mElementFaceID==3) pFace->setFaceID(3);
                        if(mElementFaceID==4) pFace->setFaceID(4);
                        break;
                    }
                    break;
                }
                numOfVert=4;
                setupNodeID_progFace(pElem, numOfVert);
            }
            mvProgFace.push_back(pFace);
        };
        mvProgFace[0]->setNode(mvConNode[0],0);
        mvProgFace[0]->setNode(mvEdgeNode[0],1);
        mvProgFace[0]->setNode(mpFaceNode,  2);
        mvProgFace[0]->setNode(mvEdgeNode[3],3);
        mvProgFace[1]->setNode(mvConNode[1],0);
        mvProgFace[1]->setNode(mvEdgeNode[1],1);
        mvProgFace[1]->setNode(mpFaceNode,  2);
        mvProgFace[1]->setNode(mvEdgeNode[0],3);
        mvProgFace[2]->setNode(mvConNode[2],0);
        mvProgFace[2]->setNode(mvEdgeNode[2],1);
        mvProgFace[2]->setNode(mpFaceNode,  2);
        mvProgFace[2]->setNode(mvEdgeNode[1],3);
        mvProgFace[3]->setNode(mvConNode[3],0);
        mvProgFace[3]->setNode(mvEdgeNode[3],1);
        mvProgFace[3]->setNode(mpFaceNode,  2);
        mvProgFace[3]->setNode(mvEdgeNode[2],3);
        break;
    case(ElementType::Triangle):
    case(ElementType::Triangle2):
        mvProgFace.reserve(3);
        for(iprog=0; iprog< 3; iprog++) {
            pFace = generateFace();//---------------------------------- 生成
            pFace->setRank(mRank);
            pFace->setLevel(mLevel+1);
            if(mnOrder==ElementOrder::First) {
                pFace->resizeNode(4);
                pFace->setShapeType(ElementType::Quad);
            }
            if(mnOrder==ElementOrder::Second) {
                pFace->resizeNode(8);
                pFace->setShapeType(ElementType::Quad2);
            }
            pFace->setID(faceID);
            faceID++;
            if(mbSelfDom) {
                pFace->markingSelf();
                pFace->setMeshID(mMeshID);
                nodeID= mvConNode[iprog]->getNodeID();
                pProgElem= pElem->getProgElem_NodeID(nodeID);
                pFace->setElementID(pProgElem->getID());

                switch(pElem->getType()) {
                case(ElementType::Tetra):
                case(ElementType::Tetra2):
                    localNum= pElem->getLocalVertNum(nodeID);
                    switch(localNum) {
                    case(0):
                        if(mElementFaceID==0) pFace->setFaceID(0);
                        if(mElementFaceID==1) pFace->setFaceID(2);
                        if(mElementFaceID==3) pFace->setFaceID(4);
                        break;
                    case(1):
                        if(mElementFaceID==0) pFace->setFaceID(0);
                        if(mElementFaceID==1) pFace->setFaceID(2);
                        if(mElementFaceID==2) pFace->setFaceID(5);
                        break;
                    case(2):
                        if(mElementFaceID==0) pFace->setFaceID(0);
                        if(mElementFaceID==2) pFace->setFaceID(3);
                        if(mElementFaceID==3) pFace->setFaceID(4);
                        break;
                    case(3):
                        if(mElementFaceID==1) pFace->setFaceID(2);
                        if(mElementFaceID==2) pFace->setFaceID(1);
                        if(mElementFaceID==3) pFace->setFaceID(4);
                        break;
                    }
                    break;
                case(ElementType::Triangle):
                case(ElementType::Triangle2):
                    pFace->setFaceID(mElementFaceID);
                    break;
                case(ElementType::Prism):
                case(ElementType::Prism2):
                    localNum= pElem->getLocalVertNum(nodeID);
                    switch(localNum) {
                    case(0):
                    case(1):
                    case(2):
                        pFace->setFaceID(0);
                        break;
                    case(3):
                    case(4):
                    case(5):
                        pFace->setFaceID(1);
                        break;
                    }
                    break;
                }
                numOfVert=3;
                setupNodeID_progFace(pElem, numOfVert);
            }
            mvProgFace.push_back(pFace);
        };
        mvProgFace[0]->setNode(mvConNode[0],0);
        mvProgFace[0]->setNode(mvEdgeNode[0],1);
        mvProgFace[0]->setNode(mpFaceNode,  2);
        mvProgFace[0]->setNode(mvEdgeNode[2],3);
        mvProgFace[1]->setNode(mvConNode[1],0);
        mvProgFace[1]->setNode(mvEdgeNode[1],1);
        mvProgFace[1]->setNode(mpFaceNode,  2);
        mvProgFace[1]->setNode(mvEdgeNode[0],3);
        mvProgFace[2]->setNode(mvConNode[2],0);
        mvProgFace[2]->setNode(mvEdgeNode[2],1);
        mvProgFace[2]->setNode(mpFaceNode,  2);
        mvProgFace[2]->setNode(mvEdgeNode[1],3);
        break;
    case(ElementType::Beam):
    case(ElementType::Beam2):
        mvProgFace.reserve(2);
        for(iprog=0; iprog< 2; iprog++) {
            pFace = generateFace();//---------------------------------- 生成
            pFace->setRank(mRank);
            pFace->setLevel(mLevel+1);
            if(mnOrder==ElementOrder::First) {
                pFace->resizeNode(2);
                pFace->setShapeType(ElementType::Beam);
            }
            if(mnOrder==ElementOrder::Second) {
                pFace->resizeNode(3);
                pFace->setShapeType(ElementType::Beam2);
            }
            pFace->setID(faceID);
            faceID++;
            if(mbSelfDom) {
                pFace->markingSelf();
                pFace->setMeshID(mMeshID);
                nodeID= mvConNode[iprog]->getNodeID();
                pProgElem= pElem->getProgElem_NodeID(nodeID);
                pFace->setElementID(pProgElem->getID());
                pFace->setFaceID(0);
                pEdgeNode= pElem->getEdgeInterNode(0);
                mvEdgeNode[0]->setNodeID(pEdgeNode->getID());
                mpFaceNode->setNodeID(pEdgeNode->getID());
            }
            mvProgFace.push_back(pFace);
        };
        mvProgFace[0]->setNode(mvConNode[0],0);
        mvProgFace[0]->setNode(mvEdgeNode[0],1);
        mvProgFace[1]->setNode(mvConNode[1],0);
        mvProgFace[1]->setNode(mvEdgeNode[0],1);
        break;
    }
}
vdouble& CSkinFace::CalcNzNormalVector()
{
    double X1= mvConNode[1]->getX() - mvConNode[0]->getX();
    double Y1= mvConNode[1]->getY() - mvConNode[0]->getY();
    double Z1= mvConNode[1]->getZ() - mvConNode[0]->getZ();
    double X2= mvConNode[2]->getX() - mvConNode[0]->getX();
    double Y2= mvConNode[2]->getY() - mvConNode[0]->getY();
    double Z2= mvConNode[2]->getZ() - mvConNode[0]->getZ();
    mvNormalVector[0]= Y1*Z2 - Z1*Y2;
    mvNormalVector[1]= Z1*X2 - X1*Z2;
    mvNormalVector[2]= X1*Y2 - Y1*X2;
    double vecLength= sqrt(mvNormalVector[0]*mvNormalVector[0]
                           + mvNormalVector[1]*mvNormalVector[1]
                           + mvNormalVector[2]*mvNormalVector[2]);
    mvNormalVector[0] /= vecLength;
    mvNormalVector[1] /= vecLength;
    mvNormalVector[2] /= vecLength;
    return mvNormalVector;
}
vdouble& CSkinFace::getNzNormalVector()
{
    return mvNormalVector;
}
void CSkinFace::addSlaveNode(CContactNode* pConNode)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn, "invalid method CSkinFace::addSlaveNode");
}
void CSkinFace::CalcSlave(const uiint& islave, const uiint& valType)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn, "invalid method CSkinFace::CalcSlave");
}
double& CSkinFace::getCoef(const uiint& slaveID, const uiint& ivert)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn, "invalid method CSkinFace::getCoef");
    return pLogger->getDDummyValue();
}
void CSkinFace::deleteProgData()
{
    vector<CContactNode*>().swap(mvEdgeNode);
    vector<CSkinFace*>().swap(mvEdgeFace);
    delete []mvbEdgeMarking;
}
void CSkinFace::replaceEdgeNode()
{
    if(mnOrder==ElementOrder::Second) {
        uiint nNumOfVert;
        if(mNumOfEdge==4) nNumOfVert=4;
        if(mNumOfEdge==3) nNumOfVert=3;
        if(mNumOfEdge==1) nNumOfVert=2;
        if(mNumOfEdge==0) nNumOfVert=1;
        for(uiint iedge=0; iedge < mNumOfEdge; iedge++) {
            mvConNode[nNumOfVert + iedge] = mvEdgeNode[iedge];
        };
    }
}
