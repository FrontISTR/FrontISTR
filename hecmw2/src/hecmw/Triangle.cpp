/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Triangle.cpp
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
#include "FaceTree.h"
#include "Triangle.h"
using namespace pmw;
uiint CTriangle::mnElemType = ElementType::Triangle;
uiint CTriangle::mnElemOrder = 1;
uiint CTriangle::mNumOfFace = 1;
uiint CTriangle::mNumOfEdge = 3;
uiint CTriangle::mNumOfNode = 3;
uiint CTriangle::mNumOfVert = 3;
CTriangle::CTriangle()
{
    ;
}
CTriangle::~CTriangle()
{
}
void CTriangle::initialize()
{
    mvNode.resize(mNumOfNode);
    mvvEdgeElement.resize(mNumOfEdge);
    mvEdgeInterNode.resize(mNumOfEdge);
    mvb_edge = new bool[mNumOfEdge];
    uiint i;
    for(i=0; i< mNumOfEdge; i++) {
        mvb_edge[i] = false;
    };
    mvFaceElement.resize(mNumOfFace);
    mvFaceNode.resize(mNumOfFace);
    mvb_face = new bool[mNumOfFace];
    for(i=0; i< mNumOfFace; i++) {
        mvb_face[i] = false;
    };
    mvProgElement.resize(mNumOfNode);
    mvbMPCFace = new bool[mNumOfFace];
    for(i=0; i< mNumOfFace; i++) {
        mvbMPCFace[i]= false;
    };
////    mvbCommEntity = new bool[mNumOfEdge];
////    for(i=0; i< mNumOfEdge; i++){
////        mvbCommEntity[i]= false;
////    };
    mvbCommLFace = new bool[mNumOfFace];
    for(i=0; i< mNumOfFace; i++) {
        mvbCommLFace[i]= false;
    };
    mvbCommLEdge = new bool[mNumOfEdge];
    for(i=0; i< mNumOfEdge; i++) {
        mvbCommLEdge[i]= false;
    };
    mvbCommLVert = new bool[mNumOfVert];
    for(i=0; i< mNumOfVert; i++) {
        mvbCommLVert[i]= false;
    }
}
const uiint& CTriangle::getType()
{
    return mnElemType;
}
const uiint& CTriangle::getOrder()
{
    return mnElemOrder;
}
bool CTriangle::IndexCheck(const uiint& propType, const uiint& index, string& method_name)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    uiint numOfProp;
    switch(propType) {
    case(ElementPropType::Face):
        numOfProp = mNumOfFace;
        break;
    case(ElementPropType::Edge):
        numOfProp = mNumOfEdge;
        break;
    case(ElementPropType::Node):
        numOfProp = mNumOfNode;
        break;
    default:
        numOfProp = 0;
        break;
    }
    if(index >= numOfProp) {
        pLogger->Info(Utility::LoggerMode::Error, "CTriangle::"+ method_name);
        return false;
    } else {
        return true;
    }
}
PairNode CTriangle::getPairNode(const uiint& iedge)
{
    PairNode pairNode;
    switch(iedge) {
    case(0):
        pairNode.first=mvNode[0];
        pairNode.second=mvNode[1];
        break;
    case(1):
        pairNode.first=mvNode[1];
        pairNode.second=mvNode[2];
        break;
    case(2):
        pairNode.first=mvNode[2];
        pairNode.second=mvNode[0];
        break;
    default:
        break;
    }
    return pairNode;
}
void CTriangle::getPairNode(vint& pairNodeIndex, const uiint& iedge)
{
    switch(iedge) {
    case(0):
        pairNodeIndex[0]= mvNode[0]->getID();
        pairNodeIndex[1]= mvNode[1]->getID();
        break;
    case(1):
        pairNodeIndex[0]= mvNode[1]->getID();
        pairNodeIndex[1]= mvNode[2]->getID();
        break;
    case(2):
        pairNodeIndex[0]= mvNode[2]->getID();
        pairNodeIndex[1]= mvNode[0]->getID();
        break;
    default:
        break;
    }
}
uiint& CTriangle::getEdgeIndex(CNode* pNode0, CNode* pNode1)
{
    uiint id0, id1;
    id0 = pNode0->getID();
    id1 = pNode1->getID();
    CEdgeTree *edgeTree = CEdgeTree::Instance();
    return edgeTree->getTriangleEdgeIndex(mmIDLocal[id0], mmIDLocal[id1]);
}
uiint& CTriangle::getEdgeIndex(const uiint& nodeID_0, const uiint& nodeID_1)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    return pEdgeTree->getTriangleEdgeIndex(mmIDLocal[nodeID_0], mmIDLocal[nodeID_1]);
}
bool CTriangle::isEdgeElem(CNode* pNode0, CNode* pNode1)
{
    uiint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);
    return mvb_edge[edgeNum];
}
void CTriangle::setBoolEdgeElem(CNode* pNode0, CNode* pNode1)
{
    uiint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);
    mvb_edge[edgeNum]=true;
}
uiint& CTriangle::getLocalFaceNum(const vuint& vLocalNodeNum)
{
    CFaceTree *pFaceTree = CFaceTree::Instance();
    uiint faceIndex = pFaceTree->getTriangleFaceIndex(vLocalNodeNum);
    if(faceIndex==0) {
        mTempo= 0;
    } else {
        Utility::CLogger *pLogger = Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"mismatch!, Local Node number@CTriangle::getLocalFaceNum");
        mTempo= 1;
    }
    return mTempo;
}
uiint& CTriangle::getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    vuint vLocalNum;
    vLocalNum.resize(3);
    vLocalNum[0]= 0;
    vLocalNum[1]= 1;
    vLocalNum[2]= 2;
    CFaceTree *faceTree = CFaceTree::Instance();
    uiint faceIndex = faceTree->getTriangleFaceIndex(vLocalNum);
    if(faceIndex==0) {
        mTempo= 0;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"mismatch!, Local Node Number @Triangle::getLocalFaceNum(arg:pNode)");
        mTempo= 1;
    }
    return mTempo;
}
uiint& CTriangle::getFaceIndex(const uiint& edge0, const uiint& edge1)
{
    CEdgeFaceTree *pTree= CEdgeFaceTree::Instance();
    return pTree->getTriangleFaceIndex(edge0, edge1);
}
vector<CNode*> CTriangle::getFaceCnvNodes(const uiint& iface)
{
    CNode *pNode;
    vector<CNode*> vFaceCnvNode;
    vFaceCnvNode.resize(3);
    uiint ivert;
    for(ivert=0; ivert< 3; ivert++) {
        pNode= mvNode[ivert];
        vFaceCnvNode[ivert]=pNode;
    };
    return vFaceCnvNode;
}
vector<CNode*> CTriangle::getConnectNode(CNode* pNode)
{
    vector<CNode*> vNode;
    vNode.reserve(2);
    uiint gID= pNode->getID();
    uiint localID= mmIDLocal[gID];
    CNodeConnectNodeTree *pConnTree = CNodeConnectNodeTree::Instance();
    vuint vLocalID;
    vLocalID= pConnTree->getTriangleConnectNode(localID);
    uiint i;
    for(i=0; i< 2; i++) {
        localID= vLocalID[i];
        vNode.push_back(mvNode[localID]);
    };
    return vNode;
}
void CTriangle::deleteProgData()
{
    uiint iedge;
    for(iedge=0; iedge < mNumOfEdge; iedge++) {
        vector<CElement*>().swap(mvvEdgeElement[iedge]);
    };
    vector<vector<CElement*> >().swap(mvvEdgeElement);
    vector<CNode*>().swap(mvEdgeInterNode);
    delete []mvb_edge;
    vector<CElement*>().swap(mvFaceElement);
    vector<CNode*>().swap(mvFaceNode);
    delete []mvb_face;
    delete []mvbMPCFace;
    delete []mvbCommLFace;
    delete []mvbCommLEdge;
    delete []mvbCommLVert;
}
