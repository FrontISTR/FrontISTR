/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Triangle.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FaceTree.h"
#include "Triangle.h"
using namespace pmw;
uint CTriangle::mnElemType = ElementType::Triangle;
uint CTriangle::mnElemType2= ElementType::Triangle2;
uint CTriangle::mNumOfFace = 1;
uint CTriangle::mNumOfEdge = 3;
uint CTriangle::mNumOfNode = 3;
CTriangle::CTriangle()
{
    mvNode.resize(mNumOfNode);
    mvvEdgeElement.resize(mNumOfEdge);
    mvEdgeInterNode.resize(mNumOfEdge);
    mvb_edge.reserve(mNumOfEdge);
    uint i;
    for(i=0; i< mNumOfEdge; i++){
      mvb_edge.push_back(false);
    };
    mvFaceElement.resize(mNumOfFace);
    mvFaceNode.resize(mNumOfFace);
    mvb_face.reserve(mNumOfFace);
    for(i=0; i< mNumOfFace; i++){
        mvb_face.push_back(false);
    };
    mvProgElement.resize(mNumOfNode);
    mvbMPCFace.resize(mNumOfFace);
    for(i=0; i< mNumOfFace; i++){
        mvbMPCFace[i]= false;
    };
    mvbCommEntity.resize(mNumOfEdge);
    for(i=0; i< mNumOfEdge; i++){
        mvbCommEntity[i]= false;
    };
}
CTriangle::~CTriangle()
{
}
const uint& CTriangle::getType()
{
    if(mvEdgeInterNode.size()==0){
        return mnElemType;
    }else{
        return mnElemType2;
    }
}
bool CTriangle::IndexCheck(const uint& propType, const uint& index, string& method_name)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    uint numOfProp;
    switch(propType){
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
    if(index >= numOfProp){
        pLogger->Info(Utility::LoggerMode::Error, "CTriangle::"+ method_name);
        return false;
    }else{
        return true;
    }
}
PairNode& CTriangle::getPairNode(const uint& edgeIndex)
{
    switch(edgeIndex){
        case(0):
            mEdgeNode.first=mvNode[0];
            mEdgeNode.second=mvNode[1];
            break;
        case(1):
            mEdgeNode.first=mvNode[1];
            mEdgeNode.second=mvNode[2];
            break;
        case(2):
            mEdgeNode.first=mvNode[2];
            mEdgeNode.second=mvNode[0];
            break;
        default:
            break;
    }
    return mEdgeNode;
}
void CTriangle::getPairNode(vint& pairNodeIndex, const uint& edgeIndex)
{
    switch(edgeIndex){
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
uint& CTriangle::getEdgeIndex(CNode* pNode0, CNode* pNode1)
{ 
    getLocalNodeNum(pNode0, pNode1);
    CEdgeTree *edgeTree = CEdgeTree::Instance();
    return edgeTree->getTriangleEdgeIndex(mvPairNodeLocalNum[0], mvPairNodeLocalNum[1]);
}
uint& CTriangle::getEdgeIndex(const uint& nodeID_0, const uint& nodeID_1)
{
    mvPairNodeLocalNum[0]= mmIDLocal[nodeID_0];
    mvPairNodeLocalNum[1]= mmIDLocal[nodeID_1];
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    return pEdgeTree->getTriangleEdgeIndex(mvPairNodeLocalNum[0], mvPairNodeLocalNum[1]);
}
bool CTriangle::isEdgeElem(CNode* pNode0, CNode* pNode1)
{
    uint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);
    return mvb_edge[edgeNum];
}
void CTriangle::setBoolEdgeElem(CNode* pNode0, CNode* pNode1)
{
    uint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);
    mvb_edge[edgeNum]=true;
}
uint& CTriangle::getLocalFaceNum(const vuint& vLocalNodeNum)
{
    CFaceTree *pFaceTree = CFaceTree::Instance();
    uint faceIndex = pFaceTree->getTriangleFaceIndex(vLocalNodeNum);
    if(faceIndex==0){
        mTempo= 0;
    }else{
        Utility::CLogger *pLogger = Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"mismatch!, Local Node number@CTriangle::getLocalFaceNum");
        mTempo= 1;
    }
    return mTempo;
}
uint& CTriangle::getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    vuint vLocalNum; vLocalNum.resize(3);
    vLocalNum[0]= 0;
    vLocalNum[1]= 1;
    vLocalNum[2]= 2;
    CFaceTree *faceTree = CFaceTree::Instance();
    uint faceIndex = faceTree->getTriangleFaceIndex(vLocalNum);
    if(faceIndex==0){
        mTempo= 0;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"mismatch!, Local Node Number @Triangle::getLocalFaceNum(arg:pNode)");
        mTempo= 1;
    }
    return mTempo;
}
uint& CTriangle::getFaceIndex(const uint& edge0, const uint& edge1)
{
    CEdgeFaceTree *pTree= CEdgeFaceTree::Instance();
    return pTree->getTriangleFaceIndex(edge0, edge1);
}
void CTriangle::setupFaceCnvNodes()
{
    CNode *pNode;
    mvvFaceCnvNode.resize(1);
    mvvFaceCnvNode[0].resize(3);
    uint ivert;
    for(ivert=0; ivert< 3; ivert++){
        pNode= mvNode[ivert];
        mvvFaceCnvNode[0][ivert]=pNode;
    };
}
vector<CNode*> CTriangle::getConnectNode(CNode* pNode)
{
    vector<CNode*> vNode;
    vNode.reserve(2);
    uint gID= pNode->getID();
    uint localID= mmIDLocal[gID];
    CNodeConnectNodeTree *pConnTree = CNodeConnectNodeTree::Instance();
    vuint vLocalID;
    vLocalID= pConnTree->getTriangleConnectNode(localID);
    uint i;
    for(i=0; i< 2; i++){
        localID= vLocalID[i];
        vNode.push_back(mvNode[localID]);
    };
    return vNode;
}
