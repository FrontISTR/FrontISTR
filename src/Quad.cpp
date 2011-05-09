/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Quad.cxx
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
#include "Quad.h"
using namespace pmw;
uint CQuad::mnElemType = ElementType::Quad;
uint CQuad::mnElemType2= ElementType::Quad2;
uint CQuad::mNumOfFace = 1;
uint CQuad::mNumOfEdge = 4;
uint CQuad::mNumOfNode = 4;
CQuad::CQuad()
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
CQuad::~CQuad()
{
}
const uint& CQuad::getType()
{
    if(mvEdgeInterNode.size()==0){
        return mnElemType;
    }else{
        return mnElemType2;
    }
}
bool CQuad::IndexCheck(const uint& propType, const uint& index, string& method_name)
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
        pLogger->Info(Utility::LoggerMode::Error, "CQuad::"+ method_name);
        return false;
    }else{
        return true;
    }
}
PairNode& CQuad::getPairNode(const uint& edgeIndex)
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
            mEdgeNode.second=mvNode[3];
            break;
        case(3):
            mEdgeNode.first=mvNode[3];
            mEdgeNode.second=mvNode[0];
        default:
            break;
    }
    return mEdgeNode;
}
void CQuad::getPairNode(vint& pairNodeIndex, const uint& edgeIndex)
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
            pairNodeIndex[1]= mvNode[3]->getID();
            break;
        case(3):
            pairNodeIndex[0]= mvNode[3]->getID();
            pairNodeIndex[1]= mvNode[0]->getID();
            break;
        default:
            break;
    }
}
uint& CQuad::getEdgeIndex(CNode* pNode0, CNode* pNode1)
{
    getLocalNodeNum(pNode0, pNode1);
    CEdgeTree *edgeTree = CEdgeTree::Instance();
    return edgeTree->getQuadEdgeIndex(mvPairNodeLocalNum[0], mvPairNodeLocalNum[1]);
}
uint& CQuad::getEdgeIndex(const uint& nodeID_0, const uint& nodeID_1)
{
    mvPairNodeLocalNum[0]= mmIDLocal[nodeID_0];
    mvPairNodeLocalNum[1]= mmIDLocal[nodeID_1];
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    return pEdgeTree->getQuadEdgeIndex(mvPairNodeLocalNum[0], mvPairNodeLocalNum[1]);
}
bool CQuad::isEdgeElem(CNode* pNode0, CNode* pNode1)
{
    uint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);
    return mvb_edge[edgeNum];
}
void CQuad::setBoolEdgeElem(CNode* pNode0, CNode* pNode1)
{
    uint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);
    mvb_edge[edgeNum]=true;
}
uint& CQuad::getLocalFaceNum(const vuint& vLocalNodeNum)
{
    CFaceTree *pFaceTree = CFaceTree::Instance();
    uint faceIndex = pFaceTree->getQuadFaceIndex(vLocalNodeNum);
    if(faceIndex==0){
        mTempo=0;
    }else{
        Utility::CLogger *pLogger = Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"mismatch!, Local Node number@CQuad::getLocalFaceNum");
        mTempo=1;
    }
    return mTempo;
}
uint& CQuad::getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    vuint vLocalNum; vLocalNum.resize(3);
    vLocalNum[0]= mmIDLocal[pNode0->getID()];
    vLocalNum[1]= mmIDLocal[pNode1->getID()];
    vLocalNum[2]= mmIDLocal[pNode2->getID()];
    CFaceTree *faceTree = CFaceTree::Instance();
    uint faceIndex = faceTree->getQuadFaceIndex(vLocalNum);
    if(faceIndex==0){
        mTempo=0;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"mismatch, Local Node number @CQuad::getLocalFaceNum(arg:pNode)");
        mTempo=1;
    }
    return mTempo;
}
uint& CQuad::getFaceIndex(const uint& edge0, const uint& edge1)
{
    CEdgeFaceTree *pTree = CEdgeFaceTree::Instance();
    return pTree->getQuadFaceIndex(edge0, edge1);
}
void CQuad::setupFaceCnvNodes()
{
    CNode *pNode;
    mvvFaceCnvNode.resize(1);
    mvvFaceCnvNode[0].resize(4);
    uint ivert;
    for(ivert=0; ivert< 4; ivert++){
        pNode= mvNode[ivert];
        mvvFaceCnvNode[0][ivert]=pNode;
    };
}
vector<CNode*> CQuad::getConnectNode(CNode* pNode)
{
    vector<CNode*> vNode;
    vNode.reserve(2);
    uint gID= pNode->getID();
    uint localID= mmIDLocal[gID];
    CNodeConnectNodeTree *pConnTree = CNodeConnectNodeTree::Instance();
    vuint vLocalID;
    vLocalID= pConnTree->getQuadConnectNode(localID);
    uint i;
    for(i=0; i< 2; i++){
        localID= vLocalID[i];
        vNode.push_back(mvNode[localID]);
    };
    return vNode;
}
