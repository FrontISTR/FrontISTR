/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Tetra.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Tetra.h"
using namespace pmw;
uint CTetra::mnElemType = ElementType::Tetra;
uint CTetra::mnElemType2= ElementType::Tetra2;
uint CTetra::mNumOfFace = 4;
uint CTetra::mNumOfEdge = 6;
uint CTetra::mNumOfNode = 4;
CTetra::CTetra(void)
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
    mvbCommEntity.resize(mNumOfFace);
    for(i=0; i< mNumOfFace; i++){
        mvbCommEntity[i]= false;
    };
}
CTetra::~CTetra(void)
{
}
const uint& CTetra::getType()
{
    if(mvEdgeInterNode.size()==0){
        return mnElemType;
    }else{
        return mnElemType2;
    }
}
bool CTetra::IndexCheck(const uint& propType, const uint& index, string& method_name)
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
        pLogger->Info(Utility::LoggerMode::Error, "CTetra::"+ method_name);
        return false;
    }else{
        return true;
    }
}
PairNode& CTetra::getPairNode(const uint& edgeIndex)
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
        case(3):
            mEdgeNode.first=mvNode[0];
            mEdgeNode.second=mvNode[3];
            break;
        case(4):
            mEdgeNode.first=mvNode[1];
            mEdgeNode.second=mvNode[3];
            break;
        case(5):
            mEdgeNode.first=mvNode[2];
            mEdgeNode.second=mvNode[3];
            break;
        default:
            break;
    }
    return mEdgeNode;
}
void CTetra::getPairNode(vint& pairNodeIndex, const uint& edgeIndex)
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
        case(3):
            pairNodeIndex[0]= mvNode[0]->getID();
            pairNodeIndex[1]= mvNode[3]->getID();
            break;
        case(4):
            pairNodeIndex[0]= mvNode[1]->getID();
            pairNodeIndex[1]= mvNode[3]->getID();
            break;
        case(5):
            pairNodeIndex[0]= mvNode[2]->getID();
            pairNodeIndex[1]= mvNode[3]->getID();
            break;
        default:
            break;
    }
}
uint& CTetra::getEdgeIndex(CNode* pNode0, CNode* pNode1)
{
    getLocalNodeNum(pNode0, pNode1);
    CEdgeTree *edgeTree = CEdgeTree::Instance();
    return edgeTree->getTetraEdgeIndex(mvPairNodeLocalNum[0], mvPairNodeLocalNum[1]);
}
uint& CTetra::getEdgeIndex(const uint& nodeID_0, const uint& nodeID_1)
{
    mvPairNodeLocalNum[0]= mmIDLocal[nodeID_0];
    mvPairNodeLocalNum[1]= mmIDLocal[nodeID_1];
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    return pEdgeTree->getTetraEdgeIndex(mvPairNodeLocalNum[0], mvPairNodeLocalNum[1]);
}
bool CTetra::isEdgeElem(CNode* pNode0, CNode* pNode1)
{
    uint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);
    return mvb_edge[edgeNum];
}
void CTetra::setBoolEdgeElem(CNode* pNode0, CNode* pNode1)
{
    uint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);
    mvb_edge[edgeNum]=true;
}
uint& CTetra::getLocalFaceNum(const vuint& vLocalNodeNum)
{
    CFaceTree *pFaceTree = CFaceTree::Instance();
    return pFaceTree->getTetraFaceIndex2(vLocalNodeNum);
}
uint& CTetra::getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    vuint vLocalNum; vLocalNum.resize(3);
    vLocalNum[0]= mmIDLocal[pNode0->getID()];
    vLocalNum[1]= mmIDLocal[pNode1->getID()];
    vLocalNum[2]= mmIDLocal[pNode2->getID()];
    CFaceTree *faceTree = CFaceTree::Instance();
    return faceTree->getTetraFaceIndex2(vLocalNum);
}
uint& CTetra::getFaceIndex(const uint& edge0, const uint& edge1)
{
    CEdgeFaceTree *pEdgeFace = CEdgeFaceTree::Instance();
    return pEdgeFace->getTetraFaceIndex(edge0, edge1);
}
void CTetra::setupFaceCnvNodes()
{
    CFaceTree *faceTree = CFaceTree::Instance();
    uint* faceConnectivity;
    CNode *pNode;
    mvvFaceCnvNode.resize(4);
    uint iface, ivert, index;
    for(iface=0; iface< 4; iface++){
        mvvFaceCnvNode[iface].resize(3);
        faceConnectivity= faceTree->getLocalNodeTetraFace(iface);
        for(ivert=0; ivert< 3; ivert++){
            index= faceConnectivity[ivert];
            pNode= mvNode[index];
            mvvFaceCnvNode[iface][ivert]=pNode;
        };
    };
}
vector<CNode*> CTetra::getConnectNode(CNode* pNode)
{
    vector<CNode*> vNode;
    vNode.reserve(3);
    uint gID= pNode->getID();
    uint localID= mmIDLocal[gID];
    CNodeConnectNodeTree *pConnTree = CNodeConnectNodeTree::Instance();
    vuint vLocalID;
    vLocalID= pConnTree->getTetraConnectNode(localID);
    uint i;
    for(i=0; i< 3; i++){
        localID= vLocalID[i];
        vNode.push_back(mvNode[localID]);
    };
    return vNode;
}
