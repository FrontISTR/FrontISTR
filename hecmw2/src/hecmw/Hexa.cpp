/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Hexa.cpp
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
#include "Hexa.h"
using namespace pmw;
uiint CHexa::mnElemType = ElementType::Hexa;
uiint CHexa::mnElemOrder = ElementOrder::First;
uiint CHexa::mNumOfFace =  6;
uiint CHexa::mNumOfEdge = 12;
uiint CHexa::mNumOfNode =  8;
uiint CHexa::mNumOfVert =  8;
#include "Logger.h"
CHexa::CHexa(void):CSolidBase()
{
    ;
}
void CHexa::initialize()
{
    mvNode.resize(mNumOfNode);
    mvvEdgeElement.resize(mNumOfEdge);
    mvEdgeInterNode.resize(mNumOfEdge);
    mvb_edge = new bool[mNumOfEdge];
    uiint i;
    for(i=0; i< mNumOfEdge; i++) {
        mvb_edge[i]=false;
    };
    mvFaceElement.resize(mNumOfFace);
    mvFaceNode.resize(mNumOfFace);
    mvb_face = new bool[mNumOfFace];
    for(i=0; i< mNumOfFace; i++) {
        mvb_face[i]=false;
    };
    mvProgElement.resize(mNumOfNode);
    mvbMPCFace = new bool[mNumOfFace];
    for(i=0; i< mNumOfFace; i++) {
        mvbMPCFace[i]= false;
    };
////    mvbCommEntity = new bool[mNumOfFace];
////    for(i=0; i< mNumOfFace; i++){
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
CHexa::~CHexa(void)
{
}
const uiint& CHexa::getType()
{
    return mnElemType;
}
const uiint& CHexa::getOrder()
{
    return mnElemOrder;
}
bool CHexa::IndexCheck(const uiint& propType, const uiint& index, string& method_name)
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
        pLogger->Info(Utility::LoggerMode::Error, "CHexa::"+ method_name);
        return false;
    } else {
        return true;
    }
}
PairNode CHexa::getPairNode(const uiint& iedge)
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
        pairNode.second=mvNode[3];
        break;
    case(3):
        pairNode.first=mvNode[3];
        pairNode.second=mvNode[0];
        break;
    case(4):
        pairNode.first=mvNode[4];
        pairNode.second=mvNode[5];
        break;
    case(5):
        pairNode.first=mvNode[5];
        pairNode.second=mvNode[6];
        break;
    case(6):
        pairNode.first=mvNode[6];
        pairNode.second=mvNode[7];
        break;
    case(7):
        pairNode.first=mvNode[7];
        pairNode.second=mvNode[4];
        break;
    case(8):
        pairNode.first=mvNode[0];
        pairNode.second=mvNode[4];
        break;
    case(9):
        pairNode.first=mvNode[1];
        pairNode.second=mvNode[5];
        break;
    case(10):
        pairNode.first=mvNode[2];
        pairNode.second=mvNode[6];
        break;
    case(11):
        pairNode.first=mvNode[3];
        pairNode.second=mvNode[7];
        break;
    default:
        break;
    }
    return pairNode;
}
void CHexa::getPairNode(vint& pairNodeIndex, const uiint& iedge)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    if(pairNodeIndex.size() != 2) {
        pLogger->Info(Utility::LoggerMode::Error,"pair node size Error@CHexa::getPairNode");
        return;
    }
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
        pairNodeIndex[1]= mvNode[3]->getID();
        break;
    case(3):
        pairNodeIndex[0]= mvNode[3]->getID();
        pairNodeIndex[1]= mvNode[0]->getID();
        break;
    case(4):
        pairNodeIndex[0]= mvNode[4]->getID();
        pairNodeIndex[1]= mvNode[5]->getID();
        break;
    case(5):
        pairNodeIndex[0]= mvNode[5]->getID();
        pairNodeIndex[1]= mvNode[6]->getID();
        break;
    case(6):
        pairNodeIndex[0]= mvNode[6]->getID();
        pairNodeIndex[1]= mvNode[7]->getID();
        break;
    case(7):
        pairNodeIndex[0]= mvNode[7]->getID();
        pairNodeIndex[1]= mvNode[4]->getID();
        break;
    case(8):
        pairNodeIndex[0]= mvNode[0]->getID();
        pairNodeIndex[1]= mvNode[4]->getID();
        break;
    case(9):
        pairNodeIndex[0]= mvNode[1]->getID();
        pairNodeIndex[1]= mvNode[5]->getID();
        break;
    case(10):
        pairNodeIndex[0]= mvNode[2]->getID();
        pairNodeIndex[1]= mvNode[6]->getID();
        break;
    case(11):
        pairNodeIndex[0]= mvNode[3]->getID();
        pairNodeIndex[1]= mvNode[7]->getID();
        break;
    default:
        break;
    }
}
uiint& CHexa::getEdgeIndex(CNode* pNode0, CNode* pNode1)
{
    uiint id0, id1;
    id0 = pNode0->getID();
    id1 = pNode1->getID();
    CEdgeTree *edgeTree = CEdgeTree::Instance();
    return edgeTree->getHexaEdgeIndex(mmIDLocal[id0], mmIDLocal[id1]);
}
uiint& CHexa::getEdgeIndex(const uiint& nodeID_0, const uiint& nodeID_1)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    return pEdgeTree->getHexaEdgeIndex(mmIDLocal[nodeID_0], mmIDLocal[nodeID_1]);
}
bool CHexa::isEdgeElem(CNode* pNode0, CNode* pNode1)
{
    uiint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);
    return mvb_edge[edgeNum];
}
void CHexa::setBoolEdgeElem(CNode* pNode0, CNode* pNode1)
{
    uiint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);
    mvb_edge[edgeNum]=true;
}
uiint& CHexa::getLocalFaceNum(const vuint& vLocalNodeNum)
{
    CFaceTree *faceTree = CFaceTree::Instance();
    return faceTree->getHexaFaceIndex2(vLocalNodeNum);
}
uiint& CHexa::getFaceIndex(const uiint& edge0, const uiint& edge1)
{
    CEdgeFaceTree *faceTree= CEdgeFaceTree::Instance();
    return faceTree->getHexaFaceIndex(edge0, edge1);
}
uiint& CHexa::getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    vuint vLocalNum;
    vLocalNum.resize(3);
    vLocalNum[0]= mmIDLocal[pNode0->getID()];
    vLocalNum[1]= mmIDLocal[pNode1->getID()];
    vLocalNum[2]= mmIDLocal[pNode2->getID()];
    CFaceTree *faceTree = CFaceTree::Instance();
    return faceTree->getHexaFaceIndex2(vLocalNum);
}
vector<CNode*> CHexa::getFaceCnvNodes(const uiint& iface)
{
    CFaceTree *faceTree = CFaceTree::Instance();
    uiint* faceConnectivity;
    CNode *pNode;
    vector<CNode*> vFaceCnvNode;
    uiint ivert, index;
    vFaceCnvNode.resize(4);
    faceConnectivity= faceTree->getLocalNodeHexaFace(iface);
    for(ivert=0; ivert< 4; ivert++) {
        index= faceConnectivity[ivert];
        pNode= mvNode[index];
        vFaceCnvNode[ivert]=pNode;
    };
    return vFaceCnvNode;
}
vector<CNode*> CHexa::getConnectNode(CNode* pNode)
{
    vector<CNode*> vNode;
    vNode.reserve(3);
    uiint gID= pNode->getID();
    uiint localID= mmIDLocal[gID];
    CNodeConnectNodeTree *pConnTree = CNodeConnectNodeTree::Instance();
    vuint vLocalID;
    vLocalID= pConnTree->getHexaConnectNode(localID);
    uiint i;
    for(i=0; i< 3; i++) {
        localID= vLocalID[i];
        vNode.push_back(mvNode[localID]);
    };
    return vNode;
}
void CHexa::deleteProgData()
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
