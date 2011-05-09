/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Pyramid.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Pyramid.h"
using namespace pmw;
uint CPyramid::mnElemType = ElementType::Pyramid;
uint CPyramid::mnElemType2= ElementType::Pyramid2;
uint CPyramid::mNumOfFace = 5;
uint CPyramid::mNumOfEdge = 8;
uint CPyramid::mNumOfNode = 5;
CPyramid::CPyramid()
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
    mvProgElement.resize(8);
    mvbMPCFace.resize(mNumOfFace);
    for(i=0; i< mNumOfFace; i++){
        mvbMPCFace[i]= false;
    };
    mvbCommEntity.resize(mNumOfFace);
    for(i=0; i< mNumOfFace; i++){
        mvbCommEntity[i]= false;
    };
}
CPyramid::~CPyramid()
{
}
const uint& CPyramid::getType()
{
    if(mvEdgeInterNode.size()==0){
        return mnElemType;
    }else{
        return mnElemType2;
    }
}
bool CPyramid::IndexCheck(const uint& propType, const uint& index, string& method_name)
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
        pLogger->Info(Utility::LoggerMode::Error, "CPyramid::"+ method_name);
        return false;
    }else{
        return true;
    }
}
PairNode& CPyramid::getPairNode(const uint& edgeIndex)
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
            break;
        case(4):
            mEdgeNode.first=mvNode[4];
            mEdgeNode.second=mvNode[1];
            break;
        case(5):
            mEdgeNode.first=mvNode[4];
            mEdgeNode.second=mvNode[2];
            break;
        case(6):
            mEdgeNode.first=mvNode[4];
            mEdgeNode.second=mvNode[3];
            break;
        case(7):
            mEdgeNode.first=mvNode[4];
            mEdgeNode.second=mvNode[0];
            break;
        default:
            break;
    }
    return mEdgeNode;
}
void CPyramid::getPairNode(vint& pairNodeIndex, const uint& edgeIndex)
{
    switch(edgeIndex){
        case(0):
            pairNodeIndex[0]=mvNode[0]->getID();
            pairNodeIndex[1]=mvNode[1]->getID();
            break;
        case(1):
            pairNodeIndex[0]=mvNode[1]->getID();
            pairNodeIndex[1]=mvNode[2]->getID();
            break;
        case(2):
            pairNodeIndex[0]=mvNode[2]->getID();
            pairNodeIndex[1]=mvNode[3]->getID();
            break;
        case(3):
            pairNodeIndex[0]=mvNode[3]->getID();
            pairNodeIndex[1]=mvNode[0]->getID();
            break;
        case(4):
            pairNodeIndex[0]=mvNode[4]->getID();
            pairNodeIndex[1]=mvNode[1]->getID();
            break;
        case(5):
            pairNodeIndex[0]=mvNode[4]->getID();
            pairNodeIndex[1]=mvNode[2]->getID();
            break;
        case(6):
            pairNodeIndex[0]=mvNode[4]->getID();
            pairNodeIndex[1]=mvNode[3]->getID();
            break;
        case(7):
            pairNodeIndex[0]=mvNode[4]->getID();
            pairNodeIndex[1]=mvNode[0]->getID();
            break;
        default:
            break;
    }
}
uint& CPyramid::getEdgeIndex(CNode* pNode0, CNode* pNode1)
{
    getLocalNodeNum(pNode0, pNode1);
    CEdgeTree *edgeTree = CEdgeTree::Instance();
    return edgeTree->getPyramidEdgeIndex(mvPairNodeLocalNum[0], mvPairNodeLocalNum[1]);
}
uint& CPyramid::getEdgeIndex(const uint& nodeID_0, const uint& nodeID_1)
{
    mvPairNodeLocalNum[0]= mmIDLocal[nodeID_0];
    mvPairNodeLocalNum[1]= mmIDLocal[nodeID_1];
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    return pEdgeTree->getPyramidEdgeIndex(mvPairNodeLocalNum[0], mvPairNodeLocalNum[1]);
}
bool CPyramid::isEdgeElem(CNode* pNode0, CNode* pNode1)
{
    uint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);
    return mvb_edge[edgeNum];
}
void CPyramid::setBoolEdgeElem(CNode* pNode0, CNode* pNode1)
{
    uint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);
    mvb_edge[edgeNum]=true;
}
uint& CPyramid::getLocalFaceNum(const vuint& vLocalNodeNum)
{
    CFaceTree *pFaceTree = CFaceTree::Instance();
    return  pFaceTree->getPyramidFaceIndex2(vLocalNodeNum);
}
uint& CPyramid::getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    vuint vLocalNum; vLocalNum.reserve(3);
    vLocalNum[0]= mmIDLocal[pNode0->getID()];
    vLocalNum[1]= mmIDLocal[pNode1->getID()];
    vLocalNum[2]= mmIDLocal[pNode2->getID()];
    CFaceTree *faceTree = CFaceTree::Instance();
    return faceTree->getPyramidFaceIndex2(vLocalNum);
}
uint& CPyramid::getFaceIndex(const uint& edge0, const uint& edge1)
{
    CEdgeFaceTree *pTree= CEdgeFaceTree::Instance();
    return pTree->getPyramidFaceIndex(edge0, edge1);
}
void CPyramid::setupFaceCnvNodes()
{
    CFaceTree *faceTree = CFaceTree::Instance();
    uint* faceConnectivity;
    CNode *pNode;
    mvvFaceCnvNode.resize(5);
    uint iface, ivert, index;
    uint numOfVert;
    for(iface=0; iface< 5; iface++){
        if(iface==0){
            numOfVert= 4;
        }else{
            numOfVert= 3;
        }
        mvvFaceCnvNode[iface].resize(numOfVert);
        faceConnectivity= faceTree->getLocalNodePyramidFace(iface);
        for(ivert=0; ivert< numOfVert; ivert++){
            index= faceConnectivity[ivert];
            pNode= mvNode[index];
            mvvFaceCnvNode[iface][ivert]=pNode;
        };
    };
}
vector<CNode*> CPyramid::getConnectNode(CNode* pNode)
{
    uint gID= pNode->getID();
    uint localID= mmIDLocal[gID];
    CNodeConnectNodeTree *pConnTree = CNodeConnectNodeTree::Instance();
    vuint vLocalID;
    vLocalID= pConnTree->getPyramidConnectNode(localID);
    vector<CNode*> vNode;
    vNode.reserve(vLocalID.size());
    uint i;
    for(i=0; i< vLocalID.size(); i++){
        localID= vLocalID[i];
        vNode.push_back(mvNode[localID]);
    };
    return vNode;
}
