/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Beam.cpp
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
#include "Beam.h"
using namespace pmw;
uiint CBeam::mnElemType = ElementType::Beam;
uiint CBeam::mnElemOrder = 1;
uiint CBeam::mNumOfFace = 0;
uiint CBeam::mNumOfEdge = 1;
uiint CBeam::mNumOfNode = 2;
uiint CBeam::mNumOfVert = 2;
CBeam::CBeam()
{
    ;
}
CBeam::~CBeam()
{
}
void CBeam::initialize()
{
    mvNode.resize(mNumOfNode);
    mvvEdgeElement.resize(mNumOfEdge);
    mvEdgeInterNode.resize(mNumOfEdge);
    mvb_edge = new bool[mNumOfEdge];
    uiint i;
    for(i=0; i< mNumOfEdge; i++) {
        mvb_edge[i] = false;
    };
    mvProgElement.resize(mNumOfNode);
////    mvbCommEntity = new bool[mNumOfNode];
////    for(i=0; i< mNumOfNode; i++){
////        mvbCommEntity[i]= false;
////    };
    mvbCommLEdge = new bool[mNumOfEdge];
    for(i=0; i< mNumOfEdge; i++) {
        mvbCommLEdge[i]= false;
    };
    mvbCommLVert = new bool[mNumOfVert];
    for(i=0; i< mNumOfVert; i++) {
        mvbCommLVert[i]= false;
    }
}
const uiint& CBeam::getType()
{
    return mnElemType;
}
const uiint& CBeam::getOrder()
{
    return mnElemOrder;
}
bool CBeam::IndexCheck(const uiint& propType, const uiint& index, string& method_name)
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
        pLogger->Info(Utility::LoggerMode::Error, "CBeam::"+ method_name);
        return false;
    } else {
        return true;
    }
}
PairNode CBeam::getPairNode(const uiint& iedge)
{
    PairNode pairNode;
    switch(iedge) {
    case(0):
        pairNode.first=mvNode[0];
        pairNode.second=mvNode[1];
        break;
    default:
        break;
    }
    return pairNode;
}
void CBeam::getPairNode(vint& pairNodeIndex, const uiint& iedge)
{
    switch(iedge) {
    case(0):
        pairNodeIndex[0]= mvNode[0]->getID();
        pairNodeIndex[1]= mvNode[1]->getID();
        break;
    default:
        break;
    }
}
uiint& CBeam::getEdgeIndex(CNode* pNode0, CNode* pNode1)
{
    uiint id0, id1;
    id0 = pNode0->getID();
    id1 = pNode1->getID();
    CEdgeTree *edgeTree = CEdgeTree::Instance();
    return edgeTree->getBeamEdgeIndex(mmIDLocal[id0], mmIDLocal[id1]);
}
uiint& CBeam::getEdgeIndex(const uiint& nodeID_0, const uiint& nodeID_1)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    return pEdgeTree->getBeamEdgeIndex(mmIDLocal[nodeID_0], mmIDLocal[nodeID_1]);
}
bool CBeam::isEdgeElem(CNode* pNode0, CNode* pNode1)
{
    return mvb_edge[0];
}
void CBeam::setBoolEdgeElem(CNode* pNode0, CNode* pNode1)
{
    mvb_edge[0]=true;
}
uiint& CBeam::getLocalFaceNum(const vuint& vLocalNodeNum)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::getLocalFaceNum");
    return pLogger->getUDummyValue();
}
uiint& CBeam::getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::getLocalFaceNum(arg:pNode)");
    return pLogger->getUDummyValue();
}
uiint& CBeam::getFaceIndex(const uiint& edge0, const uiint& edge1)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::getFaceIndex");
    return pLogger->getUDummyValue();
}
bool CBeam::isFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::BoolFaceElem");
    return false;
}
void CBeam::setBoolFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::setBoolFaceElem");
}
vector<CNode*> CBeam::getFaceCnvNodes(const uiint& iface)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::setupFaceCnvNodes");
    vector<CNode*> vFaceCnvNode;
    return vFaceCnvNode;
}
vector<CNode*> CBeam::getConnectNode(CNode* pNode)
{
    vector<CNode*> vNode;
    vNode.reserve(1);
    uiint gID= pNode->getID();
    uiint localID= mmIDLocal[gID];
    CNodeConnectNodeTree *pConnTree = CNodeConnectNodeTree::Instance();
    vuint vLocalID;
    vLocalID= pConnTree->getBeamConnectNode(localID);
    uiint i;
    for(i=0; i< 1; i++) {
        localID= vLocalID[i];
        vNode.push_back(mvNode[localID]);
    };
    return vNode;
}
void CBeam::deleteProgData()
{
    uiint iedge;
    for(iedge=0; iedge < mNumOfEdge; iedge++) {
        vector<CElement*>().swap(mvvEdgeElement[iedge]);
    };
    vector<vector<CElement*> >().swap(mvvEdgeElement);
    vector<CNode*>().swap(mvEdgeInterNode);
    delete []mvb_edge;

    delete []mvbCommLEdge;
    delete []mvbCommLVert;
}
