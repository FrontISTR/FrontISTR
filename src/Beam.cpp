/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Beam.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Beam.h"
using namespace pmw;
uint CBeam::mnElemType = ElementType::Beam;
uint CBeam::mnElemType2= ElementType::Beam2;
uint CBeam::mNumOfFace = 0;
uint CBeam::mNumOfEdge = 1;
uint CBeam::mNumOfNode = 2;
CBeam::CBeam()
{
    mvNode.resize(mNumOfNode);
    mvvEdgeElement.resize(mNumOfEdge);
    mvEdgeInterNode.resize(mNumOfEdge);
    mvb_edge.reserve(mNumOfEdge);
    uint i;
    for(i=0; i< mNumOfEdge; i++){
      mvb_edge.push_back(false);
    };
    mvProgElement.resize(mNumOfNode);
    mvbCommEntity.resize(mNumOfNode);
    for(i=0; i< mNumOfNode; i++){
        mvbCommEntity[i]= false;
    };
}
CBeam::~CBeam()
{
}
const uint& CBeam::getType()
{
    if(mvEdgeInterNode.size()==0){
        return mnElemType;
    }else{
        return mnElemType2;
    }
}
bool CBeam::IndexCheck(const uint& propType, const uint& index, string& method_name)
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
        pLogger->Info(Utility::LoggerMode::Error, "CBeam::"+ method_name);
        return false;
    }else{
        return true;
    }
}
PairNode& CBeam::getPairNode(const uint& edgeIndex)
{
    switch(edgeIndex){
        case(0):
            mEdgeNode.first=mvNode[0];
            mEdgeNode.second=mvNode[1];
            break;
        default:
            break;
    }
    return mEdgeNode;
}
void CBeam::getPairNode(vint& pairNodeIndex, const uint& edgeIndex)
{
    switch(edgeIndex){
        case(0):
            pairNodeIndex[0]= mvNode[0]->getID();
            pairNodeIndex[1]= mvNode[1]->getID();
            break;
        default:
            break;
    }
}
uint& CBeam::getEdgeIndex(CNode* pNode0, CNode* pNode1)
{
    getLocalNodeNum(pNode0, pNode1);
    CEdgeTree *edgeTree = CEdgeTree::Instance();
    return edgeTree->getBeamEdgeIndex(mvPairNodeLocalNum[0],mvPairNodeLocalNum[1]);
}
uint& CBeam::getEdgeIndex(const uint& nodeID_0, const uint& nodeID_1)
{
    mvPairNodeLocalNum[0]= mmIDLocal[nodeID_0];
    mvPairNodeLocalNum[1]= mmIDLocal[nodeID_1];
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    return pEdgeTree->getBeamEdgeIndex(mvPairNodeLocalNum[0], mvPairNodeLocalNum[1]);
}
bool CBeam::isEdgeElem(CNode* pNode0, CNode* pNode1)
{
    return mvb_edge[0];
}
void CBeam::setBoolEdgeElem(CNode* pNode0, CNode* pNode1)
{
    mvb_edge[0]=true;
}
uint& CBeam::getLocalFaceNum(const vuint& vLocalNodeNum)
{
    mDummy= 0;
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::getLocalFaceNum");
    return mDummy;
}
uint& CBeam::getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    mDummy= 0;
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::getLocalFaceNum(arg:pNode)");
    return mDummy;
}
uint& CBeam::getFaceIndex(const uint& edge0, const uint& edge1)
{
    mDummy= 0;
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::getFaceIndex");
    return mDummy;
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
void CBeam::setupFaceCnvNodes()
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::setupFaceCnvNodes");
}
vector<CNode*> CBeam::getConnectNode(CNode* pNode)
{
    vector<CNode*> vNode;
    vNode.reserve(1);
    uint gID= pNode->getID();
    uint localID= mmIDLocal[gID];
    CNodeConnectNodeTree *pConnTree = CNodeConnectNodeTree::Instance();
    vuint vLocalID;
    vLocalID= pConnTree->getBeamConnectNode(localID);
    uint i;
    for(i=0; i< 1; i++){
        localID= vLocalID[i];
        vNode.push_back(mvNode[localID]);
    };
    return vNode;
}
