//
//  Beam.cpp
//			2009.06.19
//			2008.12.01
//			k.Takeda
#include "Beam.h"
using namespace pmw;

uint CBeam::mnElemType = ElementType::Beam;
uint CBeam::mNumOfFace = 0;
uint CBeam::mNumOfEdge = 1;
uint CBeam::mNumOfNode = 2;
//
//
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

    //CommElementのprolongation用
    mvProgElement.resize(mNumOfNode);
}

CBeam::~CBeam()
{
//    //debug
//    cout << "~CBeam" << endl;
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



// out:(PariNode edgeNode)
//
PairNode& CBeam::getPairNode(const uint& edgeIndex)
{
    switch(edgeIndex){
        //surf 0
        case(0):
            mEdgeNode.first=mvNode[0];
            mEdgeNode.second=mvNode[1];
            break;

        default:
            break;
    }

    return mEdgeNode;
}

// out:(vint& pairNodeIndex) -> pair Node index num.
//
void CBeam::getPairNode(vint& pairNodeIndex, const uint& edgeIndex)
{
    switch(edgeIndex){
        //surf 0 (low surf)
        case(0):
            pairNodeIndex[0]= mvNode[0]->getID();
            pairNodeIndex[1]= mvNode[1]->getID();
            break;

        default:
            break;
    }
}

// ノードの局所番号に対応する、Edgeのインデックス番号
//
uint& CBeam::getEdgeIndex(CNode* pNode0, CNode* pNode1)
{
    getLocalNodeNum(pNode0, pNode1);//mvPairNodeLocalNumに値が入る.

    // Beam要素なので、”Edge_Index=0” 
    CEdgeTree *edgeTree = CEdgeTree::Instance();
    
    return edgeTree->getBeamEdgeIndex(mvPairNodeLocalNum[0],mvPairNodeLocalNum[1]);
}

// EdgeElement集合がセット済みか?
//
bool CBeam::isEdgeElem(CNode* pNode0, CNode* pNode1)
{
    return mvb_edge[0];
}

// EdgeElement集合がセットされていることをスタンプ
//
void CBeam::setBoolEdgeElem(CNode* pNode0, CNode* pNode1)
{
    mvb_edge[0]=true;// スタンプ
}

// Face局所番号 -> Beamには存在しない.
//
uint& CBeam::getLocalFaceNum(const vuint& vLocalNodeNum)
{
    mDummy= 0;

    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::getLocalFaceNum");

    return mDummy;
}
//
//
uint& CBeam::getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    mDummy= 0;

    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::getLocalFaceNum(arg:pNode)");

    return mDummy;
}


// 2 Edge => Face Index
//
uint& CBeam::getFaceIndex(const uint& edge0, const uint& edge1)
{
    mDummy= 0;

    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::getFaceIndex");

    return mDummy;
}

// 面(Face)に既に面ノードがセット済みかどうか？
//
bool CBeam::isFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::BoolFaceElem");

    return false;
}


// 面(Face)に面ノードがセットされたことをスタンプ
//
void CBeam::setBoolFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::setBoolFaceElem");
}

// mvvFaceCnvNodeのセットアップ
//
void CBeam::setupFaceCnvNodes()
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::setupFaceCnvNodes");
}


// Node周囲の接続先Nodeを配列で返す.
// (係数行列 作成用途, CMesh::setupAggregate)
//
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





