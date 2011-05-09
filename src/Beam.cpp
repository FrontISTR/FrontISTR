//
//  Beam.cpp
//			2009.06.19
//			2008.12.01
//			k.Takeda
#include <vector>

#include "Beam.h"
using namespace pmw;

uint CBeam::mnElemType = ElementType::Beam;
uint CBeam::mnElemOrder = 1;
uint CBeam::mNumOfFace = 0;
uint CBeam::mNumOfEdge = 1;
uint CBeam::mNumOfNode = 2;
uint CBeam::mNumOfVert = 2;
//
//
CBeam::CBeam()
{
    ;
}
CBeam::~CBeam()
{
//    //debug
//    cout << "~CBeam" << endl;
}
void CBeam::initialize()
{
    mvNode.resize(mNumOfNode);

    mvvEdgeElement.resize(mNumOfEdge);

    mvEdgeInterNode.resize(mNumOfEdge);

    mvb_edge = new bool[mNumOfEdge];
    uint i;
    for(i=0; i< mNumOfEdge; i++){
      mvb_edge[i] = false;
    };

    //CommElementのprolongation用
    mvProgElement.resize(mNumOfNode);


    //通信界面(CommMesh2)属性:点
    mvbCommEntity = new bool[mNumOfNode];
    for(i=0; i< mNumOfNode; i++){
        mvbCommEntity[i]= false;
    };
}

const uint& CBeam::getType()
{
    return mnElemType;
}
const uint& CBeam::getOrder()
{
    return mnElemOrder;
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
PairNode CBeam::getPairNode(const uint& iedge)
{
    PairNode pairNode;

    switch(iedge){
        //surf 0
        case(0):
            pairNode.first=mvNode[0];
            pairNode.second=mvNode[1];
            break;

        default:
            break;
    }

    return pairNode;
}

// out:(vint& pairNodeIndex) -> pair Node index num.
//
void CBeam::getPairNode(vint& pairNodeIndex, const uint& iedge)
{
    switch(iedge){
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
    uint id0, id1;//MeshでのノードのIndex番号
    id0 = pNode0->getID();
    id1 = pNode1->getID();
    
    // Beam要素なので,”Edge_Index=0”
    CEdgeTree *edgeTree = CEdgeTree::Instance();
    
    return edgeTree->getBeamEdgeIndex(mmIDLocal[id0], mmIDLocal[id1]);
}
uint& CBeam::getEdgeIndex(const uint& nodeID_0, const uint& nodeID_1)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();
    // Beam要素なので,”Edge_Index=0” だが,参照で返したいので,わざわざTreeを呼んでいる.
    //
    return pEdgeTree->getBeamEdgeIndex(mmIDLocal[nodeID_0], mmIDLocal[nodeID_1]);
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
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::getLocalFaceNum");

    return pLogger->getUDummyValue();
}
//
//
uint& CBeam::getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::getLocalFaceNum(arg:pNode)");

    return pLogger->getUDummyValue();
}


// 2 Edge => Face Index
//
uint& CBeam::getFaceIndex(const uint& edge0, const uint& edge1)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::getFaceIndex");

    return pLogger->getUDummyValue();
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
vector<CNode*> CBeam::getFaceCnvNodes(const uint& iface)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::setupFaceCnvNodes");

    vector<CNode*> vFaceCnvNode;

    return vFaceCnvNode;
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


//
// 1. 辺-面 Element*配列を解放
// 2. 辺-面 Node*   配列を解放 (2次要素は辺ノードを残す)
//
void CBeam::deleteProgData()
{
    // Edge
    uint iedge;
    for(iedge=0; iedge < mNumOfEdge; iedge++){
        vector<CElement*>().swap(mvvEdgeElement[iedge]);
    };
    vector<vector<CElement*> >().swap(mvvEdgeElement);
    vector<CNode*>().swap(mvEdgeInterNode);//2次要素の場合は残すこと
    delete []mvb_edge;

    delete []mvbCommEntity;
}



