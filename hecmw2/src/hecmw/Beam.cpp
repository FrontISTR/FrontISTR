//
//  Beam.cpp
//			2009.06.19
//			2008.12.01
//			k.Takeda
#include <vector>

#include "Beam.h"
using namespace pmw;

uiint CBeam::mnElemType = ElementType::Beam;
uiint CBeam::mnElemOrder = 1;
uiint CBeam::mNumOfFace = 0;
uiint CBeam::mNumOfEdge = 1;
uiint CBeam::mNumOfNode = 2;
uiint CBeam::mNumOfVert = 2;
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
    uiint i;
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
PairNode CBeam::getPairNode(const uiint& iedge)
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
void CBeam::getPairNode(vint& pairNodeIndex, const uiint& iedge)
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
uiint& CBeam::getEdgeIndex(CNode* pNode0, CNode* pNode1)
{
    uiint id0, id1;//MeshでのノードのIndex番号
    id0 = pNode0->getID();
    id1 = pNode1->getID();
    
    // Beam要素なので,”Edge_Index=0”
    CEdgeTree *edgeTree = CEdgeTree::Instance();
    
    return edgeTree->getBeamEdgeIndex(mmIDLocal[id0], mmIDLocal[id1]);
}
uiint& CBeam::getEdgeIndex(const uiint& nodeID_0, const uiint& nodeID_1)
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
uiint& CBeam::getLocalFaceNum(const vuint& vLocalNodeNum)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::getLocalFaceNum");

    return pLogger->getUDummyValue();
}
//
//
uiint& CBeam::getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method @CBeam::getLocalFaceNum(arg:pNode)");

    return pLogger->getUDummyValue();
}


// 2 Edge => Face Index
//
uiint& CBeam::getFaceIndex(const uiint& edge0, const uiint& edge1)
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
vector<CNode*> CBeam::getFaceCnvNodes(const uiint& iface)
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

    uiint gID= pNode->getID();
    uiint localID= mmIDLocal[gID];

    CNodeConnectNodeTree *pConnTree = CNodeConnectNodeTree::Instance();
    vuint vLocalID;
    vLocalID= pConnTree->getBeamConnectNode(localID);

    uiint i;
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
    uiint iedge;
    for(iedge=0; iedge < mNumOfEdge; iedge++){
        vector<CElement*>().swap(mvvEdgeElement[iedge]);
    };
    vector<vector<CElement*> >().swap(mvvEdgeElement);
    vector<CNode*>().swap(mvEdgeInterNode);//2次要素の場合は残すこと
    delete []mvb_edge;

    delete []mvbCommEntity;
}



