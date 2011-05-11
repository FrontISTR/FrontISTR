
#include <vector>

#include "FaceTree.h"

//
// Quad.cpp
//                      2009.06.19
//                      2008.12.01
//                      k.Takeda
#include "Quad.h"
using namespace pmw;

uiint CQuad::mnElemType = ElementType::Quad;
uiint CQuad::mnElemOrder = 1;
uiint CQuad::mNumOfFace = 1;
uiint CQuad::mNumOfEdge = 4;
uiint CQuad::mNumOfNode = 4;
uiint CQuad::mNumOfVert = 4;
//
//
CQuad::CQuad()
{
    ;
}

CQuad::~CQuad()
{
//    //debug
//    cout << "~CQuad" << endl;
}

void CQuad::initialize()
{
    mvNode.resize(mNumOfNode);

    mvvEdgeElement.resize(mNumOfEdge);

    mvEdgeInterNode.resize(mNumOfEdge);

    mvb_edge = new bool[mNumOfEdge];
    uiint i;
    for(i=0; i< mNumOfEdge; i++){
      mvb_edge[i] = false;
    };

    mvFaceElement.resize(mNumOfFace);
    mvFaceNode.resize(mNumOfFace);

    mvb_face = new bool[mNumOfFace];
    // Face ノード bool
    //
    for(i=0; i< mNumOfFace; i++){
        mvb_face[i] = false;
    };

    //CommElementのprolongation用
    mvProgElement.resize(mNumOfNode);


    //MPC属性
    mvbMPCFace = new bool[mNumOfFace];
    for(i=0; i< mNumOfFace; i++){
        mvbMPCFace[i]= false;
    };

    //通信界面(CommMesh2)属性:辺
    mvbCommEntity = new bool[mNumOfEdge];
    for(i=0; i< mNumOfEdge; i++){
        mvbCommEntity[i]= false;
    };
}

const uiint& CQuad::getType()
{
    return mnElemType;
}
const uiint& CQuad::getOrder()
{
    return mnElemOrder;
}


// method
// --
//
bool CQuad::IndexCheck(const uiint& propType, const uiint& index, string& method_name)
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
        pLogger->Info(Utility::LoggerMode::Error, "CQuad::"+ method_name);
        return false;
    }else{
        return true;
    }
}



// out:(PariNode edgeNode)
//
PairNode CQuad::getPairNode(const uiint& iedge)
{
    PairNode pairNode;

    switch(iedge){
        //surf 0
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
        default:
            break;
    }

    return pairNode;
}

// out:(vint& pairNodeIndex) -> pair Node index num.
//
void CQuad::getPairNode(vint& pairNodeIndex, const uiint& iedge)
{
    switch(iedge){
        //surf 0 (low surf)
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

// (局所ノード番号、局所ノード番号)に対応した辺(Edge)番号
//
uiint& CQuad::getEdgeIndex(CNode* pNode0, CNode* pNode1)
{
    uiint id0, id1;//MeshでのノードのIndex番号
    id0 = pNode0->getID();
    id1 = pNode1->getID();
    
    
    CEdgeTree *edgeTree = CEdgeTree::Instance();

    return edgeTree->getQuadEdgeIndex(mmIDLocal[id0], mmIDLocal[id1]);//Quad EdgeTree
}
uiint& CQuad::getEdgeIndex(const uiint& nodeID_0, const uiint& nodeID_1)
{
    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    return pEdgeTree->getQuadEdgeIndex(mmIDLocal[nodeID_0], mmIDLocal[nodeID_1]);
}


// EdgeElement集合がセット済みか?
//
bool CQuad::isEdgeElem(CNode* pNode0, CNode* pNode1)
{
    //    //------------------------------------------
    //    //   getLocalNodeNum(pNode0, pNode1)と同じ
    //    //------------------------------------------
    //    uint index0, index1;//MeshでのノードのIndex番号
    //    uint local0, local1;//局所ノード番号
    //
    //    index0= pNode0->getID();  index1= pNode1->getID();
    //
    //    // pNode0とpNode1のグローバルなIndex番号 => 局所番号に変換
    //    //
    //    uint i;
    //    for(i=0; i< mvNode.size(); i++){
    //        if(index0 == mvNode[i]->getID()) local0 = i;
    //        if(index1 == mvNode[i]->getID()) local1 = i;
    //    };
    //    //--------------------------------------------


    // Edgeに要素集合が作られているか? のbool値を返す
    uiint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);//Quad

    return mvb_edge[edgeNum];
}

// EdgeElement集合がセットされていることをスタンプ
//
void CQuad::setBoolEdgeElem(CNode* pNode0, CNode* pNode1)
{
    //vuint vLocalNum = getLocalNodeNum(pNode0, pNode1);

    // Edgeに要素集合が作られているか? のbool値を返す
    uiint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);//Quad edge tree

    mvb_edge[edgeNum]=true;// スタンプ
}


// 局所ノード番号から、局所Face番号へ変換
//
uiint& CQuad::getLocalFaceNum(const vuint& vLocalNodeNum)
{
    CFaceTree *pFaceTree = CFaceTree::Instance();
    
    uiint faceIndex = pFaceTree->getQuadFaceIndex(vLocalNodeNum);

    if(faceIndex==0){
        mTempo=0;
    }else{
        Utility::CLogger *pLogger = Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"mismatch!, Local Node number@CQuad::getLocalFaceNum");
        mTempo=1;
    }
    return mTempo;
}

//
//
uiint& CQuad::getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    // 面構成ノード -> 局所ノード番号
    // 局所ノード番号 -> 面番号==iface
    vuint vLocalNum; vLocalNum.resize(3);

    vLocalNum[0]= mmIDLocal[pNode0->getID()];
    vLocalNum[1]= mmIDLocal[pNode1->getID()];
    vLocalNum[2]= mmIDLocal[pNode2->getID()];

    CFaceTree *faceTree = CFaceTree::Instance();

    uiint faceIndex = faceTree->getQuadFaceIndex(vLocalNum);
    
    if(faceIndex==0){
        mTempo=0;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"mismatch, Local Node number @CQuad::getLocalFaceNum(arg:pNode)");
        mTempo=1;
    }
    return mTempo;
}


// 2 Edge => Face Index
//
uiint& CQuad::getFaceIndex(const uiint& edge0, const uiint& edge1)
{
    CEdgeFaceTree *pTree = CEdgeFaceTree::Instance();

    return pTree->getQuadFaceIndex(edge0, edge1);
}

// Face構成Node
//
vector<CNode*> CQuad::getFaceCnvNodes(const uiint& iface)
{
    // Face構成のノード・コネクティビティ
    //
    CNode *pNode;
    vector<CNode*> vFaceCnvNode;
    vFaceCnvNode.resize(4);
    uiint ivert;

    for(ivert=0; ivert< 4; ivert++){
        pNode= mvNode[ivert];
        vFaceCnvNode[ivert]=pNode;
    };

    return vFaceCnvNode;
}


// Node周囲の接続先Nodeを配列で返す.
// (係数行列 作成用途, CMesh::setupAggregate)
//
vector<CNode*> CQuad::getConnectNode(CNode* pNode)
{
    vector<CNode*> vNode;
    vNode.reserve(2);

    uiint gID= pNode->getID();
    uiint localID= mmIDLocal[gID];

    CNodeConnectNodeTree *pConnTree = CNodeConnectNodeTree::Instance();
    vuint vLocalID;
    vLocalID= pConnTree->getQuadConnectNode(localID);

    uiint i;
    for(i=0; i< 2; i++){
        localID= vLocalID[i];
        vNode.push_back(mvNode[localID]);
    };

    return vNode;
}


//
// 1. 辺-面 Element*配列を解放
// 2. 辺-面 Node*   配列を解放 (2次要素は辺ノードを残す)
//
void CQuad::deleteProgData()
{
    // Edge
    uiint iedge;
    for(iedge=0; iedge < mNumOfEdge; iedge++){
        vector<CElement*>().swap(mvvEdgeElement[iedge]);
    };
    vector<vector<CElement*> >().swap(mvvEdgeElement);
    vector<CNode*>().swap(mvEdgeInterNode);//2次要素の場合は残すこと
    delete []mvb_edge;


    // Face
    vector<CElement*>().swap(mvFaceElement);
    vector<CNode*>().swap(mvFaceNode);
    delete []mvb_face;


    delete []mvbMPCFace;   // 面番号のどれがMPCするのか
    delete []mvbCommEntity;// CommMesh2の属性
}








