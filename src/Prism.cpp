
#include "EdgeFaceTree.h"

//
//  Prism.cpp
//
//  プリズム(三角柱)
//
//                          2009.06.29
//                          2009.06.29
//                          k.Takeda
#include "Prism.h"
using namespace pmw;

uint CPrism::mnElemType = ElementType::Prism;
uint CPrism::mNumOfFace = 5;
uint CPrism::mNumOfEdge = 9;
uint CPrism::mNumOfNode = 6;

// コンストラクター&
//          デストラクター
//
CPrism::CPrism()
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
    // Face ノード bool
    //
    for(i=0; i< mNumOfFace; i++){
        mvb_face.push_back(false);
    };

    //CommElementのprolongation用
    mvProgElement.resize(mNumOfNode);

    //MPC属性
    mvbMPCFace.resize(mNumOfFace);
    for(i=0; i< mNumOfFace; i++){
        mvbMPCFace[i]= false;
    };
}

CPrism::~CPrism()
{
//    //debug
//    cout << "~CPrism" << endl;
}


// method
// --
//
bool CPrism::IndexCheck(const uint& propType, const uint& index, string& method_name)
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
        pLogger->Info(Utility::LoggerMode::Error, "CPrism::"+ method_name);
        return false;
    }else{
        return true;
    }
}


// out:(PariNode edgeNode)
//
PairNode& CPrism::getPairNode(const uint& edgeIndex)
{
    switch(edgeIndex){
        //surf 0 (low surf)
        case(0):
            mEdgeNode.first=mvNode[0];
            mEdgeNode.second=mvNode[1];
            break;
        case(1):
            mEdgeNode.first=mvNode[0];
            mEdgeNode.second=mvNode[2];
            break;
        case(2):
            mEdgeNode.first=mvNode[1];
            mEdgeNode.second=mvNode[2];
            break;
        // 3本の柱
        case(3):
            mEdgeNode.first=mvNode[0];
            mEdgeNode.second=mvNode[3];
            break;
        case(4):
            mEdgeNode.first=mvNode[1];
            mEdgeNode.second=mvNode[4];
            break;
        case(5):
            mEdgeNode.first=mvNode[2];
            mEdgeNode.second=mvNode[5];
            break;
        // surf 1 (upper surf)
        case(6):
            mEdgeNode.first=mvNode[3];
            mEdgeNode.second=mvNode[4];
            break;
        case(7):
            mEdgeNode.first=mvNode[4];
            mEdgeNode.second=mvNode[5];
            break;
        case(8):
            mEdgeNode.first=mvNode[5];
            mEdgeNode.second=mvNode[3];
            break;
        default:
            break;
    }

    return mEdgeNode;
}

// out:(vint& pairNodeIndex) -> pair Node index num.
//
void CPrism::getPairNode(vint& pairNodeIndex, const uint& edgeIndex)
{
    switch(edgeIndex){
        //surf 0 (low surf)
        case(0):
            pairNodeIndex[0]= mvNode[0]->getID();
            pairNodeIndex[1]= mvNode[1]->getID();
            break;
        case(1):
            pairNodeIndex[0]=mvNode[0]->getID();
            pairNodeIndex[1]=mvNode[2]->getID();
            break;
        case(2):
            pairNodeIndex[0]=mvNode[1]->getID();
            pairNodeIndex[1]=mvNode[2]->getID();
            break;
        // 3本の柱
        case(3):
            pairNodeIndex[0]=mvNode[0]->getID();
            pairNodeIndex[1]=mvNode[3]->getID();
            break;
        case(4):
            pairNodeIndex[0]=mvNode[1]->getID();
            pairNodeIndex[1]=mvNode[4]->getID();
            break;
        case(5):
            pairNodeIndex[0]=mvNode[2]->getID();
            pairNodeIndex[1]=mvNode[5]->getID();
            break;
        // surf 1 (upper surf)
        case(6):
            pairNodeIndex[0]=mvNode[3]->getID();
            pairNodeIndex[1]=mvNode[4]->getID();
            break;
        case(7):
            pairNodeIndex[0]=mvNode[4]->getID();
            pairNodeIndex[1]=mvNode[5]->getID();
            break;
        case(8):
            pairNodeIndex[0]=mvNode[5]->getID();
            pairNodeIndex[1]=mvNode[3]->getID();
            break;
        default:
            break;
    }
}


// (局所ノード番号、局所ノード番号)に対応した、辺(Edge)番号
//
uint& CPrism::getEdgeIndex(CNode* pNode0, CNode* pNode1)
{
    getLocalNodeNum(pNode0, pNode1);//mvPairNodeLocalNumに値が入る.

    CEdgeTree *edgeTree = CEdgeTree::Instance();

    return edgeTree->getPrismEdgeIndex(mvPairNodeLocalNum[0], mvPairNodeLocalNum[1]);//Prism Edge Tree
}
uint& CPrism::getEdgeIndex(const uint& nodeID_0, const uint& nodeID_1)
{
    mvPairNodeLocalNum[0]= mmIDLocal[nodeID_0];
    mvPairNodeLocalNum[1]= mmIDLocal[nodeID_1];

    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    return pEdgeTree->getPrismEdgeIndex(mvPairNodeLocalNum[0], mvPairNodeLocalNum[1]);
}

// EdgeElement集合がセット済みか?
//
bool CPrism::isEdgeElem(CNode* pNode0, CNode* pNode1)
{
    //
    //-- getLocalNodeNum(pNode0, pNode1)と同じ --
    //
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
    //
    //vuint vlocal = getLocalNodeNum(pNode0, pNode1);


    // Edgeに要素集合が作られているか? のbool値を返す
    uint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);//Prism edge tree

    return mvb_edge[edgeNum];
}

// EdgeElement集合がセットされていることをスタンプ
//
void CPrism::setBoolEdgeElem(CNode* pNode0, CNode* pNode1)
{
    //vuint vLocalNum = getLocalNodeNum(pNode0, pNode1);

    // Edgeに要素集合が作られているか? のbool値を返す
    uint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);//Tetra edge tree

    mvb_edge[edgeNum]=true;// スタンプ
}


// 局所ノード番号から、局所Face番号へ変換
//
uint& CPrism::getLocalFaceNum(const vuint& vLocalNodeNum)
{
    CFaceTree *pFaceTree = CFaceTree::Instance();
    
    return pFaceTree->getPrismFaceIndex2(vLocalNodeNum);
}

// 3個のノードから面番号を取得
//
uint& CPrism::getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    // 面構成ノード -> 局所ノード番号
    // 局所ノード番号 -> 面番号==iface
    vuint vLocalNum; vLocalNum.reserve(3);

    vLocalNum[0]= mmIDLocal[pNode0->getID()];
    vLocalNum[1]= mmIDLocal[pNode1->getID()];
    vLocalNum[2]= mmIDLocal[pNode2->getID()];

    CFaceTree *faceTree = CFaceTree::Instance();

    return faceTree->getPrismFaceIndex2(vLocalNum);
}


// 2 Edge -> Face Index
//
uint& CPrism::getFaceIndex(const uint& edge0, const uint& edge1)
{
    CEdgeFaceTree *pTree = CEdgeFaceTree::Instance();

    return pTree->getPrismFaceIndex(edge0, edge1);
}

// mvvFaceCnvNodeのセットアップ
//
void CPrism::setupFaceCnvNodes()
{
    // Face構成のノード・コネクティビティ
    //
    CFaceTree *faceTree = CFaceTree::Instance();
    uint* faceConnectivity;

    CNode *pNode;
    mvvFaceCnvNode.resize(5);
    uint iface, ivert, index;
    uint numOfVert;

    for(iface=0; iface< 5; iface++){
        if(iface==0 || iface==1){
            numOfVert=3;
        }else{
            numOfVert=4;
        }
        mvvFaceCnvNode[iface].resize(numOfVert);
        faceConnectivity= faceTree->getLocalNodePrismFace(iface);

        for(ivert=0; ivert< numOfVert; ivert++){
            index= faceConnectivity[ivert];
            pNode= mvNode[index];

            mvvFaceCnvNode[iface][ivert]=pNode;
        };
    };
}


// Node周囲の接続先Nodeを配列で返す.
// (係数行列 作成用途, CMesh::setupAggregate)
//
vector<CNode*> CPrism::getConnectNode(CNode* pNode)
{
    vector<CNode*> vNode;
    vNode.reserve(3);

    uint gID= pNode->getID();
    uint localID= mmIDLocal[gID];

    CNodeConnectNodeTree *pConnTree = CNodeConnectNodeTree::Instance();
    vuint vLocalID;
    vLocalID= pConnTree->getPrismConnectNode(localID);

    uint i;
    for(i=0; i< 3; i++){
        localID= vLocalID[i];
        vNode.push_back(mvNode[localID]);
    };

    return vNode;
}












