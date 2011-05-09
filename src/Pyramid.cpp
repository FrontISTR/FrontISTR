//
//  Pyramid.cpp
//  ピラミッド(四角錐)
//
//
//                          2009.06.29
//                          2009.06.29
//                          k.Takeda

#include "Pyramid.h"
using namespace pmw;

uint CPyramid::mnElemType = ElementType::Pyramid;
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
    // Edge ノード bool
    //
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
    mvProgElement.resize(8);//Hexa==4個,Pyramid==4個


    //MPC属性
    mvbMPCFace.resize(mNumOfFace);
    for(i=0; i< mNumOfFace; i++){
        mvbMPCFace[i]= false;
    };
}

CPyramid::~CPyramid()
{
//    //debug
//    cout << "~CPyramid" << endl;
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


// out:(PariNode edgeNode)
//
PairNode& CPyramid::getPairNode(const uint& edgeIndex)
{
    switch(edgeIndex){
        //surf 0 (low surf)
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
        // 4本の柱
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

// out:(vint& pairNodeIndex) -> pair Node index num.
//
void CPyramid::getPairNode(vint& pairNodeIndex, const uint& edgeIndex)
{
    switch(edgeIndex){
        //surf 0 (low surf)
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
        // 4本の柱
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


// (局所ノード番号、局所ノード番号)に対応した、辺(Edge)番号
//
uint& CPyramid::getEdgeIndex(CNode* pNode0, CNode* pNode1)
{
    getLocalNodeNum(pNode0, pNode1);//mvPairNodeLocalNumに値が入る.

    CEdgeTree *edgeTree = CEdgeTree::Instance();

    return edgeTree->getPyramidEdgeIndex(mvPairNodeLocalNum[0], mvPairNodeLocalNum[1]);//Pyramid Edge Tree
}
uint& CPyramid::getEdgeIndex(const uint& nodeID_0, const uint& nodeID_1)
{
    mvPairNodeLocalNum[0]= mmIDLocal[nodeID_0];
    mvPairNodeLocalNum[1]= mmIDLocal[nodeID_1];

    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    return pEdgeTree->getPyramidEdgeIndex(mvPairNodeLocalNum[0], mvPairNodeLocalNum[1]);
}


// EdgeElement集合がセット済みか?
//
bool CPyramid::isEdgeElem(CNode* pNode0, CNode* pNode1)
{
    //vuint vlocal = getLocalNodeNum(pNode0, pNode1);

    // Edgeに要素集合が作られているか? のbool値を返す
    uint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);//Pyramid edge tree

    return mvb_edge[edgeNum];
}

// EdgeElement集合がセットされていることをスタンプ
//
void CPyramid::setBoolEdgeElem(CNode* pNode0, CNode* pNode1)
{
    //vuint vLocalNum = getLocalNodeNum(pNode0, pNode1);

    // Edgeに要素集合が作られているか? のbool値を返す
    uint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);//Pyramid edge tree

    mvb_edge[edgeNum]=true;// スタンプ
}


// 局所ノード番号から、局所Face番号へ変換
//
uint& CPyramid::getLocalFaceNum(const vuint& vLocalNodeNum)
{
    CFaceTree *pFaceTree = CFaceTree::Instance();
    
    return  pFaceTree->getPyramidFaceIndex2(vLocalNodeNum);
}

// 3個のノードから面番号を取得
//
uint& CPyramid::getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    // 面構成ノード -> 局所ノード番号
    // 局所ノード番号 -> 面番号==iface
    vuint vLocalNum; vLocalNum.reserve(3);

    vLocalNum[0]= mmIDLocal[pNode0->getID()];
    vLocalNum[1]= mmIDLocal[pNode1->getID()];
    vLocalNum[2]= mmIDLocal[pNode2->getID()];

    CFaceTree *faceTree = CFaceTree::Instance();

    return faceTree->getPyramidFaceIndex2(vLocalNum);
}


// 2 Edge => Face Index
//
uint& CPyramid::getFaceIndex(const uint& edge0, const uint& edge1)
{
    CEdgeFaceTree *pTree= CEdgeFaceTree::Instance();

    return pTree->getPyramidFaceIndex(edge0, edge1);
}

// mvvFaceCnvNodeのセットアップ
//
void CPyramid::setupFaceCnvNodes()
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


// Node周囲の接続先Nodeを配列で返す.
// (係数行列 作成用途, CMesh::setupAggregate)
//
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
    // Pyramidは,頂点のNode接続先のみ4ノードあるので,
    //    size()で接続先の個数を拾っている.
    // --
    for(i=0; i< vLocalID.size(); i++){
        localID= vLocalID[i];
        vNode.push_back(mvNode[localID]);
    };

    return vNode;
}
















