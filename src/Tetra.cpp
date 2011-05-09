//
//  Tetra.cpp
//                  2009.06.19
//                  2008.11.27
//                  k.Takeda
#include "Tetra.h"
using namespace pmw;

uint CTetra::mnElemType = ElementType::Tetra;
uint CTetra::mNumOfFace = 4;
uint CTetra::mNumOfEdge = 6;
uint CTetra::mNumOfNode = 4;

// construct & destruct
//
CTetra::CTetra(void)
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
}

CTetra::~CTetra(void)
{
//    //debug
//    cout << "~CTetra" << endl;
}





bool CTetra::IndexCheck(const uint& propType, const uint& index, string& method_name)
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
        pLogger->Info(Utility::LoggerMode::Error, "CTetra::"+ method_name);
        return false;
    }else{
        return true;
    }
}


// out:(PariNode edgeNode)
//
PairNode& CTetra::getPairNode(const uint& edgeIndex)
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
            mEdgeNode.second=mvNode[0];
            break;
        // connect edge (3本)
        case(3):
            mEdgeNode.first=mvNode[0];
            mEdgeNode.second=mvNode[3];
            break;
        case(4):
            mEdgeNode.first=mvNode[1];
            mEdgeNode.second=mvNode[3];
            break;
        case(5):
            mEdgeNode.first=mvNode[2];
            mEdgeNode.second=mvNode[3];
            break;
        default:
            break;
    }

    return mEdgeNode;
}

// out:(vint& pairNodeIndex) -> pair Node index num.
//
void CTetra::getPairNode(vint& pairNodeIndex, const uint& edgeIndex)
{
    switch(edgeIndex){
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
            pairNodeIndex[1]= mvNode[0]->getID();
            break;

        // connect edge (3本)
        case(3):
            pairNodeIndex[0]= mvNode[0]->getID();
            pairNodeIndex[1]= mvNode[3]->getID();
            break;
        case(4):
            pairNodeIndex[0]= mvNode[1]->getID();
            pairNodeIndex[1]= mvNode[3]->getID();
            break;
        case(5):
            pairNodeIndex[0]= mvNode[2]->getID();
            pairNodeIndex[1]= mvNode[3]->getID();
            break;
        default:
            break;
    }
}

// (局所ノード番号、局所ノード番号)に対応した、辺(Edge)番号
//
uint& CTetra::getEdgeIndex(CNode* pNode0, CNode* pNode1)
{
    getLocalNodeNum(pNode0, pNode1);//mvPairNodeLocalNumに値が入る.

    CEdgeTree *edgeTree = CEdgeTree::Instance();

    ////debug
    //cout << "CTetra::getEdgeIndex =" << edgeTree->getTetraEdgeIndex(mvPairNodeLocalNum[0], mvPairNodeLocalNum[1]) << endl;

    return edgeTree->getTetraEdgeIndex(mvPairNodeLocalNum[0], mvPairNodeLocalNum[1]);//Tetra Edge Tree
}

// EdgeElement集合がセット済みか?
//
bool CTetra::isEdgeElem(CNode* pNode0, CNode* pNode1)
{
    //    //
    //    //-- getLocalNodeNum(pNode0, pNode1)と同じ --
    //    //
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
    uint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);//Tetra edge tree

    //debug
    //cout << "CTetra::isEdgeElem edgeNum == " << edgeNum << endl;

    return mvb_edge[edgeNum];
}

// EdgeElement集合がセットされていることをスタンプ
//
void CTetra::setBoolEdgeElem(CNode* pNode0, CNode* pNode1)
{
    //vuint vLocalNum = getLocalNodeNum(pNode0, pNode1);

    // Edgeに要素集合が作られているか? のbool値を返す
    uint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);//Tetra edge tree

    mvb_edge[edgeNum]=true;// スタンプ
}


// 局所ノード番号から、局所Face番号へ変換
//
uint& CTetra::getLocalFaceNum(const vuint& vLocalNodeNum)
{
    CFaceTree *pFaceTree = CFaceTree::Instance();
    
    return pFaceTree->getTetraFaceIndex2(vLocalNodeNum);
}

// 3個のノードから面番号を取得
//
uint& CTetra::getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    // 面構成ノード -> 局所ノード番号
    // 局所ノード番号 -> 面番号==iface
    vuint vLocalNum; vLocalNum.reserve(3);

    vLocalNum[0]= mmIDLocal[pNode0->getID()];
    vLocalNum[1]= mmIDLocal[pNode1->getID()];
    vLocalNum[2]= mmIDLocal[pNode2->getID()];

    CFaceTree *faceTree = CFaceTree::Instance();

    return faceTree->getTetraFaceIndex2(vLocalNum);
}


// 2 Edge -> Face Index
//
uint& CTetra::getFaceIndex(const uint& edge0, const uint& edge1)
{
    CEdgeFaceTree *pEdgeFace = CEdgeFaceTree::Instance();

    return pEdgeFace->getTetraFaceIndex(edge0, edge1);
}

// mvvFaceCnvNodeのセットアップ
//
void CTetra::setupFaceCnvNodes()
{
    // Face構成のノード・コネクティビティ
    //
    CFaceTree *faceTree = CFaceTree::Instance();
    uint* faceConnectivity;

    CNode *pNode;
    mvvFaceCnvNode.resize(4);
    uint iface, ivert, index;

    for(iface=0; iface< 4; iface++){
        mvvFaceCnvNode[iface].resize(3);
        faceConnectivity= faceTree->getLocalNodeTetraFace(iface);

        for(ivert=0; ivert< 3; ivert++){
            index= faceConnectivity[ivert];
            pNode= mvNode[index];

            mvvFaceCnvNode[iface][ivert]=pNode;
        };
    };
}


// Node周囲の接続先Nodeを配列で返す.
// (係数行列 作成用途, CMesh::setupAggregate)
//
vector<CNode*> CTetra::getConnectNode(CNode* pNode)
{
    vector<CNode*> vNode;
    vNode.reserve(3);

    uint gID= pNode->getID();
    uint localID= mmIDLocal[gID];

    CNodeConnectNodeTree *pConnTree = CNodeConnectNodeTree::Instance();
    vuint vLocalID;
    vLocalID= pConnTree->getTetraConnectNode(localID);

    uint i;
    for(i=0; i< 3; i++){
        localID= vLocalID[i];
        vNode.push_back(mvNode[localID]);
    };

    return vNode;
}






