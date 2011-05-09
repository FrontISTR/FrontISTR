
#include <vector>

//
// Hexa.cpp
//                  2009.6.19
//                  2009.4.20
//                  k.Takeda
#include "Hexa.h"
using namespace pmw;

uint CHexa::mnElemType = ElementType::Hexa;
uint CHexa::mNumOfFace =  6;
uint CHexa::mNumOfEdge = 12;
uint CHexa::mNumOfNode =  8;


#include "Logger.h"
//
//
CHexa::CHexa(void):CSolidBase()
{
    mvNode.resize(mNumOfNode);
    mvvEdgeElement.resize(mNumOfEdge);
    mvEdgeInterNode.resize(mNumOfEdge);

    // Edgeノード bool
    //
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

CHexa::~CHexa(void)
{
//    //debug
//    cout << "~CHexa" << endl;
}

//// each reserve for EdgeElement
////
//void CHexa::reserveEdgeElement(const uint& edgeIndex, const uint& numOfElem)
//{
//    Utility::CLogger *pLogger = Utility::CLogger::Instance();
//
//    if(edgeIndex >= mNumOfEdge){
//        pLogger->Info(Utility::LoggerMode::Error, "invalid EdgeIndex @CHexa::reserveEdgeElement");
//        return;
//    }
//
//    mvvEdgeElement[edgeIndex].reserve(numOfElem);
//}
//
//// set to EdgeElement
////
//void CHexa::setEdgeElement(const uint& edgeIndex, CElement* pElem)
//{
//    Utility::CLogger *pLogger = Utility::CLogger::Instance();
//
//    if(edgeIndex >= mNumOfEdge){
//        pLogger->Info(Utility::LoggerMode::Error, "invalid EdgeIndex @CHexa::setEdgeElement");
//        return;
//    }
//
//    mvvEdgeElement[edgeIndex].push_back(pElem);
//}

bool CHexa::IndexCheck(const uint& propType, const uint& index, string& method_name)
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
        pLogger->Info(Utility::LoggerMode::Error, "CHexa::"+ method_name);

        return false;
    }else{
        return true;
    }
}


// out:(PariNode edgeNode)
//
PairNode& CHexa::getPairNode(const uint& edgeIndex)
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
        //surf 1 (upper surf)
        case(4):
            mEdgeNode.first=mvNode[4];
            mEdgeNode.second=mvNode[5];
            break;
        case(5):
            mEdgeNode.first=mvNode[5];
            mEdgeNode.second=mvNode[6];
            break;
        case(6):
            mEdgeNode.first=mvNode[6];
            mEdgeNode.second=mvNode[7];
            break;
        case(7):
            mEdgeNode.first=mvNode[7];
            mEdgeNode.second=mvNode[4];
            break;
        // connect edge (4本)
        case(8):
            mEdgeNode.first=mvNode[0];
            mEdgeNode.second=mvNode[4];
            break;
        case(9):
            mEdgeNode.first=mvNode[1];
            mEdgeNode.second=mvNode[5];
            break;
        case(10):
            mEdgeNode.first=mvNode[2];
            mEdgeNode.second=mvNode[6];
            break;
        case(11):
            mEdgeNode.first=mvNode[3];
            mEdgeNode.second=mvNode[7];
            break;
        default:
            break;
    }

    return mEdgeNode;// PairNode
}

// out:(vint& pairNodeIndex) -> pair Node index num.
//
void CHexa::getPairNode(vint& pairNodeIndex, const uint& edgeIndex)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    if(pairNodeIndex.size() != 2){
        pLogger->Info(Utility::LoggerMode::Error,"pair node size Error@CHexa::getPairNode");
        return;
    }
    
    // output is pairNodeIndex[2]
    //
    switch(edgeIndex){
        //surf 0 (low surf)
        case(0):
            pairNodeIndex[0]= mvNode[0]->getID();// <= getIndex()も用意したが、IDに統一する
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
        //surf 1 (upper surf)
        case(4):
            pairNodeIndex[0]= mvNode[4]->getID();
            pairNodeIndex[1]= mvNode[5]->getID();
            break;
        case(5):
            pairNodeIndex[0]= mvNode[5]->getID();
            pairNodeIndex[1]= mvNode[6]->getID();
            break;
        case(6):
            pairNodeIndex[0]= mvNode[6]->getID();
            pairNodeIndex[1]= mvNode[7]->getID();
            break;
        case(7):
            pairNodeIndex[0]= mvNode[7]->getID();
            pairNodeIndex[1]= mvNode[4]->getID();
            break;
        // connect edge (4本)
        case(8):
            pairNodeIndex[0]= mvNode[0]->getID();
            pairNodeIndex[1]= mvNode[4]->getID();
            break;
        case(9):
            pairNodeIndex[0]= mvNode[1]->getID();
            pairNodeIndex[1]= mvNode[5]->getID();
            break;
        case(10):
            pairNodeIndex[0]= mvNode[2]->getID();
            pairNodeIndex[1]= mvNode[6]->getID();
            break;
        case(11):
            pairNodeIndex[0]= mvNode[3]->getID();
            pairNodeIndex[1]= mvNode[7]->getID();
            break;
        default:
            break;
    }
}


// (局所ノード番号、局所ノード番号)に対応した、辺(Edge)のIndex番号
//
uint& CHexa::getEdgeIndex(CNode* pNode0, CNode* pNode1)
{
    getLocalNodeNum(pNode0, pNode1);//mvPairNodeLocalNumに値が入る.
    
    CEdgeTree *edgeTree = CEdgeTree::Instance();
    
    ////debug
    //Utility::CLogger *pLogger= Utility::CLogger::Instance();
    //pLogger->Info(Utility::LoggerMode::MWDebug,"CHexa::getEdgeIndex,  edgeIndex => ",edgeTree->getHexaEdgeIndex(mvPairNodeLocalNum[0],mvPairNodeLocalNum[1]));

    return edgeTree->getHexaEdgeIndex(mvPairNodeLocalNum[0], mvPairNodeLocalNum[1]);//Hexa 辺ツリー
}
uint& CHexa::getEdgeIndex(const uint& nodeID_0, const uint& nodeID_1)
{
    mvPairNodeLocalNum[0]= mmIDLocal[nodeID_0];
    mvPairNodeLocalNum[1]= mmIDLocal[nodeID_1];

    CEdgeTree *pEdgeTree= CEdgeTree::Instance();

    return pEdgeTree->getHexaEdgeIndex(mvPairNodeLocalNum[0], mvPairNodeLocalNum[1]);
}

// EdgeElement集合がセット済みか?
//
bool CHexa::isEdgeElem(CNode* pNode0, CNode* pNode1)
{
    // Edgeに要素集合が作られているか? のbool値を返す
    uint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);//Hexa 辺ツリー

    return mvb_edge[edgeNum];
}

// EdgeElement集合がセットされていることをスタンプ
//
void CHexa::setBoolEdgeElem(CNode* pNode0, CNode* pNode1)
{
    // Edgeに要素集合が作られているか? のbool値を返す
    uint edgeNum;
    edgeNum = getEdgeIndex(pNode0, pNode1);//Hexa edgeツリー

    mvb_edge[edgeNum]=true;// スタンプ
}



// 局所ノード番号から、局所Face番号へ変換
//
uint& CHexa::getLocalFaceNum(const vuint& vLocalNodeNum)
{
    CFaceTree *faceTree = CFaceTree::Instance();

    return faceTree->getHexaFaceIndex2(vLocalNodeNum);
}



//// 局所面番号に対応する、ノード配列を提供.
////
//vector<CNode*>& CHexa::getFaceNode(const uint& faceIndex)
//{
//    return mvvFaceNode[faceIndex];
//}

//// 局所面番号に対応する、ノード・インデックス配列を提供.
////
//vuint& CHexa::getFaceNodeID(const uint& faceIndex)
//{
//    vuint vNodeIndex;
//
//    return vNodeIndex;
//}

// 2 Edge -> Face Index( 2本の辺から、面番号を取得 )
// 
uint& CHexa::getFaceIndex(const uint& edge0, const uint& edge1)
{
    CEdgeFaceTree *faceTree= CEdgeFaceTree::Instance();

    return faceTree->getHexaFaceIndex(edge0, edge1);
}

uint& CHexa::getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    // 面構成ノード -> 局所ノード番号
    // 局所ノード番号 -> 面番号==iface
    vuint vLocalNum; vLocalNum.reserve(3);

    vLocalNum[0]= mmIDLocal[pNode0->getID()];
    vLocalNum[1]= mmIDLocal[pNode1->getID()];
    vLocalNum[2]= mmIDLocal[pNode2->getID()];

    CFaceTree *faceTree = CFaceTree::Instance();

    return faceTree->getHexaFaceIndex2(vLocalNum);
}


// mvvFaceCnvNodeのセットアップ
//
void CHexa::setupFaceCnvNodes()
{
    // Face構成のノード・コネクティビティ
    //
    CFaceTree *faceTree = CFaceTree::Instance();
    uint* faceConnectivity;

    CNode *pNode;
    mvvFaceCnvNode.resize(6);
    uint iface, ivert, index;

    for(iface=0; iface< 6; iface++){
        mvvFaceCnvNode[iface].resize(4);
        faceConnectivity= faceTree->getLocalNodeHexaFace(iface);

        for(ivert=0; ivert< 4; ivert++){
            index= faceConnectivity[ivert];
            pNode= mvNode[index];

            mvvFaceCnvNode[iface][ivert]=pNode;
        };
    };
}


// Node周囲の接続先Nodeを配列で返す.
// (係数行列 作成用途, CMesh::setupAggregate)
//
vector<CNode*> CHexa::getConnectNode(CNode* pNode)
{
    vector<CNode*> vNode;
    vNode.reserve(3);

    uint gID= pNode->getID();
    uint localID= mmIDLocal[gID];

    CNodeConnectNodeTree *pConnTree = CNodeConnectNodeTree::Instance();
    vuint vLocalID;
    vLocalID= pConnTree->getHexaConnectNode(localID);

    uint i;
    for(i=0; i< 3; i++){
        localID= vLocalID[i];
        vNode.push_back(mvNode[localID]);
    };

    return vNode;
}










