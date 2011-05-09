
#include <vector>

//
//  Element.cpp
//
//			2008.12.12
//			2008.11.5
//			k.Takeda
#include "Element.h"
using namespace pmw;

// construct & destruct
//
CElement::CElement(void)
{
    mvPairNodeLocalNum.resize(2);//局所ノード番号ペアとして使用.

    mbComm =false;//デフォルトは,CommElementではない.
    mbDComm=false;//デフォルトは計算領域 => DCommElementではない.
    mbRComm=false;//デフォルトは,RCommElementではない(CommMeshに含まれているが,Nodeが全てmyRank)

    mCommID= -1;


    //MPC属性
    mbMaster=false;
    mbSlave= false;
}

CElement::~CElement(void)
{
    ////debug
    //cout << "~CElement" << endl;
}

// ノードのセット、map<>へのセット(key:ID, val:ローカル番号)
//
void CElement::setNode(CNode* pNode, const uint& local_id)
{
    mvNode[local_id] = pNode;
    
    mmIDLocal[pNode->getID()]= local_id;
}

// each reserve for EdgeElement
//
void CElement::reserveEdgeElement(const uint& edgeIndex, const uint& numOfElem)
{
    string s_msg("reserveEdgeElement");
    
    if(IndexCheck(ElementPropType::Edge, edgeIndex, s_msg)){
        mvvEdgeElement[edgeIndex].reserve(numOfElem);
    }else{
        return;
    }

    
}

// set to EdgeElement
//
void CElement::setEdgeElement(const uint& edgeIndex, CElement* pElem)
{
    string s_msg("setEdgeElement");

    if(IndexCheck(ElementPropType::Edge, edgeIndex, s_msg)){
        mvvEdgeElement[edgeIndex].push_back(pElem);
    }else{
        return;
    }
}

// 既に集められた要素集合をセットする
//
void CElement::setEdgeAggElement(const uint& edgeIndex, vector<CElement*> vElement)
{
    mvvEdgeElement[edgeIndex] = vElement;
}

// set to InterNode on Edge
//            prolongation用のノード
void CElement::setEdgeInterNode(CNode* pNode, const uint& edgeIndex)
{
    mvEdgeInterNode[edgeIndex]=pNode;
}

// NodeのグローバルIndex番号から、要素内の局所番号へ変換
//
vuint& CElement::getLocalNodeNum(CNode* pNode0, CNode* pNode1)
{
    uint id0, id1;//MeshでのノードのIndex番号
    id0 = pNode0->getID();
    id1 = pNode1->getID();
    
    mvPairNodeLocalNum[0]= mmIDLocal[id0];
    mvPairNodeLocalNum[1]= mmIDLocal[id1];

    return mvPairNodeLocalNum;
}


// 面(Face)に既に面ノードがセット済みかどうか？
//
bool CElement::isFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    uint iface= getFaceIndex(pNode0,pNode1,pNode2);

    //debug
    //cout << "isFaceElem, iface => " << iface << endl;

    return mvb_face[iface];
}


// 面(Face)に面ノードがセットされたことをスタンプ
//
void CElement::setBoolFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    uint iface = getFaceIndex(pNode0,pNode1,pNode2);

    ////debug
    //cout << "iface = " << iface << ", CElement::setBoolFaceElem" << endl;

    mvb_face[iface]=true;
}


// 子要素をNodeIDから提供する -> CSkinFace::refine()で利用
// 
CElement* CElement::getProgElem_NodeID(const uint& nodeID)
{
    uint ivert= mmIDLocal[nodeID];//NodeID -> 局所番号に変換
    
    return mvProgElement[ivert];
}









