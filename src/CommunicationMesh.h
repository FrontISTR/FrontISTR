///*
// * File:   CommunicationTable.h
// * Author: ktakeda
// *
// *  Meshに備わる,通信領域(CommID)を管理
// *
// *  CommID:通信領域番号
// *  RankID:計算領域番号
// *
// *  -- Memo --
// *  通信インデックス順位： 要素順 -> 要素内の局所番号順
// *  ----------
// *
// * Modify     2009/08/28
// * Created on 2009/06/19, 19:42
// */
#ifndef _COMMMesh_H_39e5ebb0_da94_4efc_8616_8c2c5918cc4f
#define _COMMMesh_H_39e5ebb0_da94_4efc_8616_8c2c5918cc4f

//#include <cstdlib>// div_t 利用
#include "TypeDef.h"//"CommonStd.hを含む

#include "Element.h"
#include "Node.h"

#include "CommElement.h"
#include <map>

#include "IndexBucket.h"

//// Node around Elements (CMeshが所有)
//#include "AggregateElement.h"

#include "Logger.h"

#include "QuickSort.h"//関数template

namespace pmw{
class CCommMesh{
public:
    CCommMesh(CIndexBucket* pBucket);
    virtual ~CCommMesh();
    
protected:
    uint mCommID;//Meshが管理する通信領域番号(CommMesh自身の番号:CommunicationID)
    uint mRankID;//CommMeshが属する,計算領域のID
    
    // 通信領域は,1対1になるように定義されている
    // --
    uint mTransmitRankID;//通信相手の計算領域ID:常に1対1
    
    
    // CommMeshの要素
    // --
    vector<CCommElement*> mvCommElementAll;//CommMesh本体(計算不要領域を含むCommElement全体)
    vector<CCommElement*> mvCommElement; //計算領域間で使用される共有要素
    vector<CCommElement*> mvDCommElement;//Dead CommElement(計算に使われないCommElement)


    // CommMesh全体のノード -> CCommElement::setCommNodeIndex()にてセット
    // --
    vector<CNode*>  mvNode;//CommMeshの全体ノード:prolongation時に,CommNodeIDを生成する際に必要

    // IDからIndex逆引き
    // --
    CIndexBucket *mpBucket;//Meshのバケット


    // Nodeの"rank" 一覧 (Send,RecvNode収集に使用)
    // --
    vuint mvNodeRank;


    // mvNodeから送受信ノードに振り分けた配列 -> CCommElement::setCommNodeIndex()にてセット
    // --
    vector<CNode*>    mvSendNode;//送信Node
    vector<CNode*>    mvRecvNode;//受信Node

    // CommNodeIDを(Send,Recv別に)一列にならべた配列 -> CCommElement::setCommNodeIndex()にてセット
    // --
    vuint mvSendCommNodeID;// CommNodeIDの一覧(vector配列):SendNode
    vuint mvRecvCommNodeID;// CommNodeIDの一覧(vector配列):RecvNode


    // NodeIDからCommNodeIDを取得するHash
    map<uint, uint, less<uint> > mmCommNodeIX;// NodeID => CommNodeID(Index番号)
    map<uint, uint, less<uint> > mmCommElementIX;// ElementID => CommElementID(Index番号)
    
    
    // ◎ 計算に不要なNode,Element
    //  => CMesh内で計算用に配列を入れ替えるときに利用
    //  => this->prolongation()で収集.
    // --
    vector<CNode*>    mvDNode;   //Dead Node (useless Node)
    vector<CElement*> mvDElement;//Dead Element (useless Element) <= mvDCommElementから取得


    Utility::CLogger* mpLogger;
    
public:
    // 自身の通信番号(Communication_ID)
    void setCommID(const uint& comID){ mCommID= comID;}
    uint& getCommID(){ return mCommID;}
    
    // 所属する計算領域の番号(mpiのrank)
    void setRankID(const uint& rank){ mRankID= rank;}
    uint& getRankID(){ return mRankID;}
    
    // 通信相手の計算領域番号
    void setTransmitRankID(const uint& rank){ mTransmitRankID = rank;}
    uint& getTransmitRankID(){ return mTransmitRankID;}
    
    
    //    // Node rank 一覧 
    //    //  => '09.10.02 ファイル読み込み時(CommElementのNodeRank設定)にmvNodeRankを使用
    //    // --
    //    void pushNodeRank(const uint& rank){  mvNodeRank.push_back(rank);}


    
    // 全体のCCommElement
    // --
    void reserveCommElementAll(const uint& res_size){ mvCommElementAll.reserve(res_size);}
    void setCommElementAll(CCommElement* pCommElem){ mvCommElementAll.push_back(pCommElem);}
    uint getNumOfCommElementAll(){ return mvCommElementAll.size();}
    CCommElement* getCommElementAll(const uint& all_com_index){ return mvCommElementAll[all_com_index];}
    CElement* getElementAll(const uint& all_com_index){ return mvCommElementAll[all_com_index]->getElement();}
    
    // CommElement
    // --
    // 通信に使用される,オーバーラップ要素(共有要素):CCommElement
    //
    uint getNumOfCommElement(){ return mvCommElement.size();}
    CCommElement* getCommElement(const uint& com_index){ return mvCommElement[com_index];}//通信領域(Comm)のCommIndex(i)番目の"Comm要素"を取得
    CElement* getElement(const uint& com_index){ return mvCommElement[com_index]->getElement();}//通信領域(Comm)のCommIndex(i)番目の"要素"を取得


    // ・DNode:計算に使用されないNode
    // ・DCommElement:計算にしようされないCCommElement(prolongation後に発生)
    //   { 通信にも使用されない. }
    // --
    uint getNumOfDNode(){ return mvDNode.size();}
    //vector<CNode*> getDNode(){ return mvDNode;}
    CNode* getDNode(const uint& dcom_index){ return mvDNode[dcom_index];}

    uint getNumOfDCommElement(){ return mvDCommElement.size();}
    CCommElement* getDCommElement(const uint& dcom_index){ return mvDCommElement[dcom_index];}
    //CElement* getDElement(const uint& dcom_index){ return mvDCommElement[dcom_index]->getElement();}
    CElement* getDElement(const uint& dcom_index){ return mvDElement[dcom_index];}//ソートされたDElementを使用


    // 送受信ノード全体
    // --
    uint getNumOfNode(){ return mvNode.size();}
    void reserveNode(const uint& res_size){ mvNode.reserve(res_size);}//Node配列確保(ファイル入力時,prolongation時に使用)
    void setNode(CNode* pNode){ mvNode.push_back(pNode);}             //Nodeのセット(ファイル入力時,prolongation時に使用)
    CNode* getNode(const uint& com_index){ return mvNode[com_index];} //CommMesh内のノードIndex番号(CommNodeID)のノードを提供
    
    
    // 送信ノード
    // --
    //void reserveSendNode(const uint& res_size){ mvSendNode.reserve(res_size);}//SendNode配列確保(ファイル入力時,prolongation時に使用)
    void setSendNode(CNode* pNode, const uint& comNodeID){ mvSendCommNodeID.push_back(comNodeID); mvSendNode.push_back(pNode);}    //SendNodeのセット(ファイル入力時,prolongation時に使用)
    CNode* getSendNodeIX(const uint& index){ return mvSendNode[index];}//通信領域(Comm)のSendNodeのindex番目の送信ノードを提供
    uint getNumOfSendNode(){ return mvSendNode.size();}
    uint& getSendCommNodeID(const uint& index){ return mvSendCommNodeID[index];}

    // 受信ノード
    // --
    //void reserveRecvNode(const uint& res_size){ mvRecvNode.reserve(res_size);}//RecvNode配列確保(ファイル入力時,prolongation時に使用)
    void setRecvNode(CNode* pNode, const uint& comNodeID){ mvRecvCommNodeID.push_back(comNodeID); mvRecvNode.push_back(pNode);}    //RecvNodeのセット(ファイル入力時,prolongation時に使用)
    CNode* getRecvNodeIX(const uint& index){ return mvRecvNode[index];}//通信領域(Comm)のRecvNodeのindex番目の受信ノードを提供
    uint getNumOfRecvNode(){ return mvRecvNode.size();}
    uint& getRecvCommNodeID(const uint& index){ return mvRecvCommNodeID[index];}

    // rank
    // --
    void resizeNodeRank(const uint& res_size){ mvNodeRank.resize(res_size);}
    void setNodeRank(const uint& commNodeID, const uint& rank){ mvNodeRank[commNodeID]= rank;}//最初のファイル入力向け
    uint& getNodeRank(const uint& commNodeID){ return  mvNodeRank[commNodeID];}


    // prolongation
    // --
    // 通信に使用するCommElementと,使用しないCommElementの振り分け:通信に用いないCommElementはカップラー向け
    // CommElementIDはIndex番号(通し番号)<= 全体に割り振る
    // --
    void AllocateCommElement();
    
    // Communication(CommMesh)内でのIndex生成(Global番号相当),mvNode取得,Send-Recvの取得
    // --
    void setupAggCommElement(vector<CElement*> &vElement);//Vertexの要素Index集合を利用して,CommElementAllの隣接情報をセットアップ
    void sortCommNodeIndex();// mvNodeの取得,CommNodeIDの生成: CommNodeIDはIndex番号(通し番号),Send,Recvの収集
    
    // mapデータのセットアップ
    // --
    // mmCommNodeIX : NodeID => CommNodeID(Index番号)のHash
    // mmCommElementIX : ElementID => CommElementID(Index番号)のHash
    void setupMapID2CommID();

};
}
#endif	/* _COMMMesh_H */

















