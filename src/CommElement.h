///*
// * File:   CommElement.h
// * Author: ktakeda
// *
// * Created on 2009/08/28, 13:19
// */
#include "TypeDef.h"

#include "Node.h"
#include "Element.h"

#include "NodeConnectNodeTree.h"// rank判定にTreeを使用予定.##
#include "EdgeTree.h"// rank判定にTreeを使用.
#include "FaceTree.h"// rank判定にTreeを使用.
#include "EdgeFaceTree.h"// rank判定にTreeを使用.
#include "ProgElementTree.h"//prolongationするときの親の辺・面のノードと子要素の頂点位置関係ツリー

#include "Logger.h"
#include "ElementType.h"

namespace pmw{
#ifndef _COMMELEMENT_H
#define	_COMMELEMENT_H
class CCommElement{
public:
    CCommElement();
    virtual ~CCommElement();

protected:
    uint mID;
    CElement* mpElement;//内包するMeshの要素

    //各頂点に共有CommElement集合を作る
    vector<vector<CCommElement*> > mvvAggCommElem;//<= 全てのCommElement(DCommも含む)
    vvuint mvvNeibCommElemVert;                   //共有CommElemの接続先の頂点番号
    

    // Send-Recv 局所Node番号
    // --
    //
    vbool mvbSend;//頂点番号のノードがSend？ => CommMeshのSendノード収集に使用
    vbool mvbRecv;//頂点番号のノードがRecv？ => CommMeshのRecvノード収集に使用
    vbool mvbOther;//頂点番号のノードは無関係なRank？
    vector<CNode*> mvSendNode;// myRank所属ノード
    vector<CNode*> mvRecvNode;// transmitRank所属ノード
    vector<CNode*> mvOtherNode;//myRankでもtransmitRankでもない,別のRankに所属

    // Node Index (CommMesh内でのNodeのIndex)<= Global Indexの代替
    // --
    vuint mvCommNodeIndex;//CommMesh内のグローバルIndex ,局所番号順にIndexを並べてある.

    // Communicationに用いるか否か(CommElementか,DCommElementか)
    bool mbCommunication;
    // CommElementから通常のElementに戻った(全ての頂点がmyRank)
    bool mbRCommunication;

    // Node選択されたかマーキングするためのboolean
    // --
    vector<bool> mvbNodeIXCheck;//グローバルIndex生成時に,Nodeが選択されたかマーキング
    vector<bool> mvbDNodeMarking;//DNodeとして選ばれたかマーキング

    Utility::CLogger *mpLogger;


    // --
    // Node:rank番号(領域番号)
    // --
    vuint mvNodeRank;//頂点ごとの所属計算領域(DomID <=rank)
    vuint mvEdgeRank;//辺 ごとの所属計算領域(DomID <=rank)
    vuint mvFaceRank;//面 ごとの所属計算領域(DomID <=rank)
     uint  mVolRank; //体積中心の所属計算領域(DomID <=rank)
public:
    // prolongation :Factoryで利用される
    void setupProgNodeRank(const uint& mgLevel);//ProgElemのNodeRankの決定.<= Edge,Face,Volumeのランクを決定と同義

    // ID
    void setID(const uint& id){ mID= id;}
    uint& getID(){ return mID;}


    // 要素形状
    virtual uint getShapeType()=0;
    virtual uint getBaseShapeType()=0;
    const uint& getNumOfVert(){ return mpElement->getNumOfNode();}
    const uint& getNumOfEdge(){ return mpElement->getNumOfEdge();}
    const uint& getNumOfFace(){ return mpElement->getNumOfFace();}


    // 頂点のRank設定
    void setNodeRank(const uint& ivert, const uint& rank){ mvNodeRank[ivert]=rank;}//ファイル読み込み時(Level=0時の各Nodeのランク)
    uint& getNodeRank(const uint& ivert){ return mvNodeRank[ivert];}


    // 所有Element
    // --
    void setElement(CElement* pElem){ mpElement= pElem;}
    CElement* getElement(){ return mpElement;}


    // Node Index
    // --
    // comID:CommMeshグローバルNodeIndex(comID)のカウントアップ,Index値のセット,Marking
    // vNode:CommMeshのmvNodeを引数で受ける.
    // --
    void setCommNodeIndex(const uint& ivert, uint& comID, vector<CNode*>& vCommMeshNode, vuint& vCommMeshNodeRank);
protected:
    void setNeibCommNodeIndex(const uint& ivert, const uint& comID);//隣接CommElemへのCommNodeIndex割り振り
public:
    uint& getCommNodeIndex(const uint& ivert){ return mvCommNodeIndex[ivert];}

    // DNode判定後に提供
    // => mbCommunication==false  但し-> Nodeを共有しているCommElemが,mbCommunication==trueの場合は除外
    // --
    void getDNode(const uint& ivert, vector<CNode*>& vDNode);
protected:
    void markingDNode(const uint& ivert);

public:
    bool isMarkingDNode(const uint& ivert){ return mvbDNodeMarking[ivert];}

    // Nodeを提供
    // --
    CNode* getNode(const uint& ivert){ return mpElement->getNode(ivert);}
    // Send,Recv Nodeを提供
    // --
    vector<CNode*>& getSendNode(){ return mvSendNode;}
    vector<CNode*>& getRecvNode(){ return mvRecvNode;}
    uint getNumOfSendNode(){ return mvSendNode.size();}
    uint getNumOfRecvNode(){ return mvRecvNode.size();}
    CNode* getSendNode(const uint& i){ return mvSendNode[i];}
    CNode* getRecvNode(const uint& i){ return mvRecvNode[i];}

    
    //頂点共有Element
    void setAggCommElement(const uint& ivert, CCommElement* pComElem){ mvvAggCommElem[ivert].push_back(pComElem);}
    void setNeibCommElemVert(const uint& ivert, const uint& neibVert){ mvvNeibCommElemVert[ivert].push_back(neibVert);}
    
    
    // rank選別,CommElement選択(Tetra,Triangleは頂点間の関係だけで記述可能)
    // --
    bool isCommElement(){ return mbCommunication;}//通信に用いる要素か否か
    bool isRCommElement(){ return mbRCommunication;}//計算のみの要素に復帰した要素か否か
    void sortNodeRank(const uint& myRank, const uint& transRank);//Send,Recv,Otherに選別,isCommElementのbool値設定
    
    
    //debug method 
    virtual bool isTypeCoincidence()=0;//CommElementと所有Elemの型が一致しているか？

    // prolongation CommElement のRank提供
    // --
    vuint& getEdgeRank(){ return mvEdgeRank;}
     uint& getEdgeRank(const uint& iedge){ return mvEdgeRank[iedge];}
    vuint& getFaceRank(){ return mvFaceRank;}
     uint& getFaceRank(const uint& iface){ return mvFaceRank[iface];}
     uint& getVolRank(){ return mVolRank;}
};
#endif	/* _COMMELEMENT_H */
}
