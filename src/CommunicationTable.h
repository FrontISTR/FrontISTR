///*
// * File:   CommunicationTable.h
// * Author: ktakeda
// *
// *  通信インデックス順位： 要素順 -> 要素内の局所番号順
// *
// * Modify     2009/08/21
// * Created on 2009/06/19, 19:42
// */

#ifndef _COMMTABLE_H_39e5ebb0_da94_4efc_8616_8c2c5918cc4f
#define _COMMTABLE_H_39e5ebb0_da94_4efc_8616_8c2c5918cc4f


#include "Element.h"
#include "Node.h"

#include <map>

namespace pmw{
class CCommTable{
private:
    CCommTable();
public:
    static CCommTable* Instance(){
        static CCommTable commTable;
        return &commTable;
    }
    virtual ~CCommTable();
    
protected:
    // 自身の領域番号
    uint mDomainID;
    // 通信対象の領域(複数存在する)
    vuint mvCommDomain;

    // ドメインID .-> ドメインIndex
    map<uint, uint, less<uint> > mmCommDomain;
    
    
    // -- Memo --
    // 通信インデックス：要素順 -> 要素内の局所番号順
    // --
    
    // Mesh内部のローカル要素-節点
    vector<vector<CElement*> > mvvElement;//オーバーラップ要素 => 要素を共有する場合は共有要素
    vector<vector<CNode*> >  mvvSendNode;//送信Node
    vector<vector<CNode*> >  mvvRecvNode;//受信Node

    //    // 外部通信のためのグローバルID管理 => Mesh内部とのペア
    //    vvuint mvvGlobalElementID;
    //    vvuint mvvGlobalSendNodeID;
    //    vvuint mvvGlobalRecvNodeID;

public:
    // 自身の領域番号
    void setDomainID(const uint& domID){ mDomainID= domID;}
    uint& getDomainID(){ return mDomainID;}


    // 通信する対象領域 => 通信対象領域数に応じた,通信テーブル領域確保
    // --
    void reserveCommDomain(const uint& res_size);//{要素,送受信ノード}の通信領域数に対応する配列確保も行う.
    void setCommDomain(const uint& domID);
    uint getNumOfCommDomain();
    uint& getCommDomainID(const uint& index){ return mvCommDomain[index];}//配列から通信領域IDを取得
    uint& getCommDomainIndex(const uint& domID){ return mmCommDomain[domID];}//Hashから通信領域IDを取得
    
    
    
    //    // グローバルID
    //    // --
    //    // オーバーラップ要素
    //    void reserveGlobalElementID(const uint& domID, const uint& res_size){ uint index(mmCommDomain[domID]);; mvvGlobalElementID[index].reserve(res_size);}
    //    void setGlobalElementID(const uint& domID, const uint& id){ uint index(mmCommDomain[domID]);; mvvGlobalElementID[index].push_back(id);}
    //    // 送信ノード
    //    void reserveGlobalSendNodeID(const uint& domID, const uint& res_size){ uint index(mmCommDomain[domID]);; mvvGlobalSendNodeID[index].reserve(res_size);}
    //    void setGlobalSendNodeID(const uint& domID, const uint& id){ uint index(mmCommDomain[domID]);; mvvGlobalSendNodeID[index].push_back(id);}
    //    // 受信ノード
    //    void reserveGlobalRecvNodeID(const uint& domID, const uint& res_size){ uint index(mmCommDomain[domID]);; mvvGlobalRecvNodeID[index].reserve(res_size);}
    //    void setGlobalRecvNodeID(const uint& domID, const uint& id){ uint index(mmCommDomain[domID]);; mvvGlobalRecvNodeID[index].push_back(id);}


    // -- Memo --
    // 通信インデックス：要素順 -> 要素内の局所番号順
    // --
    
    // Mesh
    // --
    // オーバーラップ要素(共有要素), { 引数:domID => 通信対象の領域ID }
    //
    void reserveElement(const uint& domID, const uint& res_size);
    void setElement(const uint& domID, CElement* pElem);
    CElement* getElement(const uint& domID, const uint& i);//指定通信領域のi番目の要素を取得
    // 送信ノード
    void reserveSendNode(const uint& domID, const uint& res_size);
    void setSendNode(const uint& domID, CNode* pNode);
    CNode* getSendNode(const uint& domID, const uint& i);//指定通信領域のi番目の送信ノードを取得
    // 受信ノード
    void reserveRecvNode(const uint& domID, const uint& res_size);
    void setRecvNode(const uint& domID, CNode* pNode);
    CNode* getRecvNode(const uint& domID, const uint& i);//指定通信領域のi番目の受信ノードを取得
};
}
#endif	/* _COMMTABLE_H */

















