
#include <vector>

//
//  CommunicationTable.cpp
//
//
//
//                      2009.06.19
//                      2009.06.19
//                      k.Takeda
#include "CommunicationTable.h"
using namespace pmw;

#include <iostream>
// construct & destruct
//
CCommTable::CCommTable()
{
    ;
}
CCommTable::~CCommTable()
{
    //debug
    std::cout << "~CCommTable" << std::endl;
}

// 通信する対象領域 => 通信対象領域数に応じた,通信テーブル領域確保
//
void CCommTable::reserveCommDomain(const uint& res_size)
{
    mvCommDomain.reserve(res_size);

    // Mesh内部のローカル要素-節点を"通信領域数"分確保(resize)
    mvvElement.resize(res_size);
    mvvSendNode.resize(res_size);
    mvvRecvNode.resize(res_size);

    //    // 外部通信のためのグローバルID管理 => Mesh内部とのペアを"通信領域数"分確保(resize)
    //    mvvGlobalElementID.resize(res_size);
    //    mvvGlobalSendNodeID.resize(res_size);
    //    mvvGlobalRecvNodeID.resize(res_size);
}

void CCommTable::setCommDomain(const uint& domID)
{
    mvCommDomain.push_back(domID);// 配列

    mmCommDomain[domID]= mvCommDomain.size()-1;// 通信領域ID => 通信領域Index を取得する,Hash
}

uint CCommTable::getNumOfCommDomain()
{
    return mvCommDomain.size();
}


// 通信要素(共有要素)
//
void CCommTable::reserveElement(const uint& domID, const uint& res_size)
{
    uint domIndex = mmCommDomain[domID];
    mvvElement[domIndex].reserve(res_size);
}

void CCommTable::setElement(const uint& domID, CElement* pElem)
{
    uint domIndex = mmCommDomain[domID];
    mvvElement[domIndex].push_back(pElem);
}

CElement* CCommTable::getElement(const uint& domID, const uint& i)
{
    uint domIndex = mmCommDomain[domID];

    return mvvElement[domIndex][i];
}

// 送信ノード
//
void CCommTable::reserveSendNode(const uint& domID, const uint& res_size)
{
    uint domIndex = mmCommDomain[domID];
    mvvSendNode[domIndex].reserve(res_size);
}

void CCommTable::setSendNode(const uint& domID, CNode* pNode)
{
    uint domIndex(mmCommDomain[domID]);
    mvvSendNode[domIndex].push_back(pNode);
}

CNode* CCommTable::getSendNode(const uint& domID, const uint& i)
{
    uint domIndex(mmCommDomain[domID]);

    return mvvSendNode[domIndex][i];
}

// 受信ノード
//
void CCommTable::reserveRecvNode(const uint& domID, const uint& res_size)
{
    uint domIndex(mmCommDomain[domID]);
    mvvRecvNode[domIndex].reserve(res_size);
}

void CCommTable::setRecvNode(const uint& domID, CNode* pNode)
{
    uint domIndex(mmCommDomain[domID]);
    mvvRecvNode[domIndex].push_back(pNode);
}

CNode* CCommTable::getRecvNode(const uint& domID, const uint& i)
{
    uint domIndex(mmCommDomain[domID]);
    return mvvRecvNode[domIndex][i];
}



















