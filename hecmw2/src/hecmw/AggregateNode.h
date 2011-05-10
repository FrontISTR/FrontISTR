/* 
 * File:   AggregateNode.h
 * Author: ktakeda
 *
 * Created on 2009/08/11, 14:57
 */
#ifndef _AGGREGATENODE_H_3dc2add2_f427_4db1_b0e3_6437869e9236
#define	_AGGREGATENODE_H_3dc2add2_f427_4db1_b0e3_6437869e9236

#include "CommonStd.h"
#include "TypeDef.h"

#include "Node.h"

namespace pmw{
class CAggregateNode{
public:
    CAggregateNode();
    virtual ~CAggregateNode();

protected:
    uint mID;//コア・ノード インデックス番号
    vector<CNode*> mvNode;

public:
    // コアになるノードのID (インデックス番号)
    void setID(const uint& index){ mID= index;}
    uint& getID(){ return mID;}

    // コア・ノード周囲のノード => Aggregate Node
    void reserveNode(const uint& res_size){ mvNode.reserve(res_size);}
    void setNode(CNode* pNode);//既に取得済みかどうかの判定の上でpush_back.

    uint getNumOfNode(){ return mvNode.size();}
    CNode* getNode(const uint& i){ return mvNode[i];}


};
}
#endif	/* _AGGREGATENODE_H */
