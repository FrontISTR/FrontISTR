/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   AggregateNode.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
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
    uint mID;
    vector<CNode*> mvNode;
public:
    void setID(const uint& index){ mID= index;}
    uint& getID(){ return mID;}
    void reserveNode(const uint& res_size){ mvNode.reserve(res_size);}
    void setNode(CNode* pNode);
    uint getNumOfNode(){ return mvNode.size();}
    CNode* getNode(const uint& i){ return mvNode[i];}
};
}
#endif	/* _AGGREGATENODE_H */
