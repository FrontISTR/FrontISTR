/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/AggregateNode.h
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
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
namespace pmw
{
class CAggregateNode
{
public:
    CAggregateNode();
    virtual ~CAggregateNode();
protected:
    uiint mID;
    vector<CNode*> mvNode;
public:
    void setID(const uiint& index) {
        mID= index;
    }
    uiint& getID() {
        return mID;
    }
    void reserveNode(const uiint& res_size) {
        mvNode.reserve(res_size);
    }
    void setNode(CNode* pNode);
    uiint getNumOfNode() {
        return mvNode.size();
    }
    CNode* getNode(const uiint& i) {
        return mvNode[i];
    }
};
}
#endif	/* _AGGREGATENODE_H */
