/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/CommNode.h
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
#include "TypeDef.h"
#include "Vertex.h"
#include "Node.h"
namespace pmw
{
#ifndef _COMMNODE_H
#define	_COMMNODE_H
class CCommNode:public CVertex
{
public:
    CCommNode();
    virtual ~CCommNode();
protected:
    CNode* mpNode;
    uiint mLevel;
public:
    void   setNode(CNode* pNode) {
        mpNode= pNode;
    }
    CNode* getNode() {
        return mpNode;
    }
    uiint&  getNodeID() {
        return mpNode->getID();
    }
    void setLevel(const uiint& level) {
        mLevel= level;
    }
    uiint& getLevel() {
        return mLevel;
    }
};
#endif	/* _COMMNODE_H */
}
