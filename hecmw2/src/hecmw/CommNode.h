/* 
 * File:   CommNode.h
 * Author: ktakeda
 *
 * Created on 2010/03/02, 17:25
 */
#include "TypeDef.h"
#include "Vertex.h"

#include "Node.h"

namespace pmw{
#ifndef _COMMNODE_H
#define	_COMMNODE_H
class CCommNode:public CVertex{
public:
    CCommNode();
    virtual ~CCommNode();

protected:
    CNode* mpNode;

    uiint mLevel;

public:
    // root class Vertex
    // --
    // 参照しているNode
    void   setNode(CNode* pNode){ mpNode= pNode;}
    CNode* getNode(){ return mpNode;}
    uiint&  getNodeID(){ return mpNode->getID();}

    void setLevel(const uiint& level){ mLevel= level;}
    uiint& getLevel(){ return mLevel;}
};
#endif	/* _COMMNODE_H */
}


