/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   CommNode.h
|
|                     Written by T.Takeda,    2010/06/01
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
namespace pmw{
#ifndef _COMMNODE_H
#define	_COMMNODE_H
class CCommNode:public CVertex{
public:
    CCommNode();
    virtual ~CCommNode();
protected:
    CNode* mpNode;
    uint mLevel;
public:
    void   setNode(CNode* pNode){ mpNode= pNode;}
    CNode* getNode(){ return mpNode;}
    uint&  getNodeID(){ return mpNode->getID();}
    void setLevel(const uint& level){ mLevel= level;}
    uint& getLevel(){ return mLevel;}
};
#endif	/* _COMMNODE_H */
}
