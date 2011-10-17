/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/Node.h
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "TypeDef.h"
#include "NodeType.h"
#include "Logger.h"
#include "Vertex.h"
namespace pmw{
#ifndef NODE_HH_
#define NODE_HH_
class CNode:public CVertex{
public:
    CNode();
    virtual ~CNode();
protected:
    vector<CNode*> mvParentNode;
    bool mbSComm;
public:
    virtual uiint& getType()=0;
    void markingSCommNode();
    bool isSCommNode();
    virtual void setScalarDOF(const uiint& nNDOF)=0;
    virtual void setVectorDOF(const uiint& nNDOF)=0;
    virtual uiint& getScalarDOF()=0;
    virtual uiint& getVectorDOF()=0;
    virtual uiint getTotalDOF()=0;
    vector<CNode*>& getParentNode(){ return mvParentNode;}
    CNode* getParentNode(const uiint& index){ return mvParentNode[index];}
    void reserveParentNode(const uiint& res_size){ mvParentNode.reserve(res_size);}
    void addParentNode(CNode* pNode){ mvParentNode.push_back(pNode);}
    uiint getNumOfParentNode(){ return mvParentNode.size();}
};
#endif
}
