/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/BoundarySBNode.h
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
#include "BndVertex.h"
#include "Node.h"
namespace pmw{
#ifndef _BOUNDARYSBNODE_H
#define	_BOUNDARYSBNODE_H
class CBoundarySBNode:public CBndVertex{
public:
    CBoundarySBNode();
    virtual ~CBoundarySBNode();
protected:
    vuint mvDOF;
    map<uiint, double, less<uiint> > mmValue;
    CNode *mpNode;
public:
    void addDOF(const uiint& dof);
    uiint& getDOF(const uiint& index);
    uiint getNumOfDOF();
    void setValue(const uiint& dof, const double& val);
    double& getValue(const uiint& dof);
    void setNode(CNode *pNode){ mpNode= pNode;}
    CNode* getNode(){ return mpNode;}
};
#endif	/* _BOUNDARYSBNODE_H */
}
