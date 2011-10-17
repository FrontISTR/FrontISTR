/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/BoundaryNode.h
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
#include "BndVertex.h"
#include "BoundaryType.h"
#include "Node.h"
#include "Logger.h"
#include <map>
namespace pmw{
#ifndef _BOUNDARYNODE_H_
#define	_BOUNDARYNODE_H_
class CBoundaryNode:public CBndVertex{
public:
    CBoundaryNode();
    virtual ~CBoundaryNode();
protected:
    uiint mMGLevel; 
    vector<map<uiint, double, less<uiint> > > mvValue;
    CNode *mpNode;
public:
    void setMGLevel(const uiint& mgLevel){ mMGLevel= mgLevel;}
    uiint& getMGLevel(){ return mMGLevel;}
    void resizeValue(const uiint& numOfDiffLevel);
    void initValue(const uiint& dof, const uiint& mgLevel);
    void setValue(const uiint& dof, const uiint& mgLevel, const double& val);
    void addValue(const uiint& dof, const uiint& mgLevel, const double& val);
    double& getValue(const uiint& dof, const uiint& mgLevel);
    void setNode(CNode *pNode){ mpNode= pNode;}
    CNode* getNode(){ return mpNode;}
};
#endif	/* _BOUNDARYNODE_H */
}
