/* 
 * File:   BoundarySBNode.h
 * Author: ktakeda
 *
 * Levelを持たないBoundaryNode
 * # DOF別の値を管理するだけ
 *
 * Created on 2010/06/25, 15:13
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
    //----
    // DOF番号管理
    //----
    vuint mvDOF;//境界のDOF番号
    map<uiint, double, less<uiint> > mmValue;//[DOF]別の境界値 :Level無し

    CNode *mpNode;//境界条件が付与されるNode(Mesh-Node)

public:
    //DOF
    void addDOF(const uiint& dof);
    uiint& getDOF(const uiint& index);
    uiint getNumOfDOF();
    
    
    //境界値
    void setValue(const uiint& dof, const double& val);//代入
    double& getValue(const uiint& dof);//提供

    //メッシュ-節点
    void setNode(CNode *pNode){ mpNode= pNode;}
    CNode* getNode(){ return mpNode;}
};
#endif	/* _BOUNDARYSBNODE_H */
}

