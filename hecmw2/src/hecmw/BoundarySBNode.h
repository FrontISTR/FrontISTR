/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundarySBNode.h
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
#include "BndVertex.h"
#include "Node.h"
namespace pmw
{
#ifndef _BOUNDARYSBNODE_H
#define	_BOUNDARYSBNODE_H
class CBoundarySBNode:public CBndVertex
{
public:
    CBoundarySBNode();
    virtual ~CBoundarySBNode();
protected:
    vuint mvDOF;
    map<uiint, double, less<uiint> > mmValue;
    map<uiint, double, less<uiint> > mmEntValue;//ディレクレ数式処理 基礎データ
    CNode *mpNode;
public:
    void addDOF(const uiint& dof);
    uiint& getDOF(const uiint& index);
    uiint getNumOfDOF();

    //境界値
    void setValue(const uiint& dof, const double& val);
    double& getValue(const uiint& dof);

    //基礎データ
    void setEntValue(const uiint& dof, const double& val);
    double& getEntValue(const uiint& dof);

    void setNode(CNode *pNode) {
        mpNode= pNode;
    }
    CNode* getNode() {
        return mpNode;
    }

    double& getX();
    double& getY();
    double& getZ();
};
#endif	/* _BOUNDARYSBNODE_H */
}
