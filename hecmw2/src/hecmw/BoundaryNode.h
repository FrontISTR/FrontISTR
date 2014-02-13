/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryNode.h
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
#include "BndVertex.h"
#include "BoundaryType.h"
#include "Node.h"
#include "Logger.h"
#include <map>
namespace pmw
{
#ifndef _BOUNDARYNODE_H_
#define	_BOUNDARYNODE_H_
class CBoundaryNode:public CBndVertex
{
public:
    CBoundaryNode();
    virtual ~CBoundaryNode();
protected:
    uiint mMGLevel;
    vector<map<uiint, double, less<uiint> > > mvValue;//----- 境界値(ノイマン,ディレクレ)，
    vector<map<uiint, double, less<uiint> > > mvEntValue;//--- ディレクレ数式計算 基礎データ

    vector<map<uiint, bool> > mvMarkingValue;//---Rcap_Refineのディレクレ分配時に利用

    CNode *mpNode;
public:
    void setMGLevel(const uiint& mgLevel) {
        mMGLevel= mgLevel;
    }
    uiint& getMGLevel() {
        return mMGLevel;
    }
    //境界値 & 基礎データ共通
    ////void resizeValue(const uiint& nNumOfDiffLevel, vuint& vDOF);//初期化(0.0)も含める
    void resizeValue(const uiint& nNumOfDiffLevel);
    void initValue(const uiint& dof, const uiint& mgLevel);//------ Neumann値のaddのための初期化(0.0)
    void initRcapBool(const uiint& dof, const uiint& mgMaxLevel);//--- Rcapのディレクレ処理のため

    //境界値
    void setValue(const uiint& dof, const uiint& mgLevel, const double& val);
    void addValue(const uiint& dof, const uiint& mgLevel, const double& val);
    double& getValue(const uiint& dof, const uiint& mgLevel);

    bool isSetupValue(const uiint& dof, const uiint& mgLevel);//---Rcap_Refineのディレクレ分配時に利用

    //基礎データ(ディレクレ数式演算の基礎データ)
    void setEntValue(const uiint& dof, const uiint& mgLevel, const double& val);
    double& getEntValue(const uiint& dof, const uiint& mgLevel);

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
#endif	/* _BOUNDARYNODE_H */
}
