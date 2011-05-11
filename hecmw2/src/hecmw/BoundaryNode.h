/* 
 * File:   BoundaryNode.h
 * Author: ktakeda
 *
 * BoundaryVolumeMesh, BoundaryFaceMesh, BoundaryEdgeMesh などに利用するBNode
 * BoundaryNodeMeshは、別クラスを利用(BoundarySBNode)
 *
 * CBndVertexを親とする派生クラス
 *
 * --------------------
 * 荷重
 * 変位
 * 速度
 * 加速度
 * 固定温度(熱応力に利用)
 * 集中熱流束
 * --------------------
 *
 * Created on 2009/05/13, 15:04
 */
#include "TypeDef.h"
#include "BndVertex.h"//親
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
    // mID => BndVertex
    // mnBndType =>  BoundaryMesh, Dirichlet | Neumann は,BoundaryMeshが知っていれば良い.

    uiint mMGLevel; // 自身の所属Level

    //    //----
    //    // DOF番号管理
    //    //----
    //    vuint mvDOF;//境界のDOF番号
    
    vector<map<uiint, double, less<uiint> > > mvValue;//[階層Level][DOF]別の境界値
    
    CNode *mpNode;//境界条件が付与されるNode(Mesh-Node)

public:
    //階層Level
    void setMGLevel(const uiint& mgLevel){ mMGLevel= mgLevel;}
    uiint& getMGLevel(){ return mMGLevel;}

    //境界値
    void resizeValue(const uiint& numOfDiffLevel);//Level数分の配列確保(計算Levelと節点MGLevelの差)

    void initValue(const uiint& dof, const uiint& mgLevel);//初期化(ゼロ・クリア)
    void setValue(const uiint& dof, const uiint& mgLevel, const double& val);//代入
    void addValue(const uiint& dof, const uiint& mgLevel, const double& val);//加算
    double& getValue(const uiint& dof, const uiint& mgLevel);//境界値の提供


    //メッシュ-節点
    void setNode(CNode *pNode){ mpNode= pNode;}
    CNode* getNode(){ return mpNode;}
};
#endif	/* _BOUNDARYNODE_H */
}



