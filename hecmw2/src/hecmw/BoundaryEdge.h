/* 
 * File:   BoundaryEdge.h
 * Author: ktakeda
 *
 * Created on 2010/04/13, 19:39
 */
#include "BoundaryParts.h"
#include <utility>//pair

typedef std::pair<pmw::CBoundaryNode*, pmw::CBoundaryNode*> PairBNode;

namespace pmw{
#ifndef _BOUNDARYEDGE_H
#define	_BOUNDARYEDGE_H
class CBoundaryEdge:public CBoundaryParts{
public:
    CBoundaryEdge();
    virtual ~CBoundaryEdge();

protected:
    uint mnElemEdgeID;//エンティティ:Mesh-Element-Edge-ID

    CBoundaryNode *mpEdgeBNode;// 辺中心のBNode

    //PairBNode mPairBNode;

    vector<CBoundaryEdge*> mvProgBEdge;

    double mLength;
public:
    // 要素辺ID(Edge ID)
    void setElementEdgeID(const uint& id){ mnElemEdgeID= id;}
    uint& getElementEdgeID(){ return mnElemEdgeID;}

    // 辺形状( Beam )
    void setBEdgeShape(const uint& elemType);
    uint& getBEdgeShape(){ return mnShapeType;}

    virtual uint getNumOfVert();

    // 辺中心のBNode
    void setEdgeBNode(CBoundaryNode *pBNode);
    CBoundaryNode* getEdgeBNode(){ return mpEdgeBNode;}

    // 辺の両端のBNode
    PairBNode getPairBNode();

    // 辺Nodeを,辺BoundaryNodeにセット
    void setupNode();

    // 線分の距離
    double& calcLength();
    double& getLength(){ return mLength;}


    // 辺の再分割
    // ----
    void refine(uint& countID, const vuint& vDOF);

    vector<CBoundaryEdge*>& getProgParts(){ return mvProgBEdge;}

    void replaceEdgeBNode();// 2次要素の場合、辺BNodeをmvBNodeに移設.

    void distDirichletVal(const uint& dof, const uint& mgLevel, const uint& nMaxMGLevel);//上位グリッドBNodeへのディレクレ値の分配
};
#endif	/* _BOUNDARYEDGE_H */
}



