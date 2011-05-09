/* 
 * File:   BoundaryFace.h
 * Author: ktakeda
 *
 * Created on 2009/05/18, 15:11
 */
#include "BoundaryParts.h"
#include <utility>//pair

#include "EdgeTree.h"

#include "ElementProperty.h"

typedef std::pair<pmw::CBoundaryNode*, pmw::CBoundaryNode*> PairBNode;

namespace pmw{
#ifndef _BOUNDARYFACE_H_c5c754e6_
#define	_BOUNDARYFACE_H_c5c754e6_   
class CBoundaryFace:public CBoundaryParts{
public:
    CBoundaryFace();
    virtual ~CBoundaryFace();

protected:
    uint mnElemFaceID;//エンティティID:Mesh-Element-Face-ID

    //圧力(垂直荷重)
    //面方向の面力(粘性:Viscous Force)
    //熱流束(面)

    vuint mvEdgeNeibFace;//辺に隣接する面(自身は含まない)
    
    //PairBNode mPairBNode;//辺の両端のBNode, テンポラリー変数

    vector<CBoundaryNode*> mvEdgeBNode;
    CBoundaryNode *mpFaceBNode;

public:
    // 要素面ID(Face ID)
    void setElementFaceID(const uint& id){ mnElemFaceID= id;}
    uint& getElementFaceID(){ return mnElemFaceID;}

    // 面形状( Quad | Triangle )
    void setBFaceShape(const uint& elemType);
    uint& getBFaceShape(){ return mnShapeType;}

    virtual uint getNumOfVert();
    
protected:
    bool* mvbMarkingEdge;// 辺が使用済みか否か.

public:
    void markingEdge(const uint& iedge);
    bool isMarkingEdge(const uint& iedge){ return mvbMarkingEdge[iedge];}
    
    void setEdgeNeibFace(const uint& iedge, const uint& neibFaceID);
    void setEdgeBNode(const uint& iedge, CBoundaryNode *pEdgeBNode);

    uint getNumOfEdge();
    PairBNode getPairBNode(const uint& iedge);
    uint& getEdgeID(PairBNode& pairBNode);

    void setFaceBNode(CBoundaryNode *pFaceBNode);

    //void setupNode();//BNodeにNodeをセット : 辺,面BNodeすべてにセット
    void setupNode_Edge();
    void setupNode_Face();


protected:
    vector<CBoundaryFace*> mvProgBFace;
    
    double mArea;//Faceの面積
    double  triArea(CNode* pNode0, CNode* pNode1, CNode* pNode2);// 外積による三角形面積(四辺形の場合は,2回使用)

    
public:
    void refine(uint& countID, const vuint& vDOF);// 自身の再分割
    
    double& calcArea();//Faceの面積計算(境界値再配分の為)
    double& getArea(){ return mArea;}

    vector<CBoundaryFace*>& getProgParts(){ return mvProgBFace;}

    void distDirichletVal(const uint& dof, const uint& mgLevel, const uint& nMaxMGLevel);//上位グリッドBNodeへのディレクレ値の分配

    void replaceEdgeBNode();//2次要素の場合、辺BNodeをmvBNodeに移設.

    void deleteProgData();//Refine時の辺BNode*の解放
};
#endif	/* _BOUNDARYFACE_H */
}



