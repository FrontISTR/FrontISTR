/* 
 * File:   BoundaryVolume.h
 * Author: ktakeda
 *
 *  重加
 *  速度力
 *  遠心力 回転軸(点1,点2),角速度
 *  発熱
 *
 * Modify     2009/05/18
 * Created on 2009/05/13, 16:08
 */
#include "EdgeTree.h"
#include "FaceTree.h"

#include "DiscreteVolume.h"//Hexa,PrismをTetraへ分割したときの頂点の並び

#include "BoundaryParts.h"

typedef std::pair<pmw::CBoundaryNode*, pmw::CBoundaryNode*> PairBNode;

namespace pmw{
#ifndef _BOUNDARYELEMENT_H_
#define	_BOUNDARYELEMENT_H_
class CBoundaryVolume:public CBoundaryParts{
public:
    CBoundaryVolume();// 形状( Tetra | Hexa | Prism )初期化
    virtual ~CBoundaryVolume();

protected:
    // Mesh-ElementのエンティティID
    // => 常に0

    bool* mvbMarkingEdge;
    bool* mvbMarkingFace;
    
    vector<CBoundaryNode*> mvEdgeBNode;// 辺-BNode
    vector<CBoundaryNode*> mvFaceBNode;// 面-BNode
    CBoundaryNode          *mpVolBNode;// 体-BNode

    vvuint mvEdgeNeibVol;// 辺に隣接するVolume(複数)
    vuint  mvFaceNeibVol;// 面に隣接するVolume(１つ)

public:
    
    void markingEdge(const uiint& iedge){ mvbMarkingEdge[iedge]=true;}
    void markingFace(const uiint& iface){ mvbMarkingFace[iface]=true;}

    bool isMarkingEdge(const uiint& iedge){ return mvbMarkingEdge[iedge];}
    bool isMarkingFace(const uiint& iface){ return mvbMarkingFace[iface];}

    virtual uiint getElemType()=0;
    virtual uiint getNumOfEdge()=0;
    virtual uiint getNumOfFace()=0;
    virtual uiint getNumOfNode()=0;
    //Vertex ==> BoundaryParts::getNumOfVert()=0

    virtual void setOrder(const uiint& order)=0;//{ mnOrder=order;};//Hexa,Tetra,Prismの1次・2次
    
    virtual PairBNode getPairBNode(const uiint& iedge)=0;//辺の両端のBNode
    virtual uiint& getEdgeID(PairBNode& pairBNode)=0;

    virtual vector<CBoundaryNode*> getFaceCnvNodes(const uiint& iface)=0;//面を構成するBNode
    virtual uiint& getFaceID(vector<CBoundaryNode*>& vBNode)=0;


    void setEdgeNeibVol(const uiint& iedge, const uiint& neibVolID);
    void setFaceNeibVol(const uiint& iface, const uiint& neibVolID);
    vuint& getEdgeNeibVolID(const uiint& iedge){ return mvEdgeNeibVol[iedge];}
    uiint&  getFaceNeibVolID(const uiint& iface){ return mvFaceNeibVol[iface];}

    void setEdgeBNode(const uiint& iedge, CBoundaryNode *pBNode);
    void setFaceBNode(const uiint& iface, CBoundaryNode *pBNode);
    void setVolBNode(CBoundaryNode *pBNode);
    CBoundaryNode* getEdgeBNode(const uiint& iedge){ return mvEdgeBNode[iedge];}
    CBoundaryNode* getFaceBNode(const uiint& iface){ return mvFaceBNode[iface];}
    CBoundaryNode* getVolBNode(){ return mpVolBNode;}

protected:
    
    vector<CBoundaryVolume*> mvProgVolume;

    double mCubicVolume;// BoundaryVolumeの体積

    double tetraVolume(CNode* pNode0, CNode* pNode1, CNode* pNode2, CNode* pNode3);//四面体の体積


    uiint dividHexa(const uiint& iprog, CBoundaryVolume* pProgVol);//return は、子供がぶらさがる場所(BNodeの頂点番号)
    uiint dividTetra(const uiint& iprog, CBoundaryVolume* pProgVol);
    uiint dividPrism(const uiint& iprog, CBoundaryVolume* pProgVol);

    void distValue(CBoundaryVolume *pProgVol, const double& coef, const vuint& vDOF);//要素の境界値の分配

    virtual uiint* getLocalNode_Edge(const uiint& iedge)=0;//setupNode_Edgeで使用
    virtual uiint* getLocalNode_Face(const uiint& iface)=0;//setupNode_Faceで使用
public:
    // Refine 再分割
    virtual void refine(uiint& countID, const vuint& vDOF)=0;

    vector<CBoundaryVolume*>& getProgParts(){ return mvProgVolume;}

    virtual double& calcVolume()=0;// BoundaryVolumeの体積計算(計算結果はmCubicVolumeへ)
    double& getCubicVolume(){ return mCubicVolume;};

    virtual void distDirichletVal(const uiint& dof, const uiint& mgLevel, const uiint& nMaxMGLevel)=0;//上位グリッドBNodeへのディレクレ値の分配

    virtual void replaceEdgeBNode(const uiint& iedge)=0;//2次要素の場合に辺BNodeをmvBNodeへ移設.

    // Refine 後処理 : 辺-面 BNode vectorの解放
    // ----
    virtual void deleteProgData()=0;

};
#endif	/* _BOUNDARYELEMENT_H */
}


