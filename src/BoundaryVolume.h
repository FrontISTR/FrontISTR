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
    
    void markingEdge(const uint& iedge){ mvbMarkingEdge[iedge]=true;}
    void markingFace(const uint& iface){ mvbMarkingFace[iface]=true;}

    bool isMarkingEdge(const uint& iedge){ return mvbMarkingEdge[iedge];}
    bool isMarkingFace(const uint& iface){ return mvbMarkingFace[iface];}

    virtual uint getElemType()=0;
    virtual uint getNumOfEdge()=0;
    virtual uint getNumOfFace()=0;

    
    virtual PairBNode getPairBNode(const uint& iedge)=0;//辺の両端のBNode
    virtual uint& getEdgeID(PairBNode& pairBNode)=0;

    virtual vector<CBoundaryNode*> getFaceCnvNodes(const uint& iface)=0;//面を構成するBNode
    virtual uint& getFaceID(vector<CBoundaryNode*>& vBNode)=0;


    void setEdgeNeibVol(const uint& iedge, const uint& neibVolID);
    void setFaceNeibVol(const uint& iface, const uint& neibVolID);
    vuint& getEdgeNeibVolID(const uint& iedge){ return mvEdgeNeibVol[iedge];}
    uint&  getFaceNeibVolID(const uint& iface){ return mvFaceNeibVol[iface];}

    void setEdgeBNode(const uint& iedge, CBoundaryNode *pBNode);
    void setFaceBNode(const uint& iface, CBoundaryNode *pBNode);
    void setVolBNode(CBoundaryNode *pBNode);
    CBoundaryNode* getEdgeBNode(const uint& iedge){ return mvEdgeBNode[iedge];}
    CBoundaryNode* getFaceBNode(const uint& iface){ return mvFaceBNode[iface];}
    CBoundaryNode* getVolBNode(){ return mpVolBNode;}

protected:
    
    vector<CBoundaryVolume*> mvProgVolume;

    double mCubicVolume;// BoundaryVolumeの体積

    double tetraVolume(CNode* pNode0, CNode* pNode1, CNode* pNode2, CNode* pNode3);//四面体の体積


    uint dividHexa(const uint& iprog, CBoundaryVolume* pProgVol);//return は、子供がぶらさがる場所(BNodeの頂点番号)
    uint dividTetra(const uint& iprog, CBoundaryVolume* pProgVol);
    uint dividPrism(const uint& iprog, CBoundaryVolume* pProgVol);

    void distValue(CBoundaryVolume *pProgVol, const double& coef, const vuint& vDOF);//要素の境界値の分配

    virtual uint* getLocalNode_Edge(const uint& iedge)=0;//setupNode_Edgeで使用
    virtual uint* getLocalNode_Face(const uint& iface)=0;//setupNode_Faceで使用
public:
    // Refine 再分割
    virtual void refine(uint& countID, const vuint& vDOF)=0;

    vector<CBoundaryVolume*>& getProgParts(){ return mvProgVolume;}

    virtual double& calcVolume()=0;// BoundaryVolumeの体積計算(計算結果はmCubicVolumeへ)
    double& getCubicVolume(){ return mCubicVolume;};

    virtual void distDirichletVal(const uint& dof, const uint& mgLevel)=0;//上位グリッドBNodeへのディレクレ値の分配


    // Refine 後処理 : 辺-面 BNode vectorの解放
    // ----
    virtual void deleteProgData()=0;

};
#endif	/* _BOUNDARYELEMENT_H */
}


