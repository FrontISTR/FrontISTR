//
//  Beam : BeamBase : Element
//
//				2009.06.19
//				2008.12.01
//				k.Takeda
#ifndef BEAM_HH_785972AE_015C_439f_AA7A_0464A768D304
#define BEAM_HH_785972AE_015C_439f_AA7A_0464A768D304


//#include "Element.h"
#include "BeamBase.h"

namespace pmw{
class CBeam:public CBeamBase
{
public:
    CBeam(void);
    virtual ~CBeam(void);

private:
    static uiint mnElemType;
    static uiint mnElemOrder;
    static uiint mNumOfFace;
    static uiint mNumOfEdge;
    static uiint mNumOfNode;
    static uiint mNumOfVert;

    ////API 形状関数・導関数
    //pmw::CShapeLine     *mpShapeLine;

protected:
    virtual bool IndexCheck(const uiint& propType, const uiint& index, string& method_name);
    virtual uiint& getLocalFaceNum(const vuint& vLocalNodeNum);
    //    virtual uint& getLocalFaceNum(CNode* pNode0, CNode* pNode1, CNode* pNode2);

public:
    virtual void initialize();

public:
    // Property
    virtual const uiint& getType();
    virtual const uiint& getOrder();
    virtual const uiint& getNumOfFace(){ return mNumOfFace;}
    virtual const uiint& getNumOfEdge(){ return mNumOfEdge;}
    virtual const uiint& getNumOfNode(){ return mNumOfNode;}
    virtual const uiint& getNumOfVert(){ return mNumOfVert;}


    // EdgeNode(Pair Node)
    virtual PairNode getPairNode(const uiint& iedge);
    virtual void getPairNode(vint& pairNodeIndex, const uiint& iedge);
    virtual bool isEdgeElem(CNode* pNode0, CNode* pNode1);
    virtual void setBoolEdgeElem(CNode* pNode0, CNode* pNode1);
    // ノードの局所番号に対応した、辺(Edge)の番号
    virtual uiint& getEdgeIndex(CNode* pNode0, CNode* pNode1);
    virtual uiint& getEdgeIndex(const uiint& nodeID_0, const uiint& nodeID_1);

    // Edge -> Face Index
    virtual uiint& getFaceIndex(const uiint& edge0, const uiint& edge1);
    // 3 Nodes -> Face Index
    virtual uiint& getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2);

    // Faceノード
    virtual bool isFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2);   //Faceにノードがセットされているか？
    virtual void setBoolFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2);//Faceノードがセットされたことをスタンプ

    // Face構成Node
    //
    virtual vector<CNode*> getFaceCnvNodes(const uiint& iface);


    // Node周囲の接続先Nodeを配列で返す. (係数行列 作成用途, CMesh::setupAggregate)
    //
    virtual vector<CNode*> getConnectNode(CNode* pNode);


    // Refine後
    // 1. 辺-面 Element*配列を解放
    // 2. 辺-面 Node*   配列を解放 (2次要素は辺ノードを残す)
    // --
    virtual void deleteProgData();
};
}
#endif
