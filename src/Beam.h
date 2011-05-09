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
    static uint mnElemType;
    static uint mnElemType2;//2次要素
    static uint mNumOfFace;
    static uint mNumOfEdge;
    static uint mNumOfNode;

    ////API 形状関数・導関数
    //pmw::CShapeLine     *mpShapeLine;

protected:
    virtual bool IndexCheck(const uint& propType, const uint& index, string& method_name);
    virtual uint& getLocalFaceNum(const vuint& vLocalNodeNum);
    //    virtual uint& getLocalFaceNum(CNode* pNode0, CNode* pNode1, CNode* pNode2);

public:
    // Property
    virtual const uint& getType();
    virtual const uint& getNumOfFace(){ return mNumOfFace;}
    virtual const uint& getNumOfEdge(){ return mNumOfEdge;}
    virtual const uint& getNumOfNode(){ return mNumOfNode;}

    // EdgeNode(Pair Node)
    virtual PairNode& getPairNode(const uint& edgeIndex);
    virtual void getPairNode(vint& pairNodeIndex, const uint& edgeIndex);
    virtual bool isEdgeElem(CNode* pNode0, CNode* pNode1);
    virtual void setBoolEdgeElem(CNode* pNode0, CNode* pNode1);
    // ノードの局所番号に対応した、辺(Edge)の番号
    virtual uint& getEdgeIndex(CNode* pNode0, CNode* pNode1);
    virtual uint& getEdgeIndex(const uint& nodeID_0, const uint& nodeID_1);

    // Edge -> Face Index
    virtual uint& getFaceIndex(const uint& edge0, const uint& edge1);
    // 3 Nodes -> Face Index
    virtual uint& getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2);

    // Face Element, Node
    virtual void setFaceElement(CElement* pElem, const uint& faceIndex);
    virtual void setFaceNode(CNode* pNode, const uint& faceIndex);
    virtual CNode* getFaceNode(const uint& faceIndex);


    // Faceノード
    virtual bool isFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2);   //Faceにノードがセットされているか？
    virtual void setBoolFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2);//Faceノードがセットされたことをスタンプ

    // mvvFaceCnvNodeのセットアップ
    virtual void setupFaceCnvNodes();

    // Node周囲の接続先Nodeを配列で返す.
    // (係数行列 作成用途, CMesh::setupAggregate)
    virtual vector<CNode*> getConnectNode(CNode* pNode);


    // prolongation後の後始末
    // --
    virtual void deleteProgData();
};
}
#endif
