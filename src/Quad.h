//
//  Quad.h 
//  Quad : ShellBase : Element
//
//				2008.06.19
//				2008.12.01
//				k.Takeda
#ifndef QUAD_HH_6BF7327A_6740_43e8_ADA1_B90414CCA216
#define QUAD_HH_6BF7327A_6740_43e8_ADA1_B90414CCA216

#include "ShellBase.h"
//#include "Element.h"

namespace pmw{
class CQuad:public CShellBase
{
public:
    CQuad(void);
    virtual ~CQuad(void);

private:
    static uint mnElemType;
    static uint mnElemType2;//2次要素
    static uint mNumOfFace;
    static uint mNumOfEdge;
    static uint mNumOfNode;

    uint mTempo;

    ////API 形状関数・導関数
    //pmw::CShapeQuad     *mpShapeQuad;

protected:
    virtual bool IndexCheck(const uint& propType, const uint& index, string& method_name);
    virtual uint& getLocalFaceNum(const vuint& vLocalNodeNum);//局所ノード番号が正しければ"0"
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
    // (局所ノード番号、局所ノード番号)に対応した辺(Edge)番号
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
    virtual bool isFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2);
    virtual void setBoolFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2);


    // Faceノード
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
