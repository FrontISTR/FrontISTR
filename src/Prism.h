/* 
 * File:   Prism.h
 * Author: ktakeda
 *
 * Created on 2009/06/29, 14:19
 */

#ifndef _PRISM_H_4a9dd9b8_e22d_43b7_b22e_67c0612b96ba
#define	_PRISM_H_4a9dd9b8_e22d_43b7_b22e_67c0612b96ba

#include "SolidBase.h"

namespace pmw{
class CPrism:public CSolidBase
{
public:
    CPrism();
    virtual ~CPrism();

private:
    static uint mnElemType;
    static uint mnElemType2;//2次要素
    static uint mNumOfFace;
    static uint mNumOfEdge;
    static uint mNumOfNode;

    ////API 形状関数・導関数
    //pmw::CShapePrism    *mpShapePrism;

protected:
    virtual bool IndexCheck(const uint& propType, const uint& index, string& method_name);
    virtual uint& getLocalFaceNum(const vuint& vLocalNodeNum);
    //virtual uint& getLocalFaceNum(CNode* pNode0, CNode* pNode1, CNode* pNode2);

public:
    // Property
    virtual const uint& getType();
    virtual const uint& getNumOfFace(){ return mNumOfFace;}
    virtual const uint& getNumOfEdge(){ return mNumOfEdge;}
    virtual const uint& getNumOfNode(){ return mNumOfNode;}

    // EdgeNode(Pair Node)
    // --
    virtual PairNode& getPairNode(const uint& edgeIndex);
    virtual void getPairNode(vint& pairNodeIndex, const uint& edgeIndex);
    virtual bool isEdgeElem(CNode* pNode0, CNode* pNode1);
    virtual void setBoolEdgeElem(CNode* pNode0, CNode* pNode1);
    //
    // (局所ノード番号,局所ノード番号)に対応した辺の、Index番号
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

#endif	/* _PRISM_H */

