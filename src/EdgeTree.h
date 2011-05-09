/* 
 * File:   EdgeTree.h
 * Author: ktakeda
 *
 * 要素辺の端点の局所ノード番号から、辺(Edge)番号を取得
 * (各種の要素形状を一まとめに関数として記述）
 *
 * Created on 2009/06/23, 15:55
 */

#ifndef _EDGETREE_H_32180645_1eef_4b49_8624_c3bf4ed245a1
#define	_EDGETREE_H_32180645_1eef_4b49_8624_c3bf4ed245a1

#include "CommonStd.h"
#include "TypeDef.h"

namespace pmw{
class CEdgeTree{
private:
    CEdgeTree();
    uint mEdgeNum;  //Edge番号
    //vuint mvNodeNum;//Edge番号の端点の局所ノード番号

    // Edgeを構成する局所ノード番号
    uint mHexaEdgeIndex[12][2];
    uint mTetraEdgeIndex[6][2];
    uint mPrismEdgeIndex[9][2];
    uint mPyramidEdgeIndex[8][2];
    uint mQuadEdgeIndex[4][2];
    uint mTriangleEdgeIndex[3][2];
    uint mBeamEdgeIndex[1][2];

    //Nodeが接続するEdge番号
    uint mHexaConnEdge[8][3];
    uint mTetraConnEdge[4][3];
    uint mPrismConnEdge[6][3];
    uint mPyramidConnEdge[5][4];
    uint mQuadConnEdge[4][2];
    uint mTriangleConnEdge[3][2];
    uint mBeamConnEdge[2][1];

public:
    static CEdgeTree* Instance(){
        static CEdgeTree edgeTree;
        return &edgeTree;
    }
    virtual ~CEdgeTree();

public:
    // 辺の端点の局所ノード番号から,辺(Edge)番号を提供
    uint& getHexaEdgeIndex(const uint& localNum0, const uint& localNum1);
    uint& getTetraEdgeIndex(const uint& localNum0, const uint& localNum1);
    uint& getPrismEdgeIndex(const uint& localNum0, const uint& localNum1);
    uint& getPyramidEdgeIndex(const uint& localNum0, const uint& localNum1);
    uint& getQuadEdgeIndex(const uint& localNum0, const uint& localNum1);
    uint& getTriangleEdgeIndex(const uint& localNum0, const uint& localNum1);
    uint& getBeamEdgeIndex(const uint& localNum0, const uint& localNum1);
    uint& getDisagTypeEdgeIndex(const uint& localNum0, const uint& localNum1);
    
    // 辺(Edge)番号から,端点の局所ノード番号を提供
    uint* getHexaLocalNodeNum(const uint& edgeNum){ return mHexaEdgeIndex[edgeNum];}
    uint* getTetraLocalNodeNum(const uint& edgeNum){ return mTetraEdgeIndex[edgeNum];}
    uint* getPrismLocalNodeNum(const uint& edgeNum){ return mPrismEdgeIndex[edgeNum];}
    uint* getPyramidLocalNodeNum(const uint& edgeNum){ return mPyramidEdgeIndex[edgeNum];}
    uint* getQuadLocalNodeNum(const uint& edgeNum){ return mQuadEdgeIndex[edgeNum];}
    uint* getTriangleLocalNodeNum(const uint& edgeNum){ return mTriangleEdgeIndex[edgeNum];}
    uint* getBeamLocalNodeNum(const uint& edgeNum){ return mBeamEdgeIndex[edgeNum];}
    
    // 一点に接続する,辺の番号を提供
    uint* getHexaConnEdge(const uint& ivert){ return mHexaConnEdge[ivert];}
    uint* getTetraConnEdge(const uint& ivert){ return mTetraConnEdge[ivert];}
    uint* getPrismConnEdge(const uint& ivert){ return mPrismConnEdge[ivert];}
    uint* getPyramidConnEdge(const uint& ivert){ return mPyramidConnEdge[ivert];}
    uint* getQuadConnEdge(const uint& ivert){ return mQuadConnEdge[ivert];}
    uint* getTriangleConnEdge(const uint& ivert){ return mTriangleConnEdge[ivert];}
    uint* getBeamConnEdge(const uint& ivert){ return mBeamConnEdge[ivert];}
};
}

#endif	/* _EDGETREE_H */

