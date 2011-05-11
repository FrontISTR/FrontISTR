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

#include "Logger.h"

namespace pmw{
class CEdgeTree{
private:
    CEdgeTree();
    uiint mEdgeNum;  //Edge番号
    //vuint mvNodeNum;//Edge番号の端点の局所ノード番号

    Utility::CLogger *mpLogger;

    // Edgeを構成する局所ノード番号
    uiint mHexaEdgeIndex[12][2];
    uiint mTetraEdgeIndex[6][2];
    uiint mPrismEdgeIndex[9][2];
    uiint mPyramidEdgeIndex[8][2];
    uiint mQuadEdgeIndex[4][2];
    uiint mTriangleEdgeIndex[3][2];
    uiint mBeamEdgeIndex[1][2];

    //Nodeが接続するEdge番号
    uiint mHexaConnEdge[8][3];
    uiint mTetraConnEdge[4][3];
    uiint mPrismConnEdge[6][3];
    uiint mPyramidConnEdge[5][4];
    uiint mQuadConnEdge[4][2];
    uiint mTriangleConnEdge[3][2];
    uiint mBeamConnEdge[2][1];

public:
    static CEdgeTree* Instance(){
        static CEdgeTree edgeTree;
        return &edgeTree;
    }
    virtual ~CEdgeTree();

public:
    // 辺の端点の局所ノード番号から,辺(Edge)番号を提供
    uiint& getHexaEdgeIndex(const uiint& localNum0, const uiint& localNum1);
    uiint& getTetraEdgeIndex(const uiint& localNum0, const uiint& localNum1);
    uiint& getPrismEdgeIndex(const uiint& localNum0, const uiint& localNum1);
    uiint& getPyramidEdgeIndex(const uiint& localNum0, const uiint& localNum1);
    uiint& getQuadEdgeIndex(const uiint& localNum0, const uiint& localNum1);
    uiint& getTriangleEdgeIndex(const uiint& localNum0, const uiint& localNum1);
    uiint& getBeamEdgeIndex(const uiint& localNum0, const uiint& localNum1);

    uiint& getDisagTypeEdgeIndex(const uiint& localNum0, const uiint& localNum1);
    
    // 辺(Edge)番号から,端点の局所ノード番号を提供
    uiint* getHexaLocalNodeNum(const uiint& edgeNum){ return mHexaEdgeIndex[edgeNum];}
    uiint* getTetraLocalNodeNum(const uiint& edgeNum){ return mTetraEdgeIndex[edgeNum];}
    uiint* getPrismLocalNodeNum(const uiint& edgeNum){ return mPrismEdgeIndex[edgeNum];}
    uiint* getPyramidLocalNodeNum(const uiint& edgeNum){ return mPyramidEdgeIndex[edgeNum];}
    uiint* getQuadLocalNodeNum(const uiint& edgeNum){ return mQuadEdgeIndex[edgeNum];}
    uiint* getTriangleLocalNodeNum(const uiint& edgeNum){ return mTriangleEdgeIndex[edgeNum];}
    uiint* getBeamLocalNodeNum(const uiint& edgeNum){ return mBeamEdgeIndex[edgeNum];}
    // 統合版
    uiint* getLocalNodeNum(const uiint& elemType, const uiint& iedge);

    
    // 一点に接続する,辺の番号を提供
    uiint* getHexaConnEdge(const uiint& ivert){ return mHexaConnEdge[ivert];}
    uiint* getTetraConnEdge(const uiint& ivert){ return mTetraConnEdge[ivert];}
    uiint* getPrismConnEdge(const uiint& ivert){ return mPrismConnEdge[ivert];}
    uiint* getPyramidConnEdge(const uiint& ivert){ return mPyramidConnEdge[ivert];}
    uiint* getQuadConnEdge(const uiint& ivert){ return mQuadConnEdge[ivert];}
    uiint* getTriangleConnEdge(const uiint& ivert){ return mTriangleConnEdge[ivert];}
    uiint* getBeamConnEdge(const uiint& ivert){ return mBeamConnEdge[ivert];}
};
}

#endif	/* _EDGETREE_H */

