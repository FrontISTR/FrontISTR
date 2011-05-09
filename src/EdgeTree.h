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
    uint mEdgeNum;
public:
    static CEdgeTree* Instance(){
        static CEdgeTree edgeTree;
        return &edgeTree;
    }
    virtual ~CEdgeTree();

public:
    uint& getHexaEdgeIndex(const uint& localNum0, const uint& localNum1);
    uint& getTetraEdgeIndex(const uint& localNum0, const uint& localNum1);
    uint& getPrismEdgeIndex(const uint& localNum0, const uint& localNum1);
    uint& getPyramidEdgeIndex(const uint& localNum0, const uint& localNum1);
    uint& getQuadEdgeIndex(const uint& localNum0, const uint& localNum1);
    uint& getTriangleEdgeIndex(const uint& localNum0, const uint& localNum1);
    uint& getBeamEdgeIndex(const uint& localNum0, const uint& localNum1);
};
}

#endif	/* _EDGETREE_H */

