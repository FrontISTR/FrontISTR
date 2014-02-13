/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/EdgeTree.h
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _EDGETREE_H_32180645_1eef_4b49_8624_c3bf4ed245a1
#define	_EDGETREE_H_32180645_1eef_4b49_8624_c3bf4ed245a1
#include "CommonStd.h"
#include "TypeDef.h"
#include "Logger.h"
namespace pmw
{
class CEdgeTree
{
private:
    CEdgeTree();
    uiint mEdgeNum;
    Utility::CLogger *mpLogger;

    uiint mHexaEdgeIndex[12][2];
    uiint mTetraEdgeIndex[6][2];
    uiint mPrismEdgeIndex[9][2];
    uiint mPyramidEdgeIndex[8][2];
    uiint mQuadEdgeIndex[4][2];
    uiint mTriangleEdgeIndex[3][2];
    uiint mBeamEdgeIndex[1][2];
    uiint mHexaConnEdge[8][3];
    uiint mTetraConnEdge[4][3];
    uiint mPrismConnEdge[6][3];
    uiint mPyramidConnEdge[5][4];
    uiint mQuadConnEdge[4][2];
    uiint mTriangleConnEdge[3][2];
    uiint mBeamConnEdge[2][1];

    uiint mDummyEdgeIndex;
public:
    static CEdgeTree* Instance() {
        static CEdgeTree edgeTree;
        return &edgeTree;
    }
    virtual ~CEdgeTree();
public:
    uiint& getDummyEdgeIndex();

    uiint& getHexaEdgeIndex(const uiint& localNum0, const uiint& localNum1);
    uiint& getTetraEdgeIndex(const uiint& localNum0, const uiint& localNum1);
    uiint& getPrismEdgeIndex(const uiint& localNum0, const uiint& localNum1);
    uiint& getPyramidEdgeIndex(const uiint& localNum0, const uiint& localNum1);
    uiint& getQuadEdgeIndex(const uiint& localNum0, const uiint& localNum1);
    uiint& getTriangleEdgeIndex(const uiint& localNum0, const uiint& localNum1);
    uiint& getBeamEdgeIndex(const uiint& localNum0, const uiint& localNum1);
    uiint& getDisagTypeEdgeIndex(const uiint& localNum0, const uiint& localNum1);
    uiint* getHexaLocalNodeNum(const uiint& edgeNum) {
        return mHexaEdgeIndex[edgeNum];
    }
    uiint* getTetraLocalNodeNum(const uiint& edgeNum) {
        return mTetraEdgeIndex[edgeNum];
    }
    uiint* getPrismLocalNodeNum(const uiint& edgeNum) {
        return mPrismEdgeIndex[edgeNum];
    }
    uiint* getPyramidLocalNodeNum(const uiint& edgeNum) {
        return mPyramidEdgeIndex[edgeNum];
    }
    uiint* getQuadLocalNodeNum(const uiint& edgeNum) {
        return mQuadEdgeIndex[edgeNum];
    }
    uiint* getTriangleLocalNodeNum(const uiint& edgeNum) {
        return mTriangleEdgeIndex[edgeNum];
    }
    uiint* getBeamLocalNodeNum(const uiint& edgeNum) {
        return mBeamEdgeIndex[edgeNum];
    }
    uiint* getLocalNodeNum(const uiint& elemType, const uiint& iedge);
    uiint* getHexaConnEdge(const uiint& ivert) {
        return mHexaConnEdge[ivert];
    }
    uiint* getTetraConnEdge(const uiint& ivert) {
        return mTetraConnEdge[ivert];
    }
    uiint* getPrismConnEdge(const uiint& ivert) {
        return mPrismConnEdge[ivert];
    }
    uiint* getPyramidConnEdge(const uiint& ivert) {
        return mPyramidConnEdge[ivert];
    }
    uiint* getQuadConnEdge(const uiint& ivert) {
        return mQuadConnEdge[ivert];
    }
    uiint* getTriangleConnEdge(const uiint& ivert) {
        return mTriangleConnEdge[ivert];
    }
    uiint* getBeamConnEdge(const uiint& ivert) {
        return mBeamConnEdge[ivert];
    }
};
}
#endif	/* _EDGETREE_H */
