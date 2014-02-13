/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/EdgeFaceTree.h
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
#ifndef _EDGEFACETREE_H_1a9fde10_f2e1_4a72_88ac_e9b7696c73c4
#define	_EDGEFACETREE_H_1a9fde10_f2e1_4a72_88ac_e9b7696c73c4
#include "CommonStd.h"
#include "TypeDef.h"
#include "Logger.h"
namespace pmw
{
class CEdgeFaceTree
{
public:
    static CEdgeFaceTree* Instance() {
        static CEdgeFaceTree edgeFaceTree;
        return &edgeFaceTree;
    }
    virtual ~CEdgeFaceTree();
private:
    CEdgeFaceTree();
    uiint mFaceIndex;
    uiint mHexaFaceConnEdge[6][4];
    uiint mTetraFaceConnEdge[4][3];
    uiint mPrismFaceConnEdge[5][4];
    uiint mPyramidFaceConnEdge[5][4];
    uiint mQuadFaceConnEdge[1][4];
    uiint mTriangleFaceConnEdge[1][3];
public:
    uiint& getHexaFaceIndex(const uiint& edge0, const uiint& edge1);
    uiint& getTetraFaceIndex(const uiint& edge0, const uiint& edge1);
    uiint& getPrismFaceIndex(const uiint& edge0, const uiint& edge1);
    uiint& getPyramidFaceIndex(const uiint& edge0, const uiint& edge1);
    uiint& getQuadFaceIndex(const uiint& edge0, const uiint& edge1);
    uiint& getTriangleFaceIndex(const uiint& edge0, const uiint& edge1);
    uiint* getHexaFaceConnEdge(const uiint& iface) {
        return mHexaFaceConnEdge[iface];
    }
    uiint* getTetraFaceConnEdge(const uiint& iface) {
        return mTetraFaceConnEdge[iface];
    }
    uiint* getPrismFaceConnEdge(const uiint& iface) {
        return mPrismFaceConnEdge[iface];
    }
    uiint* getPyramidFaceConnEdge(const uiint& iface) {
        return mPyramidFaceConnEdge[iface];
    }
    uiint* getQuadFaceConnEdge(const uiint& iface) {
        return mQuadFaceConnEdge[iface];
    }
    uiint* getTriangleFaceConnEdge(const uiint& iface) {
        return mTriangleFaceConnEdge[iface];
    }
};
}
#endif	/* _EDGEFACETREE_H */
