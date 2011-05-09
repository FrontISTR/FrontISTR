/* 
 * File:   EdgeFaceTree.h
 * Author: ktakeda
 *
 * Created on 2009/06/30, 13:57
 */

#ifndef _EDGEFACETREE_H_1a9fde10_f2e1_4a72_88ac_e9b7696c73c4
#define	_EDGEFACETREE_H_1a9fde10_f2e1_4a72_88ac_e9b7696c73c4

#include "CommonStd.h"
#include "TypeDef.h"

#include "Logger.h"

namespace pmw{
class CEdgeFaceTree{
public:
    static CEdgeFaceTree* Instance(){
        static CEdgeFaceTree edgeFaceTree;
        return &edgeFaceTree;
    }
    virtual ~CEdgeFaceTree();
    
private:
    CEdgeFaceTree();
    // Face 番号
    uint mFaceIndex;

    // Faceを構成するEdge番号
    uint mHexaFaceConnEdge[6][4];
    uint mTetraFaceConnEdge[4][3];
    uint mPrismFaceConnEdge[5][4];
    uint mPyramidFaceConnEdge[5][4];
    uint mQuadFaceConnEdge[1][4];
    uint mTriangleFaceConnEdge[1][3];

public:
    uint& getHexaFaceIndex(const uint& edge0, const uint& edge1);
    uint& getTetraFaceIndex(const uint& edge0, const uint& edge1);
    uint& getPrismFaceIndex(const uint& edge0, const uint& edge1);
    uint& getPyramidFaceIndex(const uint& edge0, const uint& edge1);
    uint& getQuadFaceIndex(const uint& edge0, const uint& edge1);
    uint& getTriangleFaceIndex(const uint& edge0, const uint& edge1);

    // Faceを構成するEdge番号配列の提供
    uint* getHexaFaceConnEdge(const uint& iface){ return mHexaFaceConnEdge[iface];}
    uint* getTetraFaceConnEdge(const uint& iface){ return mTetraFaceConnEdge[iface];}
    uint* getPrismFaceConnEdge(const uint& iface){ return mPrismFaceConnEdge[iface];}
    uint* getPyramidFaceConnEdge(const uint& iface){ return mPyramidFaceConnEdge[iface];}
    uint* getQuadFaceConnEdge(const uint& iface){ return mQuadFaceConnEdge[iface];}
    uint* getTriangleFaceConnEdge(const uint& iface){ return mTriangleFaceConnEdge[iface];}
};
}

#endif	/* _EDGEFACETREE_H */

