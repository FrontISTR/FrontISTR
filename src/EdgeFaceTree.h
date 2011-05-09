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
    uint mFaceIndex;

public:
    uint& getHexaFaceIndex(const uint& edge0, const uint& edge1);
    uint& getTetraFaceIndex(const uint& edge0, const uint& edge1);
    uint& getPrismFaceIndex(const uint& edge0, const uint& edge1);
    uint& getPyramidFaceIndex(const uint& edge0, const uint& edge1);

    uint& getQuadFaceIndex(const uint& edge0, const uint& edge1);
    uint& getTriangleFaceIndex(const uint& edge0, const uint& edge1);
};
}

#endif	/* _EDGEFACETREE_H */

