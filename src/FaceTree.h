/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FaceTree.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _FACETREE_H_d4a90ba4_2f6f_42de_9b13_f78ae0961894
#define	_FACETREE_H_d4a90ba4_2f6f_42de_9b13_f78ae0961894
#include "CommonStd.h"
#include "TypeDef.h"
namespace pmw{
class CFaceTree{
private:
    CFaceTree();
    uint  mFaceIndex;
    uint mHexaFaceIndex[6][4];
    uint mTetraFaceIndex[4][3];
    uint mPrismFaceIndex[5][4];
    uint mPyramidFaceIndex[5][4];
    uint mQuadFaceIndex[4];
    uint mTriangleFaceIndex[3];
    uint mHexaConnFace[8][3];
    uint mTetraConnFace[4][3];
    uint mPrismConnFace[6][3];
    uint mPyramidConnFace[5][4];
    uint mQuadConnFace[4][1];
    uint mTriangleConnFace[3][1];
public:
    static CFaceTree* Instance(){
        static CFaceTree faceTree;
        return &faceTree;
    }
    virtual ~CFaceTree();
    uint& getHexaFaceIndex(const vuint& vLocalNodeIndex);
    uint& getHexaFaceIndex2(const vuint& vLocalNodeIndex);
    uint& getTetraFaceIndex(const vuint& vLocalNodeIndex);
    uint& getTetraFaceIndex2(const vuint& vLocalNodeIndex);
    uint& getPrismFaceIndex(const vuint& vLocalNodeIndex);
    uint& getPrismFaceIndex2(const vuint& vLocalNodeIndex);
    uint& getPyramidFaceIndex(const vuint& vLocalNodeIndex);
    uint& getPyramidFaceIndex2(const vuint& vLocalNodeIndex);
    uint& getQuadFaceIndex(const vuint& vLocalNodeIndex);
    uint& getTriangleFaceIndex(const vuint& vLocalNodeIndex);
    uint* getLocalNodeHexaFace(const uint& faceIndex);
    uint* getLocalNodeTetraFace(const uint& faceIndex);
    uint* getLocalNodePrismFace(const uint& faceIndex);
    uint* getLocalNodePyramidFace(const uint& faceIndex);
    uint* getLocalNodeQuadFace();
    uint* getLocalNodeTriangleFace();
    uint* getHexaConnFace(const uint& ivert){ return mHexaConnFace[ivert];}
    uint* getTetraConnFace(const uint& ivert){ return mTetraConnFace[ivert];}
    uint* getPrismConnFace(const uint& ivert){ return mPrismConnFace[ivert];}
    uint* getPyramidConnFace(const uint& ivert){ return mPyramidConnFace[ivert];}
    uint* getQuadConnFace(const uint& ivert){ return mQuadConnFace[ivert];}
    uint* getTriangleConnFace(const uint& ivert){ return mTriangleConnFace[ivert];}
};
}
#endif	/* _FACETREE_H */
