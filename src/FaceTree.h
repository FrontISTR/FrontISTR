/* 
 * File:   FaceTree.h
 * Author: ktakeda
 *
 * 4面体,6面体のVertex Index -> Face番号
 * (局所ノード番号から、Face番号を提供)
 *
 * Created on 2009/06/25, 17:59
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
    
    //uint QuadVertSize,TriVertSize;//Faceの頂点数(四辺形=4,三角形=3)
    //uint HexFaceSize, TetFaceSize;//Face数(Hexa=6, Tetra=4)
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

    uint& getQuadFaceIndex(const vuint& vLocalNodeIndex);//ノードが一致すれば,面番号=”0”を返す.
    uint& getTriangleFaceIndex(const vuint& vLocalNodeIndex);//ノードが一致すれば,面番号=”0”を返す.

    // 局所Face番号 -> 局所ノード番号
    uint* getLocalNodeHexaFace(const uint& faceIndex);
    uint* getLocalNodeTetraFace(const uint& faceIndex);
    uint* getLocalNodePrismFace(const uint& faceIndex);
    uint* getLocalNodePyramidFace(const uint& faceIndex);
    uint* getLocalNodeQuadFace();
    uint* getLocalNodeTriangleFace();
};
}

#endif	/* _FACETREE_H */

