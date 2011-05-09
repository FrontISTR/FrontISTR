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

    // Faceを構成する頂点(局所ノード番号)
    uint mHexaFaceIndex[6][4];
    uint mTetraFaceIndex[4][3];
    uint mPrismFaceIndex[5][4];
    uint mPyramidFaceIndex[5][4];
    uint mQuadFaceIndex[4];
    uint mTriangleFaceIndex[3];
    
    //uint QuadVertSize,TriVertSize;//Faceの頂点数(四辺形=4,三角形=3)
    //uint HexFaceSize, TetFaceSize;//Face数(Hexa=6, Tetra=4)

    // 頂点に接続する面番号
    uint mHexaConnFace[8][3];
    uint mTetraConnFace[4][3];
    uint mPrismConnFace[6][3];
    uint mPyramidConnFace[5][4];
    uint mQuadConnFace[4][1];
    uint mTriangleConnFace[3][1];


    // 面番号ごとの頂点数
    uint mHexaFaceNumOfVert[6];
    uint mTetraFaceNumOfVert[4];
    uint mPrismFaceNumOfVert[5];
    uint mPyramidFaceNumOfVert[5];
    uint mQuadFaceNumOfVert[1];
    uint mTriangleFaceNumOfVert[1];


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

    uint& getDefaultFaceIndex();//switch文のdefault用途

    // Face番号 -> 局所ノード番号
    uint* getLocalNodeHexaFace(const uint& faceIndex);
    uint* getLocalNodeTetraFace(const uint& faceIndex);
    uint* getLocalNodePrismFace(const uint& faceIndex);
    uint* getLocalNodePyramidFace(const uint& faceIndex);
    uint* getLocalNodeQuadFace();
    uint* getLocalNodeTriangleFace();
    // Face番号 -> 局所ノード番号(引数:要素タイプ)
    uint* getLocalNodeFace(const uint& elemType, const uint& iface);

    // 局所ノード(Vertex)番号に接続する面(Face)番号
    uint* getHexaConnFace(const uint& ivert){ return mHexaConnFace[ivert];}
    uint* getTetraConnFace(const uint& ivert){ return mTetraConnFace[ivert];}
    uint* getPrismConnFace(const uint& ivert){ return mPrismConnFace[ivert];}
    uint* getPyramidConnFace(const uint& ivert){ return mPyramidConnFace[ivert];}
    uint* getQuadConnFace(const uint& ivert){ return mQuadConnFace[ivert];}
    uint* getTriangleConnFace(const uint& ivert){ return mTriangleConnFace[ivert];}

    //面番号毎の頂点数
    uint& getHexaFaceNumOfVert(const uint& iface){ return mHexaFaceNumOfVert[iface];}
    uint& getTetraFaceNumOfVert(const uint& iface){ return mTetraFaceNumOfVert[iface];}
    uint& getPrismFaceNumOfVert(const uint& iface){ return mPrismFaceNumOfVert[iface];}
    uint& getPyramidFaceNumOfVert(const uint& iface){ return mPyramidFaceNumOfVert[iface];}
    uint& getQuadFaceNumOfVert(const uint& iface){ return mQuadFaceNumOfVert[iface];}
    uint& getTriangleFaceNumOfVert(const uint& iface){ return mTriangleFaceNumOfVert[iface];}
};
}

#endif	/* _FACETREE_H */

