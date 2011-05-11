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
    
    uiint  mFaceIndex;

    // Faceを構成する頂点(局所ノード番号)
    uiint mHexaFaceIndex[6][4];
    uiint mTetraFaceIndex[4][3];
    uiint mPrismFaceIndex[5][4];
    uiint mPyramidFaceIndex[5][4];
    uiint mQuadFaceIndex[4];
    uiint mTriangleFaceIndex[3];
    
    //uint QuadVertSize,TriVertSize;//Faceの頂点数(四辺形=4,三角形=3)
    //uint HexFaceSize, TetFaceSize;//Face数(Hexa=6, Tetra=4)

    // 頂点に接続する面番号
    uiint mHexaConnFace[8][3];
    uiint mTetraConnFace[4][3];
    uiint mPrismConnFace[6][3];
    uiint mPyramidConnFace[5][4];
    uiint mQuadConnFace[4][1];
    uiint mTriangleConnFace[3][1];


    // 面番号ごとの頂点数
    uiint mHexaFaceNumOfVert[6];
    uiint mTetraFaceNumOfVert[4];
    uiint mPrismFaceNumOfVert[5];
    uiint mPyramidFaceNumOfVert[5];
    uiint mQuadFaceNumOfVert[1];
    uiint mTriangleFaceNumOfVert[1];


public:
    static CFaceTree* Instance(){
        static CFaceTree faceTree;
        return &faceTree;
    }
    virtual ~CFaceTree();

    uiint& getHexaFaceIndex(const vuint& vLocalNodeIndex);
    uiint& getHexaFaceIndex2(const vuint& vLocalNodeIndex);

    uiint& getTetraFaceIndex(const vuint& vLocalNodeIndex);
    uiint& getTetraFaceIndex2(const vuint& vLocalNodeIndex);

    uiint& getPrismFaceIndex(const vuint& vLocalNodeIndex);
    uiint& getPrismFaceIndex2(const vuint& vLocalNodeIndex);

    uiint& getPyramidFaceIndex(const vuint& vLocalNodeIndex);
    uiint& getPyramidFaceIndex2(const vuint& vLocalNodeIndex);

    uiint& getQuadFaceIndex(const vuint& vLocalNodeIndex);//ノードが一致すれば,面番号=”0”を返す.
    uiint& getTriangleFaceIndex(const vuint& vLocalNodeIndex);//ノードが一致すれば,面番号=”0”を返す.

    uiint& getDefaultFaceIndex();//switch文のdefault用途

    // Face番号 -> 局所ノード番号
    uiint* getLocalNodeHexaFace(const uiint& faceIndex);
    uiint* getLocalNodeTetraFace(const uiint& faceIndex);
    uiint* getLocalNodePrismFace(const uiint& faceIndex);
    uiint* getLocalNodePyramidFace(const uiint& faceIndex);
    uiint* getLocalNodeQuadFace();
    uiint* getLocalNodeTriangleFace();
    // Face番号 -> 局所ノード番号(引数:要素タイプ)
    uiint* getLocalNodeFace(const uiint& elemType, const uiint& iface);

    // 局所ノード(Vertex)番号に接続する面(Face)番号
    uiint* getHexaConnFace(const uiint& ivert){ return mHexaConnFace[ivert];}
    uiint* getTetraConnFace(const uiint& ivert){ return mTetraConnFace[ivert];}
    uiint* getPrismConnFace(const uiint& ivert){ return mPrismConnFace[ivert];}
    uiint* getPyramidConnFace(const uiint& ivert){ return mPyramidConnFace[ivert];}
    uiint* getQuadConnFace(const uiint& ivert){ return mQuadConnFace[ivert];}
    uiint* getTriangleConnFace(const uiint& ivert){ return mTriangleConnFace[ivert];}

    //面番号毎の頂点数
    uiint& getHexaFaceNumOfVert(const uiint& iface){ return mHexaFaceNumOfVert[iface];}
    uiint& getTetraFaceNumOfVert(const uiint& iface){ return mTetraFaceNumOfVert[iface];}
    uiint& getPrismFaceNumOfVert(const uiint& iface){ return mPrismFaceNumOfVert[iface];}
    uiint& getPyramidFaceNumOfVert(const uiint& iface){ return mPyramidFaceNumOfVert[iface];}
    uiint& getQuadFaceNumOfVert(const uiint& iface){ return mQuadFaceNumOfVert[iface];}
    uiint& getTriangleFaceNumOfVert(const uiint& iface){ return mTriangleFaceNumOfVert[iface];}
};
}

#endif	/* _FACETREE_H */

