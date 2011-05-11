/* 
 * File:   BoundaryFaceMesh.h
 * Author: ktakeda
 *
 * BoundaryFace管理クラス
 *
 * Created on 2010/04/12, 16:59
 */
#include "TypeDef.h"
#include <map>
#include <utility>//pair

#include "BoundaryFace.h"//個別の面(面自体が値を持つ)

// Equivalent Nodal Force(等価節点力)
#include "ShapeQuad.h"
#include "ShapeTriangle.h"



#include "BoundaryMesh.h"

typedef std::pair<pmw::CBoundaryNode*, pmw::CBoundaryNode*> PairBNode;

namespace pmw{
#ifndef _BOUNDARYFACEMESH_H
#define	_BOUNDARYFACEMESH_H
class CBoundaryFaceMesh:public CBoundaryMesh{
public:
    CBoundaryFaceMesh();
    virtual ~CBoundaryFaceMesh();

protected:
    vector<CBoundaryFace*> mvBFace;
    map<uiint, uiint, less<uiint> > mmBFaceID2Index;
    
    vvuint mvAggregateFace;// BNodeの面-集合

    vector<CBoundaryNode*> mvBEdgeBNode;// 辺ノード
    vector<CBoundaryNode*> mvBFaceBNode;// 面ノード

public:   
    //境界面 :=> BFace
    void resizeBFace(const uiint& res_size){ mvBFace.resize(res_size);}
    void setBFace(const uiint& index, CBoundaryFace *pBFace);
    void addBFace(CBoundaryFace *pBFace);
    uiint getNumOfBFace(){ return mvBFace.size();}
    CBoundaryFace* getBFaceIX(const uiint& index){ return mvBFace[index];}
    CBoundaryFace* getBFaceID(const uiint& id){
        uiint index= mmBFaceID2Index[id];
        return mvBFace[index];
    }
    uiint& getBFaceIndex(const uiint& id){ return mmBFaceID2Index[id];}


public:
    // 面-集合
    // ----
    void resizeAggFace();
    vuint& getAggFace(const uiint& ibnode){ return mvAggregateFace[ibnode];}
    void setupAggFace();// 点-周囲の面集合

    void setAggFace(const uiint& ibnode, const uiint& nFaceID){ mvAggregateFace[ibnode].push_back(nFaceID);}


    // refine準備
    // ----
    void GeneEdgeBNode();// 辺ノード生成 { 辺-周囲の面集合処理を含む }
    void GeneFaceBNode();// 面ノード生成

//    uint getNumOfBEdge(){ return mvBEdgeBNode.size();}
//    CBoundaryNode* getEdgeBNode(const uint& iedge){ return mvBEdgeBNode[iedge];}// 辺ノード提供
//    CBoundaryNode* getFaceBNode(const uint& iface){ return mvBFaceBNode[iface];}// 面ノード提供


    // 境界要素の再分割 => progBMeshにセット
    // ----
    void refine(CBoundaryFaceMesh *pProgBFaceMesh);

    // Refine時のデータの解放
    // ----
    // 各Faceの辺BNodeの解放
    // mvBEdgeBNodeの解放
    // mvBFaceBNodeの解放
    void deleteProgData();
    
protected:
    //Neumann条件の節点分配(形状関数による等価節点力:EquivalentNodeForce)
    //----
    virtual void distNeumannValue();

////    //Dirichlet条件の節点分配(BNode周囲の面集合平均値)
////    //----
////    virtual void distDirichletValue_at_CGrid();//Level==0(coarse grid)
////    virtual void distDirichletValue_at_FGrid();//Level>=1(fine grid)
    virtual void distDirichletValue();
    
};
#endif	/* _BOUNDARYFACEMESH_H */
}


