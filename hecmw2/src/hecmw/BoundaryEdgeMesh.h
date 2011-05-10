/* 
 * File:   BoundaryEdgeMesh.h
 * Author: ktakeda
 *
 * Created on 2010/04/13, 19:41
 */
#include "TypeDef.h"
#include "BoundaryEdge.h"

#include "BoundaryMesh.h"

// Equivalent Nodal Force(等価節点力)
#include "ShapeLine.h"

#include "Node.h"
#include "Element.h"

namespace pmw{
#ifndef _BOUNDARYEDGEMESH_H
#define	_BOUNDARYEDGEMESH_H
class CBoundaryEdgeMesh:public CBoundaryMesh{
public:
    CBoundaryEdgeMesh();
    virtual ~CBoundaryEdgeMesh();

protected:
    vector<CBoundaryEdge*> mvBEdge;
    map<uint, uint, less<uint> > mmBEdgeID2Index;

    vvuint mvAggregateEdge;// BNodeの辺-集合

    vector<CBoundaryNode*> mvBEdgeBNode;// 辺ノード

public:
    // 境界辺 :=> BEdge
    void resizeEdge(const uint& res_size);
    uint getNumOfEdge(){ return mvBEdge.size();}

    void setBEdge(const uint& index, CBoundaryEdge *pBEdge);
    void addBEdge(CBoundaryEdge *pBEdge);

    CBoundaryEdge* getBEdgeIX(const uint& index){ return mvBEdge[index];}
    CBoundaryEdge* getBEdgeID(const uint& id){
        uint index= mmBEdgeID2Index[id];
        return mvBEdge[index];
    }
    uint& getBEdgeIndex(const uint& id){ return mmBEdgeID2Index[id];}

    // 辺-集合
    void resizeAggEdge();
    vuint& getAggEdge(const uint& ibnode){ return mvAggregateEdge[ibnode];}
    void setupAggEdge();
    void setAggEdge(const uint& ibnode, const uint& nEdgeID){ mvAggregateEdge[ibnode].push_back(nEdgeID);}

    // 辺ノード生成
    void GeneEdgeBNode();
    // 辺ノード提供
    CBoundaryNode* getEdgeBNode(const uint& iedge){ return mvBEdgeBNode[iedge];}


    // Refine 辺の再分割  => progBMeshにセット
    void refine(CBoundaryEdgeMesh *pProgEdgeMesh);

    // Refine時のmvBEdgeBNodeの解放
    void deleteProgData();

protected:
    //Neumann条件の節点分配(形状関数による等価節点力:EquivalentNodeForce)
    //----
    virtual void distNeumannValue();

    //Dirichlet条件の節点分配(BNode周囲の面集合平均値)
    //----
    virtual void distDirichletValue_at_CGrid();//Level==0(coarse grid)
    virtual void distDirichletValue_at_FGrid();//Level>=1(fine grid)
};
#endif	/* _BOUNDARYEDGEMESH_H */
}


