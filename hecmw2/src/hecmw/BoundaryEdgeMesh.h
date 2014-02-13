/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryEdgeMesh.h
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
#include "TypeDef.h"
#include "BoundaryEdge.h"
#include "BoundaryMesh.h"
#include "ShapeLine.h"
#include "Node.h"
#include "Element.h"
namespace pmw
{
#ifndef _BOUNDARYEDGEMESH_H
#define	_BOUNDARYEDGEMESH_H
class CBoundaryEdgeMesh:public CBoundaryMesh
{
public:
    CBoundaryEdgeMesh();
    virtual ~CBoundaryEdgeMesh();
protected:
    vector<CBoundaryEdge*> mvBEdge;
    map<uiint, uiint, less<uiint> > mmBEdgeID2Index;
    vvuint mvAggregateEdge;
    vector<CBoundaryNode*> mvBEdgeBNode;
public:
    void resizeEdge(const uiint& res_size);
    uiint getNumOfEdge() {
        return mvBEdge.size();
    }
    void setBEdge(const uiint& index, CBoundaryEdge *pBEdge);
    void addBEdge(CBoundaryEdge *pBEdge);
    CBoundaryEdge* getBEdgeIX(const uiint& index) {
        return mvBEdge[index];
    }
    CBoundaryEdge* getBEdgeID(const uiint& id) {
        uiint index= mmBEdgeID2Index[id];
        return mvBEdge[index];
    }
    uiint& getBEdgeIndex(const uiint& id) {
        return mmBEdgeID2Index[id];
    }
    void resizeAggEdge();
    vuint& getAggEdge(const uiint& ibnode) {
        return mvAggregateEdge[ibnode];
    }
    void setupAggEdge();
    void setAggEdge(const uiint& ibnode, const uiint& nEdgeID) {
        mvAggregateEdge[ibnode].push_back(nEdgeID);
    }
    void GeneEdgeBNode();
    CBoundaryNode* getEdgeBNode(const uiint& iedge) {
        return mvBEdgeBNode[iedge];
    }
    void refine(CBoundaryEdgeMesh *pProgEdgeMesh);
    void deleteProgData();
protected:
    virtual void distNeumannValue();
    virtual void distDirichletValue();
};
#endif	/* _BOUNDARYEDGEMESH_H */
}
