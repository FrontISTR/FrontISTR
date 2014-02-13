/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryFaceMesh.h
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
#include <map>
#include <utility>
#include "BoundaryFace.h"
#include "ShapeQuad.h"
#include "ShapeTriangle.h"
#include "BoundaryMesh.h"
typedef std::pair<pmw::CBoundaryNode*, pmw::CBoundaryNode*> PairBNode;

namespace pmw
{

#ifndef _BOUNDARYFACEMESH_H
#define	_BOUNDARYFACEMESH_H

class CBoundaryFaceMesh:public CBoundaryMesh
{
public:
    CBoundaryFaceMesh();
    virtual ~CBoundaryFaceMesh();
protected:
    vector<CBoundaryFace*> mvBFace;
    map<uiint, uiint, less<uiint> > mmBFaceID2Index;
    vvuint mvAggregateFace;
    vector<CBoundaryNode*> mvBEdgeBNode;
    vector<CBoundaryNode*> mvBFaceBNode;
public:
    void resizeBFace(const uiint& res_size) {
        mvBFace.resize(res_size);
    }
    void setBFace(const uiint& index, CBoundaryFace *pBFace);
    void addBFace(CBoundaryFace *pBFace);
    uiint getNumOfBFace() {
        return mvBFace.size();
    }
    CBoundaryFace* getBFaceIX(const uiint& index) {
        return mvBFace[index];
    }
    CBoundaryFace* getBFaceID(const uiint& id) {
        uiint index= mmBFaceID2Index[id];
        return mvBFace[index];
    }
    uiint& getBFaceIndex(const uiint& id) {
        return mmBFaceID2Index[id];
    }
public:
    void resizeAggFace();
    vuint& getAggFace(const uiint& ibnode) {
        return mvAggregateFace[ibnode];
    }
    void setupAggFace();
    void setAggFace(const uiint& ibnode, const uiint& nFaceID) {
        mvAggregateFace[ibnode].push_back(nFaceID);
    }
    void GeneEdgeBNode();
    void GeneFaceBNode();
    void refine(CBoundaryFaceMesh *pProgBFaceMesh);
    void deleteProgData();
protected:
    virtual void distNeumannValue();
    virtual void distDirichletValue();
};
#endif	/* _BOUNDARYFACEMESH_H */
}
