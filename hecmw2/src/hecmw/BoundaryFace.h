/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryFace.h
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
#include "BoundaryParts.h"
#include <utility>
#include "EdgeTree.h"
#include "ElementProperty.h"
typedef std::pair<pmw::CBoundaryNode*, pmw::CBoundaryNode*> PairBNode;

namespace pmw
{
#ifndef _BOUNDARYFACE_H_c5c754e6_
#define	_BOUNDARYFACE_H_c5c754e6_

class CBoundaryFace:public CBoundaryParts
{
public:
    CBoundaryFace();
    virtual ~CBoundaryFace();
protected:
    uiint mnElemFaceID;
    vuint mvEdgeNeibFace;
    vector<CBoundaryNode*> mvEdgeBNode;
    CBoundaryNode *mpFaceBNode;
public:
    void setElementFaceID(const uiint& id) {
        mnElemFaceID= id;
    }
    uiint& getElementFaceID() {
        return mnElemFaceID;
    }
    void setBFaceShape(const uiint& elemType);
    uiint& getBFaceShape() {
        return mnShapeType;
    }
    virtual uiint getNumOfVert();
protected:
    bool* mvbMarkingEdge;
public:
    void markingEdge(const uiint& iedge);
    bool isMarkingEdge(const uiint& iedge) {
        return mvbMarkingEdge[iedge];
    }
    void setEdgeNeibFace(const uiint& iedge, const uiint& neibFaceID);
    void setEdgeBNode(const uiint& iedge, CBoundaryNode *pEdgeBNode);
    uiint getNumOfEdge();
    PairBNode getPairBNode(const uiint& iedge);
    uiint& getEdgeID(PairBNode& pairBNode);
    void setFaceBNode(CBoundaryNode *pFaceBNode);
    void setupNode_Edge();
    void setupNode_Face();
protected:
    vector<CBoundaryFace*> mvProgBFace;
    double mArea;
    double  triArea(CNode* pNode0, CNode* pNode1, CNode* pNode2);
public:
    void refine(uiint& countID, const vuint& vDOF);
    double& calcArea();
    double& getArea() {
        return mArea;
    }
    vector<CBoundaryFace*>& getProgParts() {
        return mvProgBFace;
    }

    void distDirichletVal(const uiint& dof, const uiint& mgLevel, const uiint& nMaxMGLevel, vuint& vPolandDOF, map<uiint,CPoland*>& mPoland);

    void replaceEdgeBNode();
    void deleteProgData();

};
#endif	/* _BOUNDARYFACE_H */
}
