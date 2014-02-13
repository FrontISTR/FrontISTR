/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryEdge.h
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
typedef std::pair<pmw::CBoundaryNode*, pmw::CBoundaryNode*> PairBNode;
namespace pmw
{
#ifndef _BOUNDARYEDGE_H
#define	_BOUNDARYEDGE_H
class CBoundaryEdge:public CBoundaryParts
{
public:
    CBoundaryEdge();
    virtual ~CBoundaryEdge();
protected:
    uiint mnElemEdgeID;
    CBoundaryNode *mpEdgeBNode;
    vector<CBoundaryEdge*> mvProgBEdge;
    double mLength;
public:
    void setElementEdgeID(const uiint& id) {
        mnElemEdgeID= id;
    }
    uiint& getElementEdgeID() {
        return mnElemEdgeID;
    }
    void setBEdgeShape(const uiint& elemType);
    uiint& getBEdgeShape() {
        return mnShapeType;
    }
    virtual uiint getNumOfVert();
    void setEdgeBNode(CBoundaryNode *pBNode);
    CBoundaryNode* getEdgeBNode() {
        return mpEdgeBNode;
    }
    PairBNode getPairBNode();
    void setupNode();
    double& calcLength();
    double& getLength() {
        return mLength;
    }
    void refine(uiint& countID, const vuint& vDOF);
    vector<CBoundaryEdge*>& getProgParts() {
        return mvProgBEdge;
    }
    void replaceEdgeBNode();

    void distDirichletVal(const uiint& dof, const uiint& mgLevel, const uiint& nMaxMGLevel, vuint& vPolandDOF, map<uiint,CPoland*>& mPoland);

};
#endif	/* _BOUNDARYEDGE_H */
}
