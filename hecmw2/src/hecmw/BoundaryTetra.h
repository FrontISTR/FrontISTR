/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryTetra.h
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
#include "BoundaryVolume.h"
#include "BoundaryHexa.h"
namespace pmw
{
#ifndef _BOUNDARYTETRA_H
#define	_BOUNDARYTETRA_H
class CBoundaryTetra:public CBoundaryVolume
{
public:
    CBoundaryTetra();
    virtual ~CBoundaryTetra();
protected:
    virtual uiint* getLocalNode_Edge(const uiint& iedge);
    virtual uiint* getLocalNode_Face(const uiint& iface);
public:
    virtual uiint getElemType();
    virtual uiint getNumOfEdge();
    virtual uiint getNumOfFace();
    virtual uiint getNumOfNode();
    virtual uiint getNumOfVert();
    virtual void setOrder(const uiint& order);
    virtual PairBNode getPairBNode(const uiint& iedge);
    virtual uiint& getEdgeID(PairBNode& pairBNode);
    virtual vector<CBoundaryNode*> getFaceCnvNodes(const uiint& iface);
    virtual uiint& getFaceID(vector<CBoundaryNode*>& vBNode);
    virtual void refine(uiint& countID, const vuint& vDOF);
    virtual double& calcVolume();
    virtual void distDirichletVal(const uiint& dof, const uiint& mgLevel, const uiint& nMaxMGLevel, vuint& vPolandDOF, map<uiint,CPoland*>& mPoland);
    virtual void replaceEdgeBNode(const uiint& iedge);
    virtual void deleteProgData();
};
#endif	/* _BOUNDARYTETRA_H */
}
