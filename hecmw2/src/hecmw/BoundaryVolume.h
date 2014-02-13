/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryVolume.h
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
#include "EdgeTree.h"
#include "FaceTree.h"
#include "DiscreteVolume.h"
#include "BoundaryParts.h"

#include "Poland.h"

typedef std::pair<pmw::CBoundaryNode*, pmw::CBoundaryNode*> PairBNode;
namespace pmw
{
#ifndef _BOUNDARYELEMENT_H_
#define	_BOUNDARYELEMENT_H_
class CBoundaryVolume:public CBoundaryParts
{
public:
    CBoundaryVolume();
    virtual ~CBoundaryVolume();
protected:
    bool* mvbMarkingEdge;
    bool* mvbMarkingFace;
    vector<CBoundaryNode*> mvEdgeBNode;
    vector<CBoundaryNode*> mvFaceBNode;
    CBoundaryNode          *mpVolBNode;
    vvuint mvEdgeNeibVol;
    vuint  mvFaceNeibVol;
public:
    void markingEdge(const uiint& iedge) {
        mvbMarkingEdge[iedge]=true;
    }
    void markingFace(const uiint& iface) {
        mvbMarkingFace[iface]=true;
    }
    bool isMarkingEdge(const uiint& iedge) {
        return mvbMarkingEdge[iedge];
    }
    bool isMarkingFace(const uiint& iface) {
        return mvbMarkingFace[iface];
    }
    virtual uiint getElemType()=0;
    virtual uiint getNumOfEdge()=0;
    virtual uiint getNumOfFace()=0;
    virtual uiint getNumOfNode()=0;
    virtual void setOrder(const uiint& order)=0;
    virtual PairBNode getPairBNode(const uiint& iedge)=0;
    virtual uiint& getEdgeID(PairBNode& pairBNode)=0;
    virtual vector<CBoundaryNode*> getFaceCnvNodes(const uiint& iface)=0;
    virtual uiint& getFaceID(vector<CBoundaryNode*>& vBNode)=0;
    void setEdgeNeibVol(const uiint& iedge, const uiint& neibVolID);
    void setFaceNeibVol(const uiint& iface, const uiint& neibVolID);
    vuint& getEdgeNeibVolID(const uiint& iedge) {
        return mvEdgeNeibVol[iedge];
    }
    uiint&  getFaceNeibVolID(const uiint& iface) {
        return mvFaceNeibVol[iface];
    }
    void setEdgeBNode(const uiint& iedge, CBoundaryNode *pBNode);
    void setFaceBNode(const uiint& iface, CBoundaryNode *pBNode);
    void setVolBNode(CBoundaryNode *pBNode);
    CBoundaryNode* getEdgeBNode(const uiint& iedge) {
        return mvEdgeBNode[iedge];
    }
    CBoundaryNode* getFaceBNode(const uiint& iface) {
        return mvFaceBNode[iface];
    }
    CBoundaryNode* getVolBNode() {
        return mpVolBNode;
    }
protected:
    vector<CBoundaryVolume*> mvProgVolume;
    double mCubicVolume;
    double tetraVolume(CNode* pNode0, CNode* pNode1, CNode* pNode2, CNode* pNode3);
    uiint dividHexa(const uiint& iprog, CBoundaryVolume* pProgVol);
    uiint dividTetra(const uiint& iprog, CBoundaryVolume* pProgVol);
    uiint dividPrism(const uiint& iprog, CBoundaryVolume* pProgVol);
    void distValue(CBoundaryVolume *pProgVol, const double& coef, const vuint& vDOF);
    virtual uiint* getLocalNode_Edge(const uiint& iedge)=0;
    virtual uiint* getLocalNode_Face(const uiint& iface)=0;
public:
    virtual void refine(uiint& countID, const vuint& vDOF)=0;
    vector<CBoundaryVolume*>& getProgParts() {
        return mvProgVolume;
    }
    virtual double& calcVolume()=0;
    double& getCubicVolume() {
        return mCubicVolume;
    };
    virtual void distDirichletVal(const uiint& dof, const uiint& mgLevel, const uiint& nMaxMGLevel, vuint& vPolandDOF, map<uiint,CPoland*>& mPoland)=0;
    virtual void replaceEdgeBNode(const uiint& iedge)=0;
    virtual void deleteProgData()=0;
};
#endif	/* _BOUNDARYELEMENT_H */
}
