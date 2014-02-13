/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/SkinFace.h
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
#include "ContactNode.h"
#include "ElementType.h"
#include "ElementProperty.h"
#include "Node.h"
#include "Element.h"
#include "EdgeTree.h"
#include <utility>
typedef std::pair<pmw::CContactNode*, pmw::CContactNode*> PairConNode;
#include "Logger.h"
namespace pmw
{
#ifndef _SKINFACE_H
#define	_SKINFACE_H
class CSkinFace
{
public:
    CSkinFace();
    virtual ~CSkinFace();
protected:
    uiint mID;
    uiint mRank;
    uiint mLevel;
    bool mbSelfDom;
    uiint mMeshID;
    uiint mElementID;
    uiint mElementFaceID;
    vector<CContactNode*> mvConNode;
    vector<CContactNode*> mvEdgeNode;
    bool* mvbEdgeMarking;
    vector<CSkinFace*>    mvEdgeFace;
    uiint& getEdgeIndex(PairConNode& pairConNode);
    CContactNode *mpFaceNode;
    vector<CSkinFace*> mvProgFace;
    uiint mNumOfEdge;
    uiint mShapeType;
    uiint mnOrder;
    virtual CSkinFace* generateFace();
    CSkinFace* mpOtherFace;
    vdouble mvNormalVector;
public:
    virtual const char* getName() {
        return "SkinFace";
    }

    void setLevel(const uiint& level) {
        mLevel= level;
    }
    uiint& getLevel() {
        return mLevel;
    }

    void setRank(const uiint& rank) {
        mRank= rank;
    }
    uiint& getRank() {
        return mRank;
    }

    void setID(const uiint& id) {
        mID= id;
    }
    uiint& getID() {
        return mID;
    }

    void markingSelf() {
        mbSelfDom=true;
    }
    bool isSelf() {
        return mbSelfDom;
    }

    void setMeshID(const uiint& meshID) {
        mMeshID= meshID;
    }
    uiint& getMeshID() {
        return mMeshID;
    }

    void  setElementID(const uiint& id) {
        mElementID = id;
    }
    uiint& getElementID() {
        return mElementID;
    }

    void  setFaceID(const uiint& faceID) {
        mElementFaceID= faceID;
    }
    uiint& getFaceID() {
        return mElementFaceID;
    }

    void resizeNode(const uiint& size) {
        mvConNode.resize(size);
    }
    void addNode(CContactNode* pConNode) {
        mvConNode.push_back(pConNode);
    }
    void setNode(CContactNode* pConNode, const uiint& index) {
        mvConNode[index] = pConNode;
    }
    void setNode(vector<CContactNode*>& vConNode) {
        mvConNode = vConNode;
    }
    void operator=(vector<CContactNode*>& vConNode) {
        mvConNode = vConNode;
    }
    uiint   getNumOfNode() {
        return mvConNode.size();
    }
    vector<CContactNode*>&  getNodes() {
        return mvConNode;
    }
    CContactNode*  getNode(const uiint& index) {
        return mvConNode[index];
    }

    void setShapeType(const uiint& shapeType);
    uiint& getShapeType() {
        return mShapeType;
    }
    uiint& getNumOfEdge() {
        return mNumOfEdge;
    }
    uiint  getNumOfVert();
    uiint& getOrder() {
        return mnOrder;
    }

    virtual void addSlaveNode(CContactNode* pConNode);
    virtual CContactNode* getSlaveNode(const uiint& index) {
        return NULL;
    }
    virtual uiint getNumOfSlaveNode() {
        return 0;
    }

    virtual void CalcSlave(const uiint& islave, const uiint& valType);

    virtual double& getCoef(const uiint& slaveID, const uiint& ivert);

    PairConNode getEdgePairNode(const uiint& iedge);
    void setEdgeFace(CSkinFace* pFace, PairConNode& pairConNode);
    void setEdgeFace(CSkinFace* pFace, const uiint& iedge);
    void setEdgeConNode(CContactNode* pEdgeConNode, PairConNode& pairConNode);
    void setEdgeConNode(CContactNode* pEdgeConNode, const uiint& iedge);
    bool isEdgeNodeMarking(const uiint& iedge) {
        return mvbEdgeMarking[iedge];
    }
    bool isEdgeNodeMarking(PairConNode& pairConNode);
    void markingEdgeNode(const uiint& iedge) {
        mvbEdgeMarking[iedge]=true;
    }
    void markingEdgeNode(PairConNode& pairConNode);
    void setFaceConNode(CContactNode* pFaceConNode) {
        mpFaceNode= pFaceConNode;
    }
    CContactNode* getEdgeConNode(const uiint& iedge) {
        return mvEdgeNode[iedge];
    }
    CContactNode* getEdgeConNode(PairConNode& pairConNode);
    CContactNode* getFaceConNode() {
        return mpFaceNode;
    }
    vector<CSkinFace*>& getProgFace() {
        return mvProgFace;
    }
    CSkinFace* getProgFace(const uiint& ivert) {
        return mvProgFace[ivert];
    }

    void refine(CElement* pElem, uiint& faceID);

protected:
    void setupNodeID_progFace(CElement* pElem, const uiint& numOfVert);
    void setupEdgeNodeID(CElement* pElem, const uiint& numOfVert);
public:
    void setupNodeID_2nd_LastLevel(CElement* pElem);

public:
    vdouble& CalcNzNormalVector();
    vdouble& getNzNormalVector();

    void replaceEdgeNode();
    void deleteProgData();
};
#endif	/* _SKINFACE_H */
}
