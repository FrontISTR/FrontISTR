/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/CommFace.h
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
#include <utility>
#include "EdgeTree.h"
#include "CommNode.h"
#include "Element.h"
#include "Logger.h"
namespace pmw
{
typedef pair<CCommNode*,CCommNode*> PairCommNode;
#ifndef _COMMFACE_H
#define	_COMMFACE_H
class CCommFace
{
public:
    CCommFace();
    virtual ~CCommFace();
private:
    uiint mID;
    uiint mMGLevel;
    uiint mMeshID;
    uiint mElementID;
////    uiint mElementFaceID;
    uiint mNumOfEdge;
    uiint mFaceType;

    vector<CCommNode*>  mvCommNode;
    vector<CCommNode*>  mvEdgeCommNode;
    CCommNode*          mpFaceCommNode;
    vector<CCommFace*>  mvEdgeCommFace;

    bool*  mvbEdgeMarking;

    uiint& getEdgeIndex(PairCommNode& pairCommNode);
    vector<CCommFace*> mvProgCommFace;

public:
    void initialize(const uiint& numOfVert, const uiint& numOfEdge, const uiint& nOrder);
    void  setID(const uiint& id) {
        mID= id;
    }
    uiint& getID() {
        return mID;
    }

    void setMGLevel(const uiint& mgLevel) {
        mMGLevel= mgLevel;
    }
    uiint& getMGLevel() {
        return mMGLevel;
    }

    uiint& getType() {
        return mFaceType;
    }
    uiint& getNumOfEdge() {
        return mNumOfEdge;
    }
    uiint  getNumOfVert();
    uiint getOrder();

    void setMeshID(const uiint& meshID) {
        mMeshID= meshID;
    }
    void setElementID(const uiint& elemID) {
        mElementID= elemID;
    }
////////    void setElementFaceID(const uiint& iface){ mElementFaceID= iface;}

    uiint& getMeshID() {
        return mMeshID;
    }
    uiint& getElementID() {
        return mElementID;
    }
    uiint& getElementFaceID(CElement* pElem);//Factoryの面中心ノード取得で使用

    void setCommNode(const uiint& index, CCommNode* pCommNode) {
        mvCommNode[index]= pCommNode;
    }
    uiint getCommNodeSize() {
        return mvCommNode.size();
    }
    CCommNode* getCommNode(const uiint& index) {
        return mvCommNode[index];
    }
    vector<CCommNode*>& getCommNode() {
        return mvCommNode;
    }

    PairCommNode getEdgePairCommNode(const uiint& iedge);
    CCommNode* getEdgeCommNode(const uiint& iedge) {
        return  mvEdgeCommNode[iedge];
    }
    void setEdgeCommFace(CCommFace* pNeibFace, const uiint& iedge);
    void setEdgeCommFace(CCommFace* pNeibFace, PairCommNode& pairCommNode);
    void setEdgeCommNode(CCommNode* pEdgeCommNode, const uiint& iedge);
    void setEdgeCommNode(CCommNode* pEdgeCommNode, PairCommNode& pairCommNode);
    void markingEdgeNode(const uiint& iedge);
    void markingEdgeNode(PairCommNode& pairCommNode);
    bool isEdgeNodeMarking(const uiint& iedge) {
        return mvbEdgeMarking[iedge];
    }
    void replaceEdgeCommNode();

    void setFaceCommNode(CCommNode* pFaceCommNode) {
        mpFaceCommNode= pFaceCommNode;
    }
    CCommNode* getFaceCommNode() {
        return mpFaceCommNode;
    }
    vector<CCommFace*>& refine(CElement* pElement);

    void deleteProgData();
};
#endif	/* _COMMFACE_H */
}
