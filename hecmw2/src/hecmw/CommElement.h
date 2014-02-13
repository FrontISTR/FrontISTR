/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/CommElement.h
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
#include "Node.h"
#include "Element.h"
#include "NodeConnectNodeTree.h"
#include "EdgeTree.h"
#include "FaceTree.h"
#include "EdgeFaceTree.h"
#include "ProgElementTree.h"
#include "Logger.h"
#include "ElementType.h"
namespace pmw
{
#ifndef _COMMELEMENT_H
#define	_COMMELEMENT_H
class CCommElement
{
public:
    CCommElement();
    virtual ~CCommElement();
protected:
    uiint mID;
    CElement* mpElement;
    vector<vector<CCommElement*> > mvvAggCommElem;
    vvuint mvvNeibCommElemVert;
    bool* mvbSend;
    bool* mvbRecv;
    bool* mvbOther;
    vector<CNode*> mvSendNode;
    vector<CNode*> mvRecvNode;
    vector<CNode*> mvOtherNode;
    vuint mvCommNodeIndex;
    bool mbCommunication;
    bool mbRCommunication;
    bool* mvbNodeIXCheck;
    bool* mvbDNodeMarking;
    Utility::CLogger *mpLogger;
    vuint mvNodeRank;
    vuint mvEdgeRank;
    vuint mvFaceRank;
    uiint  mVolRank;
public:
    void setupProgNodeRank(const uiint& mgLevel);
    void setID(const uiint& id) {
        mID= id;
    }
    uiint& getID() {
        return mID;
    }
    virtual uiint getShapeType()=0;
    virtual uiint getBaseShapeType()=0;
    const uiint& getNumOfVert() {
        return mpElement->getNumOfNode();
    }
    const uiint& getNumOfEdge() {
        return mpElement->getNumOfEdge();
    }
    const uiint& getNumOfFace() {
        return mpElement->getNumOfFace();
    }
    void setNodeRank(const uiint& ivert, const uiint& rank) {
        mvNodeRank[ivert]=rank;
    }
    uiint& getNodeRank(const uiint& ivert) {
        return mvNodeRank[ivert];
    }
    void setElement(CElement* pElem) {
        mpElement= pElem;
    }
    CElement* getElement() {
        return mpElement;
    }
    void setCommNodeIndex(const uiint& ivert, uiint& comID, vector<CNode*>& vCommMeshNode, vuint& vCommMeshNodeRank);
protected:
    void setNeibCommNodeIndex(const uiint& ivert, const uiint& comID);
public:
    uiint& getCommNodeIndex(const uiint& ivert) {
        return mvCommNodeIndex[ivert];
    }
    void getDNode(const uiint& ivert, vector<CNode*>& vDNode);
protected:
    void markingDNode(const uiint& ivert);
public:
    bool isMarkingDNode(const uiint& ivert) {
        return mvbDNodeMarking[ivert];
    }
    CNode* getNode(const uiint& ivert) {
        return mpElement->getNode(ivert);
    }
    vector<CNode*>& getSendNode() {
        return mvSendNode;
    }
    vector<CNode*>& getRecvNode() {
        return mvRecvNode;
    }
    uiint getNumOfSendNode() {
        return mvSendNode.size();
    }
    uiint getNumOfRecvNode() {
        return mvRecvNode.size();
    }
    CNode* getSendNode(const uiint& i) {
        return mvSendNode[i];
    }
    CNode* getRecvNode(const uiint& i) {
        return mvRecvNode[i];
    }
    void setAggCommElement(const uiint& ivert, CCommElement* pComElem) {
        mvvAggCommElem[ivert].push_back(pComElem);
    }
    void setNeibCommElemVert(const uiint& ivert, const uiint& neibVert) {
        mvvNeibCommElemVert[ivert].push_back(neibVert);
    }
    bool isCommElement() {
        return mbCommunication;
    }
    bool isRCommElement() {
        return mbRCommunication;
    }
    void sortNodeRank(const uiint& myRank, const uiint& transRank);
    virtual bool isTypeCoincidence()=0;
    vuint& getEdgeRank() {
        return mvEdgeRank;
    }
    uiint& getEdgeRank(const uiint& iedge) {
        return mvEdgeRank[iedge];
    }
    vuint& getFaceRank() {
        return mvFaceRank;
    }
    uiint& getFaceRank(const uiint& iface) {
        return mvFaceRank[iface];
    }
    uiint& getVolRank() {
        return mVolRank;
    }
};
#endif	/* _COMMELEMENT_H */
}
