/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   CommElement.h
|
|                     Written by T.Takeda,    2010/06/01
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
namespace pmw{
#ifndef _COMMELEMENT_H
#define	_COMMELEMENT_H
class CCommElement{
public:
    CCommElement();
    virtual ~CCommElement();
protected:
    uint mID;
    CElement* mpElement;
    vector<vector<CCommElement*> > mvvAggCommElem;
    vvuint mvvNeibCommElemVert;                   
    vbool mvbSend;
    vbool mvbRecv;
    vbool mvbOther;
    vector<CNode*> mvSendNode;
    vector<CNode*> mvRecvNode;
    vector<CNode*> mvOtherNode;
    vuint mvCommNodeIndex;
    bool mbCommunication;
    bool mbRCommunication;
    vector<bool> mvbNodeIXCheck;
    vector<bool> mvbDNodeMarking;
    Utility::CLogger *mpLogger;
    vuint mvNodeRank;
    vuint mvEdgeRank;
    vuint mvFaceRank;
     uint  mVolRank; 
public:
    void setupProgNodeRank(const uint& mgLevel);
    void setID(const uint& id){ mID= id;}
    uint& getID(){ return mID;}
    virtual uint getShapeType()=0;
    virtual uint getBaseShapeType()=0;
    const uint& getNumOfVert(){ return mpElement->getNumOfNode();}
    const uint& getNumOfEdge(){ return mpElement->getNumOfEdge();}
    const uint& getNumOfFace(){ return mpElement->getNumOfFace();}
    void setNodeRank(const uint& ivert, const uint& rank){ mvNodeRank[ivert]=rank;}
    uint& getNodeRank(const uint& ivert){ return mvNodeRank[ivert];}
    void setElement(CElement* pElem){ mpElement= pElem;}
    CElement* getElement(){ return mpElement;}
    void setCommNodeIndex(const uint& ivert, uint& comID, vector<CNode*>& vCommMeshNode, vuint& vCommMeshNodeRank);
protected:
    void setNeibCommNodeIndex(const uint& ivert, const uint& comID);
public:
    uint& getCommNodeIndex(const uint& ivert){ return mvCommNodeIndex[ivert];}
    void getDNode(const uint& ivert, vector<CNode*>& vDNode);
protected:
    void markingDNode(const uint& ivert);
public:
    bool isMarkingDNode(const uint& ivert){ return mvbDNodeMarking[ivert];}
    CNode* getNode(const uint& ivert){ return mpElement->getNode(ivert);}
    vector<CNode*>& getSendNode(){ return mvSendNode;}
    vector<CNode*>& getRecvNode(){ return mvRecvNode;}
    uint getNumOfSendNode(){ return mvSendNode.size();}
    uint getNumOfRecvNode(){ return mvRecvNode.size();}
    CNode* getSendNode(const uint& i){ return mvSendNode[i];}
    CNode* getRecvNode(const uint& i){ return mvRecvNode[i];}
    void setAggCommElement(const uint& ivert, CCommElement* pComElem){ mvvAggCommElem[ivert].push_back(pComElem);}
    void setNeibCommElemVert(const uint& ivert, const uint& neibVert){ mvvNeibCommElemVert[ivert].push_back(neibVert);}
    bool isCommElement(){ return mbCommunication;}
    bool isRCommElement(){ return mbRCommunication;}
    void sortNodeRank(const uint& myRank, const uint& transRank);
    virtual bool isTypeCoincidence()=0;
    vuint& getEdgeRank(){ return mvEdgeRank;}
     uint& getEdgeRank(const uint& iedge){ return mvEdgeRank[iedge];}
    vuint& getFaceRank(){ return mvFaceRank;}
     uint& getFaceRank(const uint& iface){ return mvFaceRank[iface];}
     uint& getVolRank(){ return mVolRank;}
};
#endif	/* _COMMELEMENT_H */
}
