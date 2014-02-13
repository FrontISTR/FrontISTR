/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/CommunicationMesh.h
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
#ifndef _COMMMesh_H_39e5ebb0_da94_4efc_8616_8c2c5918cc4f
#define _COMMMesh_H_39e5ebb0_da94_4efc_8616_8c2c5918cc4f
#include "TypeDef.h"
#include "Element.h"
#include "Node.h"
#include "CommElement.h"
#include <map>
#include "IndexBucket.h"
#include "Logger.h"
#include "QuickSort.h"
namespace pmw
{
class CCommMesh
{
public:
    CCommMesh(CIndexBucket* pBucket);
    virtual ~CCommMesh();
protected:
    uiint mCommID;
    uiint mRankID;
    uiint mTransmitRankID;
    vector<CCommElement*> mvCommElementAll;
    vector<CCommElement*> mvCommElement;
    vector<CCommElement*> mvDCommElement;
    vector<CNode*>  mvNode;
    CIndexBucket *mpBucket;
    vuint mvNodeRank;
    vector<CNode*>    mvSendNode;
    vector<CNode*>    mvRecvNode;
    vuint mvSendCommNodeID;
    vuint mvRecvCommNodeID;
    map<uiint, uiint, less<uiint> > mmCommNodeIX;
    map<uiint, uiint, less<uiint> > mmCommElementIX;
    vector<CNode*>    mvDNode;
    vector<CElement*> mvDElement;
    Utility::CLogger* mpLogger;
public:
    void setCommID(const uiint& comID) {
        mCommID= comID;
    }
    uiint& getCommID() {
        return mCommID;
    }
    void setRankID(const uiint& rank) {
        mRankID= rank;
    }
    uiint& getRankID() {
        return mRankID;
    }
    void setTransmitRankID(const uiint& rank) {
        mTransmitRankID = rank;
    }
    uiint& getTransmitRankID() {
        return mTransmitRankID;
    }
    void reserveCommElementAll(const uiint& res_size) {
        mvCommElementAll.reserve(res_size);
    }
    void setCommElementAll(CCommElement* pCommElem) {
        mvCommElementAll.push_back(pCommElem);
    }
    uiint getNumOfCommElementAll() {
        return mvCommElementAll.size();
    }
    CCommElement* getCommElementAll(const uiint& all_com_index) {
        return mvCommElementAll[all_com_index];
    }
    CElement* getElementAll(const uiint& all_com_index) {
        return mvCommElementAll[all_com_index]->getElement();
    }
    uiint getNumOfCommElement() {
        return mvCommElement.size();
    }
    CCommElement* getCommElement(const uiint& com_index) {
        return mvCommElement[com_index];
    }
    CElement* getElement(const uiint& com_index) {
        return mvCommElement[com_index]->getElement();
    }
    uiint getNumOfDNode() {
        return mvDNode.size();
    }
    CNode* getDNode(const uiint& dcom_index) {
        return mvDNode[dcom_index];
    }
    uiint getNumOfDCommElement() {
        return mvDCommElement.size();
    }
    CCommElement* getDCommElement(const uiint& dcom_index) {
        return mvDCommElement[dcom_index];
    }
    CElement* getDElement(const uiint& dcom_index) {
        return mvDElement[dcom_index];
    }
    uiint getNumOfNode() {
        return mvNode.size();
    }
    void reserveNode(const uiint& res_size) {
        mvNode.reserve(res_size);
    }
    void setNode(CNode* pNode) {
        mvNode.push_back(pNode);
    }
    CNode* getNode(const uiint& com_index) {
        return mvNode[com_index];
    }
    void setSendNode(CNode* pNode, const uiint& comNodeID) {
        mvSendCommNodeID.push_back(comNodeID);
        mvSendNode.push_back(pNode);
    }
    CNode* getSendNodeIX(const uiint& index) {
        return mvSendNode[index];
    }
    uiint getNumOfSendNode() {
        return mvSendNode.size();
    }
    uiint& getSendCommNodeID(const uiint& index) {
        return mvSendCommNodeID[index];
    }
    void setRecvNode(CNode* pNode, const uiint& comNodeID) {
        mvRecvCommNodeID.push_back(comNodeID);
        mvRecvNode.push_back(pNode);
    }
    CNode* getRecvNodeIX(const uiint& index) {
        return mvRecvNode[index];
    }
    uiint getNumOfRecvNode() {
        return mvRecvNode.size();
    }
    uiint& getRecvCommNodeID(const uiint& index) {
        return mvRecvCommNodeID[index];
    }
    void resizeNodeRank(const uiint& res_size) {
        mvNodeRank.resize(res_size);
    }
    void setNodeRank(const uiint& commNodeID, const uiint& rank) {
        mvNodeRank[commNodeID]= rank;
    }
    uiint& getNodeRank(const uiint& commNodeID) {
        return  mvNodeRank[commNodeID];
    }
    void AllocateCommElement();
    void setupAggCommElement(vector<CElement*> &vElement);
    void sortCommNodeIndex();
    void setupMapID2CommID();
};
}
#endif	/* _COMMMesh_H */
