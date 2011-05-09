/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   CommunicationMesh.h
|
|                     Written by T.Takeda,    2010/06/01
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
namespace pmw{
class CCommMesh{
public:
    CCommMesh(CIndexBucket* pBucket);
    virtual ~CCommMesh();
protected:
    uint mCommID;
    uint mRankID;
    uint mTransmitRankID;
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
    map<uint, uint, less<uint> > mmCommNodeIX;
    map<uint, uint, less<uint> > mmCommElementIX;
    vector<CNode*>    mvDNode;   
    vector<CElement*> mvDElement;
    Utility::CLogger* mpLogger;
public:
    void setCommID(const uint& comID){ mCommID= comID;}
    uint& getCommID(){ return mCommID;}
    void setRankID(const uint& rank){ mRankID= rank;}
    uint& getRankID(){ return mRankID;}
    void setTransmitRankID(const uint& rank){ mTransmitRankID = rank;}
    uint& getTransmitRankID(){ return mTransmitRankID;}
    void reserveCommElementAll(const uint& res_size){ mvCommElementAll.reserve(res_size);}
    void setCommElementAll(CCommElement* pCommElem){ mvCommElementAll.push_back(pCommElem);}
    uint getNumOfCommElementAll(){ return mvCommElementAll.size();}
    CCommElement* getCommElementAll(const uint& all_com_index){ return mvCommElementAll[all_com_index];}
    CElement* getElementAll(const uint& all_com_index){ return mvCommElementAll[all_com_index]->getElement();}
    uint getNumOfCommElement(){ return mvCommElement.size();}
    CCommElement* getCommElement(const uint& com_index){ return mvCommElement[com_index];}
    CElement* getElement(const uint& com_index){ return mvCommElement[com_index]->getElement();}
    uint getNumOfDNode(){ return mvDNode.size();}
    CNode* getDNode(const uint& dcom_index){ return mvDNode[dcom_index];}
    uint getNumOfDCommElement(){ return mvDCommElement.size();}
    CCommElement* getDCommElement(const uint& dcom_index){ return mvDCommElement[dcom_index];}
    CElement* getDElement(const uint& dcom_index){ return mvDElement[dcom_index];}
    uint getNumOfNode(){ return mvNode.size();}
    void reserveNode(const uint& res_size){ mvNode.reserve(res_size);}
    void setNode(CNode* pNode){ mvNode.push_back(pNode);}             
    CNode* getNode(const uint& com_index){ return mvNode[com_index];} 
    void setSendNode(CNode* pNode, const uint& comNodeID){ mvSendCommNodeID.push_back(comNodeID); mvSendNode.push_back(pNode);}    
    CNode* getSendNodeIX(const uint& index){ return mvSendNode[index];}
    uint getNumOfSendNode(){ return mvSendNode.size();}
    uint& getSendCommNodeID(const uint& index){ return mvSendCommNodeID[index];}
    void setRecvNode(CNode* pNode, const uint& comNodeID){ mvRecvCommNodeID.push_back(comNodeID); mvRecvNode.push_back(pNode);}    
    CNode* getRecvNodeIX(const uint& index){ return mvRecvNode[index];}
    uint getNumOfRecvNode(){ return mvRecvNode.size();}
    uint& getRecvCommNodeID(const uint& index){ return mvRecvCommNodeID[index];}
    void resizeNodeRank(const uint& res_size){ mvNodeRank.resize(res_size);}
    void setNodeRank(const uint& commNodeID, const uint& rank){ mvNodeRank[commNodeID]= rank;}
    uint& getNodeRank(const uint& commNodeID){ return  mvNodeRank[commNodeID];}
    void AllocateCommElement();
    void setupAggCommElement(vector<CElement*> &vElement);
    void sortCommNodeIndex();
    void setupMapID2CommID();
};
}
#endif	/* _COMMMesh_H */
