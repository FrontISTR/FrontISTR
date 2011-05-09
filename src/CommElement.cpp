/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   CommElement.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include <vector>
#include "Element.h"
#include "CommElement.h"
using namespace pmw;
CCommElement::CCommElement()
{
    mpLogger = Utility::CLogger::Instance();
    mbCommunication = false;
    mbRCommunication= false;
}
CCommElement::~CCommElement()
{
    ;
}
void CCommElement::setupProgNodeRank(const uint& mgLevel)
{
    CEdgeTree *pEdgeTree = CEdgeTree::Instance();
    CEdgeFaceTree *pEdgeFaceTree = CEdgeFaceTree::Instance();
    uint numOfEdge= mvEdgeRank.size();
    uint numOfFace= mvFaceRank.size();
    uint iedge,iface;
    uint *localNodes, *localEdges;
    int   rankDiff;
    for(iedge=0; iedge< numOfEdge; iedge++){
        localNodes= pEdgeTree->getHexaLocalNodeNum(iedge);
        rankDiff= mvNodeRank[localNodes[0]] - mvNodeRank[localNodes[1]];
        if(mgLevel==1) rankDiff *= -1; 
        if(rankDiff <  0) mvEdgeRank[iedge]= mvNodeRank[localNodes[0]];
        if(rankDiff >= 0) mvEdgeRank[iedge]= mvNodeRank[localNodes[1]];
    };
    for(iface=0; iface< numOfFace; iface++){
        localEdges= pEdgeFaceTree->getHexaFaceConnEdge(iface);
        rankDiff= mvEdgeRank[localEdges[0]] - mvEdgeRank[localEdges[2]];
        if(mgLevel==1) rankDiff *= -1;
        if(rankDiff <  0) mvFaceRank[iface]= mvEdgeRank[localEdges[0]];
        if(rankDiff >= 0) mvFaceRank[iface]= mvEdgeRank[localEdges[2]];
    };
    if(numOfFace > 4){
        rankDiff= mvFaceRank[0] - mvFaceRank[1];
        if(mgLevel==1) rankDiff *= -1;
        if(rankDiff < 0) mVolRank= mvFaceRank[0];
        if(rankDiff >=0) mVolRank= mvFaceRank[1];
    }else if(numOfFace == 1){
        mVolRank= mvFaceRank[0];
    }else{
        mVolRank= mvEdgeRank[0];
    }
}
void CCommElement::sortNodeRank(const uint& myRank, const uint& transRank)
{
    mvSendNode.clear(); mvRecvNode.clear(); mvOtherNode.clear();
    CNode* pNode;
    uint numOfVert(mpElement->getNumOfNode());
    uint ivert;
    for(ivert=0; ivert< numOfVert; ivert++){
        pNode= mpElement->getNode(ivert);
        if(mvNodeRank[ivert]==myRank){
            mvbSend[ivert]=true;
            mvSendNode.push_back(pNode);
        }
        else if(mvNodeRank[ivert]==transRank){
            mvbRecv[ivert]=true;
            mvRecvNode.push_back(pNode);
        }
        else{
            mvbOther[ivert]=true;
            mvOtherNode.push_back(pNode);
        }
    };
    if(mvOtherNode.size()==0 && mvSendNode.size() > 0 && mvRecvNode.size() > 0){
        mbCommunication= true;
    }
    if(mvSendNode.size() == numOfVert){
        mbRCommunication= true;
    }
}
void CCommElement::setNeibCommNodeIndex(const uint& ivert, const uint& comID)
{
    if(!mvbNodeIXCheck[ivert]){
        mvCommNodeIndex[ivert]= comID;
        mvbNodeIXCheck[ivert]= true;
    }
}
void CCommElement::setCommNodeIndex(const uint& ivert, uint& comID, vector<CNode*>& vNode, vuint& vCommMeshNodeRank)
{
    CNode* pNode;
    if(!mvbNodeIXCheck[ivert]){
        mvCommNodeIndex[ivert]= comID;
        mvbNodeIXCheck[ivert] = true;
        pNode = mpElement->getNode(ivert);
        vNode.push_back(pNode);
        vCommMeshNodeRank.push_back(mvNodeRank[ivert]);
        vector<CCommElement*> vNeibCommElem= mvvAggCommElem[ivert];
        vuint vNeibVert = mvvNeibCommElemVert[ivert];
        CCommElement* neibCommElem;
        uint neibVert;
        uint numOfAgg= vNeibCommElem.size();
        uint iagg;
        for(iagg=0; iagg< numOfAgg; iagg++){
            neibCommElem= vNeibCommElem[iagg];
            neibVert= vNeibVert[iagg];
            neibCommElem->setNeibCommNodeIndex(neibVert, comID);
        };
        comID += 1;
    }
}
void CCommElement::getDNode(const uint& ivert, vector<CNode*>& vDNode)
{
    CNode* pDNode;
    if(!mbCommunication){
    if(!mvbDNodeMarking[ivert]){
        vector<CCommElement*> vNeibCommElem= mvvAggCommElem[ivert];
        vuint vNeibVert = mvvNeibCommElemVert[ivert];
        CCommElement* neibCommElem;
        uint neibVert;
        uint numOfAgg= vNeibCommElem.size();
        uint iagg;
        bool bCommElem(false);
        for(iagg=0; iagg< numOfAgg; iagg++){
            neibCommElem= vNeibCommElem[iagg];
            if(neibCommElem->isCommElement()){
                bCommElem=true;
            }
        };
        if(!bCommElem){
            pDNode = mpElement->getNode(ivert);
            markingDNode(ivert);
            vDNode.push_back(pDNode);
            for(iagg=0; iagg< numOfAgg; iagg++){
                neibCommElem= vNeibCommElem[iagg];
                neibVert= vNeibVert[iagg];
                if(!neibCommElem->isMarkingDNode(neibVert)){
                    neibCommElem->markingDNode(neibVert);
                }
            };
        }
    }
    }
}
void CCommElement::markingDNode(const uint& ivert)
{
    mvbDNodeMarking[ivert]=true;
}
