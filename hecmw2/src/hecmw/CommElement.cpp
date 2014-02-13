/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/CommElement.cpp
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
    delete[] mvbSend;
    delete[] mvbRecv;
    delete[] mvbOther;
    delete[] mvbNodeIXCheck;
    delete[] mvbDNodeMarking;
}
void CCommElement::setupProgNodeRank(const uiint& mgLevel)
{
    CEdgeTree *pEdgeTree = CEdgeTree::Instance();
    CEdgeFaceTree *pEdgeFaceTree = CEdgeFaceTree::Instance();
    uiint numOfEdge= mvEdgeRank.size();
    uiint numOfFace= mvFaceRank.size();
    uiint iedge,iface;
    uiint *localNodes, *localEdges;
    int   rankDiff;
    for(iedge=0; iedge< numOfEdge; iedge++) {
        localNodes= pEdgeTree->getHexaLocalNodeNum(iedge);
        rankDiff= mvNodeRank[localNodes[0]] - mvNodeRank[localNodes[1]];
        if(mgLevel==1) rankDiff *= -1;
        if(rankDiff <  0) mvEdgeRank[iedge]= mvNodeRank[localNodes[0]];
        if(rankDiff >= 0) mvEdgeRank[iedge]= mvNodeRank[localNodes[1]];
    };
    for(iface=0; iface< numOfFace; iface++) {
        localEdges= pEdgeFaceTree->getHexaFaceConnEdge(iface);
        rankDiff= mvEdgeRank[localEdges[0]] - mvEdgeRank[localEdges[2]];
        if(mgLevel==1) rankDiff *= -1;
        if(rankDiff <  0) mvFaceRank[iface]= mvEdgeRank[localEdges[0]];
        if(rankDiff >= 0) mvFaceRank[iface]= mvEdgeRank[localEdges[2]];
    };
    if(numOfFace > 4) {
        rankDiff= mvFaceRank[0] - mvFaceRank[1];
        if(mgLevel==1) rankDiff *= -1;
        if(rankDiff < 0) mVolRank= mvFaceRank[0];
        if(rankDiff >=0) mVolRank= mvFaceRank[1];
    } else if(numOfFace == 1) {
        mVolRank= mvFaceRank[0];
    } else {
        mVolRank= mvEdgeRank[0];
    }
}
void CCommElement::sortNodeRank(const uiint& myRank, const uiint& transRank)
{
    mvSendNode.clear();
    mvRecvNode.clear();
    mvOtherNode.clear();
    CNode* pNode;
    uiint numOfVert(mpElement->getNumOfNode());
    uiint ivert;
    for(ivert=0; ivert< numOfVert; ivert++) {
        pNode= mpElement->getNode(ivert);
        if(mvNodeRank[ivert]==myRank) {
            mvbSend[ivert]=true;
            mvSendNode.push_back(pNode);
        } else if(mvNodeRank[ivert]==transRank) {
            mvbRecv[ivert]=true;
            mvRecvNode.push_back(pNode);
        } else {
            mvbOther[ivert]=true;
            mvOtherNode.push_back(pNode);
        }
    };
    if(mvOtherNode.size()==0 && mvSendNode.size() > 0 && mvRecvNode.size() > 0) {
        mbCommunication= true;
    }
    if(mvSendNode.size() == numOfVert) {
        mbRCommunication= true;
    }
}
void CCommElement::setNeibCommNodeIndex(const uiint& ivert, const uiint& comID)
{
    if(!mvbNodeIXCheck[ivert]) {
        mvCommNodeIndex[ivert]= comID;
        mvbNodeIXCheck[ivert]= true;
    }
}
void CCommElement::setCommNodeIndex(const uiint& ivert, uiint& comID, vector<CNode*>& vNode, vuint& vCommMeshNodeRank)
{
    CNode* pNode;
    if(!mvbNodeIXCheck[ivert]) {
        mvCommNodeIndex[ivert]= comID;
        mvbNodeIXCheck[ivert] = true;
        pNode = mpElement->getNode(ivert);
        vNode.push_back(pNode);
        vCommMeshNodeRank.push_back(mvNodeRank[ivert]);
        vector<CCommElement*> vNeibCommElem= mvvAggCommElem[ivert];
        vuint vNeibVert = mvvNeibCommElemVert[ivert];
        CCommElement* neibCommElem;
        uiint neibVert;
        uiint numOfAgg= vNeibCommElem.size();
        uiint iagg;
        for(iagg=0; iagg< numOfAgg; iagg++) {
            neibCommElem= vNeibCommElem[iagg];
            neibVert= vNeibVert[iagg];
            neibCommElem->setNeibCommNodeIndex(neibVert, comID);
        };
        comID += 1;
    }
}
void CCommElement::getDNode(const uiint& ivert, vector<CNode*>& vDNode)
{
    CNode* pDNode;
    if(!mbCommunication) {
        if(!mvbDNodeMarking[ivert]) {
            vector<CCommElement*> vNeibCommElem= mvvAggCommElem[ivert];
            vuint vNeibVert = mvvNeibCommElemVert[ivert];
            CCommElement* neibCommElem;
            uiint neibVert;
            uiint numOfAgg= vNeibCommElem.size();
            uiint iagg;
            bool bCommElem(false);
            for(iagg=0; iagg< numOfAgg; iagg++) {
                neibCommElem= vNeibCommElem[iagg];
                if(neibCommElem->isCommElement()) {
                    bCommElem=true;
                }
            };
            if(!bCommElem) {
                pDNode = mpElement->getNode(ivert);
                markingDNode(ivert);
                vDNode.push_back(pDNode);
                for(iagg=0; iagg< numOfAgg; iagg++) {
                    neibCommElem= vNeibCommElem[iagg];
                    neibVert= vNeibVert[iagg];
                    if(!neibCommElem->isMarkingDNode(neibVert)) {
                        neibCommElem->markingDNode(neibVert);
                    }
                };
            }
        }
    }
}
void CCommElement::markingDNode(const uiint& ivert)
{
    mvbDNodeMarking[ivert]=true;
}
