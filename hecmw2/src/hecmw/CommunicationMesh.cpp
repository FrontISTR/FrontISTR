/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/CommunicationMesh.cpp
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
#include "CommunicationMesh.h"
using namespace pmw;
#include <iostream>
CCommMesh::CCommMesh(CIndexBucket* pBucket)
{
    mpBucket = pBucket;
    mpLogger= Utility::CLogger::Instance();
}
CCommMesh::~CCommMesh()
{
    for_each(mvCommElementAll.begin(), mvCommElementAll.end(), DeleteObject());
    std::cout << "~CCommMesh" << std::endl;
}
void CCommMesh::AllocateCommElement()
{
    uiint numOfCommElemAll= mvCommElementAll.size();
    CCommElement* pCommElem;
    CElement* pElem;
    uiint icom;
    for(icom=0; icom< numOfCommElemAll; icom++) {
        pCommElem= mvCommElementAll[icom];
        pCommElem->sortNodeRank(mRankID, mTransmitRankID);
        pElem= pCommElem->getElement();
        if(pCommElem->isCommElement()) {
            mvCommElement.push_back(pCommElem);
            pElem->interCommElem();
        } else if(!pCommElem->isRCommElement()) {
            mvDCommElement.push_back(pCommElem);
            pElem->interDCommElem();
        } else {
            pElem->interRCommElem();
        }
        pCommElem->setID(icom);
        pElem->setCommID(icom);
    };
}
void CCommMesh::setupAggCommElement(vector<CElement*> &vElement)
{
    uiint numOfCommElem= mvCommElementAll.size();
    uiint numOfVert;
    uiint numOfAggElem;
    CElement* pElem;
    CCommElement* pCommElem;
    CCommElement* pNeibCommElem;
    uiint neibComID, neibVert;
    CNode* pNode;
    uiint icome,ivert,iagg;
    uiint elemID,elemIndex;
    for(icome=0; icome< numOfCommElem; icome++) {
        pCommElem= mvCommElementAll[icome];
        numOfVert= pCommElem->getNumOfVert();
        for(ivert=0; ivert< numOfVert; ivert++) {
            pNode= pCommElem->getNode(ivert);
            numOfAggElem= pNode->getNumOfAggElem();
            for(iagg=0; iagg< numOfAggElem; iagg++) {
                elemID= pNode->getAggElemID(iagg);
                elemIndex= mpBucket->getIndexElement(elemID);
                pElem= vElement[elemIndex];
                if(pElem->isInterCommElem()||pElem->isInterDCommElem()||pElem->isInterRCommElem()) {
                    neibComID= pElem->getCommID();
                    pNeibCommElem= mvCommElementAll[neibComID];
                    pCommElem->setAggCommElement(ivert, pNeibCommElem);
                    neibVert= pNode->getNeibElemIDVert(elemID);
                    pCommElem->setNeibCommElemVert(ivert,neibVert);
                }
            };
        };
    };
}
void CCommMesh::sortCommNodeIndex()
{
    uiint numOfCommElem= mvCommElementAll.size();
    uiint numOfVert;
    CCommElement* pCommElem;
    uiint icome, ivert;
    uiint comNodeIndex(0);
    for(icome=0; icome< numOfCommElem; icome++) {
        pCommElem= mvCommElementAll[icome];
        numOfVert= pCommElem->getNumOfVert();
        for(ivert=0; ivert< numOfVert; ivert++) {
            pCommElem->setCommNodeIndex(ivert, comNodeIndex, mvNode, mvNodeRank);
        };
    };
    mpLogger->Info(Utility::LoggerMode::Debug, "CommMesh::sortCommNodeIndex mvNode.size => ", (uiint)mvNode.size());
    bool* vIndexCheck;
    uiint numOfNode= mvNode.size();
    vIndexCheck = new bool[numOfNode];
    for(uiint i=0; i< numOfNode; i++) {
        vIndexCheck[i] = false;
    }
    uiint nIndex, rank;
    numOfCommElem= mvCommElement.size();
    for(icome=0; icome< numOfCommElem; icome++) {
        pCommElem= mvCommElement[icome];
        numOfVert= pCommElem->getNumOfVert();
        for(ivert=0; ivert< numOfVert; ivert++) {
            nIndex= pCommElem->getCommNodeIndex(ivert);
            if(!vIndexCheck[nIndex]) {
                rank= mvNodeRank[nIndex];
                if(mRankID==rank) {
                    mvSendNode.push_back(mvNode[nIndex]);
                    mvSendCommNodeID.push_back(nIndex);
                }
                if(mTransmitRankID==rank) {
                    mvRecvNode.push_back(mvNode[nIndex]);
                    mvRecvCommNodeID.push_back(nIndex);
                }
                vIndexCheck[nIndex]= true;
            }
        };
    };
    delete []vIndexCheck;
    uiint numOfDCommElem= mvDCommElement.size();
    CCommElement* pDCommElem;
    CElement* pDElem;
    for(icome=0; icome< numOfDCommElem; icome++) {
        pDCommElem= mvDCommElement[icome];
        pDElem = pDCommElem->getElement();
        mvDElement.push_back(pDElem);
    };
    mpLogger->Info(Utility::LoggerMode::Debug, "CommMesh::sortCommNodeIndex mvDElement.size => ", (uiint)mvDElement.size());
    for(icome=0; icome< numOfDCommElem; icome++) {
        pDCommElem= mvDCommElement[icome];
        numOfVert= pDCommElem->getNumOfVert();
        for(ivert=0; ivert< numOfVert; ivert++) {
            pDCommElem->getDNode(ivert, mvDNode);
        };
    };
    mpLogger->Info(Utility::LoggerMode::Debug, "CommMesh::sortCommNodeIndex mvDNode.size => ", (uiint)mvDNode.size());
    uiint maxIndex;
    if(mvDNode.size() > 0) {
        sortID<CNode*>(mvDNode, mvDNode.size());
    }
    if(mvDElement.size() > 0) {
        sortID<CElement*>(mvDElement, mvDElement.size());
    }
}
void CCommMesh::setupMapID2CommID()
{
    uiint numOfCommElem= mvCommElementAll.size();
    uiint numOfCommNode= mvNode.size();
    CNode* pNode;
    CElement* pElem;
    CCommElement* pCommElem;
    uiint index;
    for(index=0; index< numOfCommNode; index++) {
        pNode= mvNode[index];
        mmCommNodeIX[pNode->getID()]= index;
    };
    for(index=0; index< numOfCommElem; index++) {
        pCommElem= mvCommElementAll[index];
        pElem= pCommElem->getElement();
        mmCommElementIX[pElem->getID()]= index;
    };
}
