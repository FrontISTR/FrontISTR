/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   CommunicationMesh.cxx
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
    uint numOfCommElemAll= mvCommElementAll.size();
    CCommElement* pCommElem;
    CElement* pElem;
    uint icom;
    for(icom=0; icom< numOfCommElemAll; icom++){
        pCommElem= mvCommElementAll[icom];
        pCommElem->sortNodeRank(mRankID, mTransmitRankID);
        pElem= pCommElem->getElement();
        if(pCommElem->isCommElement()){
            mvCommElement.push_back(pCommElem);
            pElem->interCommElem(); 
        }else if(!pCommElem->isRCommElement()){
            mvDCommElement.push_back(pCommElem);
            pElem->interDCommElem();
        }else{
            pElem->interRCommElem();
        }
        pCommElem->setID(icom);
        pElem->setCommID(icom);
    };
}
void CCommMesh::setupAggCommElement(vector<CElement*> &vElement)
{
    uint numOfCommElem= mvCommElementAll.size();
    uint numOfVert;
    uint numOfAggElem;
    CElement* pElem;
    CCommElement* pCommElem;
    CCommElement* pNeibCommElem;
    uint neibComID, neibVert;
    CNode* pNode;
    uint icome,ivert,iagg;
    uint elemID,elemIndex;
    for(icome=0; icome< numOfCommElem; icome++){
        pCommElem= mvCommElementAll[icome];
        numOfVert= pCommElem->getNumOfVert();
        for(ivert=0; ivert< numOfVert; ivert++){
            pNode= pCommElem->getNode(ivert);
            numOfAggElem= pNode->getNumOfAggElem();
            for(iagg=0; iagg< numOfAggElem; iagg++){
                elemID= pNode->getAggElemID(iagg);
                elemIndex= mpBucket->getIndexElement(elemID);
                pElem= vElement[elemIndex];
                if(pElem->isInterCommElem()||pElem->isInterDCommElem()||pElem->isInterRCommElem()){
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
    uint numOfCommElem= mvCommElementAll.size();
    uint numOfVert;
    CCommElement* pCommElem;
    uint icome, ivert;
    uint comNodeIndex(0);
    for(icome=0; icome< numOfCommElem; icome++){
        pCommElem= mvCommElementAll[icome];
        numOfVert= pCommElem->getNumOfVert();
        for(ivert=0; ivert< numOfVert; ivert++){
            pCommElem->setCommNodeIndex(ivert, comNodeIndex, mvNode, mvNodeRank);
        };
    };
    mpLogger->Info(Utility::LoggerMode::Debug, "CommMesh::sortCommNodeIndex mvNode.size => ", (uint)mvNode.size());
    vbool vIndexCheck;  uint numOfNode= mvNode.size();
    vIndexCheck.reserve(numOfNode);
    for(uint i=0; i< numOfNode; i++){
        vIndexCheck.push_back(false);
    }
    uint nIndex, rank;
    numOfCommElem= mvCommElement.size();
    for(icome=0; icome< numOfCommElem; icome++){
        pCommElem= mvCommElement[icome];
        numOfVert= pCommElem->getNumOfVert();
        for(ivert=0; ivert< numOfVert; ivert++){
            nIndex= pCommElem->getCommNodeIndex(ivert);
            if(!vIndexCheck[nIndex]){
                rank= mvNodeRank[nIndex];
                if(mRankID==rank){
                    mvSendNode.push_back(mvNode[nIndex]);
                    mvSendCommNodeID.push_back(nIndex);
                }
                if(mTransmitRankID==rank){
                    mvRecvNode.push_back(mvNode[nIndex]);
                    mvRecvCommNodeID.push_back(nIndex);
                }
                vIndexCheck[nIndex]= true;
            }
        };
    };
    uint numOfDCommElem= mvDCommElement.size();
    CCommElement* pDCommElem;
    CElement* pDElem;
    for(icome=0; icome< numOfDCommElem; icome++){
        pDCommElem= mvDCommElement[icome];
        pDElem = pDCommElem->getElement();
        mvDElement.push_back(pDElem);
    };
    mpLogger->Info(Utility::LoggerMode::Debug, "CommMesh::sortCommNodeIndex mvDElement.size => ", (uint)mvDElement.size());
    for(icome=0; icome< numOfDCommElem; icome++){
        pDCommElem= mvDCommElement[icome];
        numOfVert= pDCommElem->getNumOfVert();
        for(ivert=0; ivert< numOfVert; ivert++){
            pDCommElem->getDNode(ivert, mvDNode);
        };
    };
    mpLogger->Info(Utility::LoggerMode::Debug, "CommMesh::sortCommNodeIndex mvDNode.size => ", (uint)mvDNode.size());
    uint maxIndex;
    if(mvDNode.size() > 0){
        sortID<CNode*>(mvDNode, mvDNode.size());
    }
    if(mvDElement.size() > 0){
        sortID<CElement*>(mvDElement, mvDElement.size());
    }
}
void CCommMesh::setupMapID2CommID()
{
    uint numOfCommElem= mvCommElementAll.size();
    uint numOfCommNode= mvNode.size();
    CNode* pNode;
    CElement* pElem;
    CCommElement* pCommElem;
    uint index;
    for(index=0; index< numOfCommNode; index++){
        pNode= mvNode[index];
        mmCommNodeIX[pNode->getID()]= index;
    };
    for(index=0; index< numOfCommElem; index++){
        pCommElem= mvCommElementAll[index];
        pElem= pCommElem->getElement();
        mmCommElementIX[pElem->getID()]= index;
    };
}
