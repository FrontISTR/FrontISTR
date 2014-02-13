/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/CommMesh2.cpp
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
#include "CommMesh2.h"
using namespace pmw;
CCommMesh2::CCommMesh2()
{
    ;
}
CCommMesh2::~CCommMesh2()
{
    uiint numOfNode= mvCommNode.size();
    CCommNode *pCommNode;
    uiint inode;
    for(inode=0; inode< numOfNode; inode++) {
        pCommNode= mvCommNode[inode];
        if(pCommNode->getLevel()==mMGLevel) {
            delete pCommNode;
        }
    };
    uiint numOfFace= mvCommFace.size();
    CCommFace *pCommFace;
    uiint iface;
    for(iface=0; iface< numOfFace; iface++) {
        pCommFace= mvCommFace[iface];
        delete pCommFace;
    };
}
void CCommMesh2::addCommNode(CCommNode* pCommNode)
{
    mvCommNode.push_back(pCommNode);
    mmCommNodeID2Index[pCommNode->getID()]= mvCommNode.size()-1;
}
void CCommMesh2::addCommFace(CCommFace* pCommFace)
{
    mvCommFace.push_back(pCommFace);
    mmCommFaceID2Index[pCommFace->getID()]= mvCommFace.size()-1;
}
CCommNode* CCommMesh2::getCommNode(const uiint& id)
{
    uiint index;
    index= mmCommNodeID2Index[id];
    return mvCommNode[index];
}
CCommFace* CCommMesh2::getCommFace(const uiint& id)
{
    uiint index;
    index= mmCommFaceID2Index[id];
    return mvCommFace[index];
}
void CCommMesh2::setupCommNode(CCommMesh2* pProgCommMesh)
{
    uiint numOfNode;
    CCommNode *pCommNode;
    numOfNode= mvCommNode.size();
    pProgCommMesh->reserveCommNode(numOfNode);
    uiint inode;
    for(inode=0; inode< numOfNode; inode++) {
        pCommNode= mvCommNode[inode];
        pProgCommMesh->addCommNode(pCommNode);
    };
}
void CCommMesh2::setupAggFace()
{
    CCommFace *pCommFace;
    CCommNode *pCommNode;
    uiint numOfCommFace(mvCommFace.size());
    uiint numOfCommNode;
    uiint iface,inode;
    numOfCommNode= mvCommNode.size();
    for(inode=0; inode< numOfCommNode; inode++) {
        pCommNode= mvCommNode[inode];
        pCommNode->clearAggElemID();
        pCommNode->clearNeibElemVert();
    };
    for(iface=0; iface< numOfCommFace; iface++) {
        pCommFace= mvCommFace[iface];
        uiint nNumOfVert= pCommFace->getNumOfVert();
        for(inode=0; inode< nNumOfVert; inode++) {
            pCommNode= pCommFace->getCommNode(inode);
            pCommNode->setAggElemID(pCommFace->getID());
            pCommNode->setNeibElemVert(pCommFace->getID(), inode);
        };
    };
}
void CCommMesh2::setupEdgeCommNode(CCommMesh2 *pProgCommMesh, const uiint& nLevel)
{
    PairCommNode pairCommNode;
    CCommFace *pFace, *pNeibFace;
    CCommNode* pEdgeCommNode;
    vdouble vCoord;
    vCoord.resize(3);
    uiint numOfFace, numOfEdge;
    uiint iface,iedge;
    uiint countID= mvCommNode.size();
    mEdgeCount= 0;
    numOfFace= mvCommFace.size();
    for(iface=0; iface< numOfFace; iface++) {
        pFace= mvCommFace[iface];
        bool bExec(false);
        if(nLevel > 0 && pFace->getOrder()==ElementOrder::Second) bExec=true;
        if(pFace->getOrder()==ElementOrder::First) bExec=true;
        if(bExec) {
            numOfEdge= pFace->getNumOfEdge();
            for(iedge=0; iedge< numOfEdge; iedge++) {
                if(!pFace->isEdgeNodeMarking(iedge)) {
                    pairCommNode= pFace->getEdgePairCommNode(iedge);
                    uiint numOfAggFaceA= pairCommNode.first->getNumOfAggElem();
                    uiint numOfAggFaceB= pairCommNode.second->getNumOfAggElem();
                    uiint faceIDa,faceIDb;
                    uiint iagg,jagg;
                    bool bfind(false);
                    for(iagg=0; iagg< numOfAggFaceA; iagg++) {
                        faceIDa= pairCommNode.first->getAggElemID(iagg);
                        for(jagg=0; jagg< numOfAggFaceB; jagg++) {
                            faceIDb= pairCommNode.second->getAggElemID(jagg);
                            if(faceIDa==faceIDb && faceIDa!=pFace->getID()) {
                                uiint index= mmCommFaceID2Index[faceIDa];
                                pNeibFace= mvCommFace[index];
                                pEdgeCommNode= new CCommNode;//---------------------------生成
                                vCoord[0]= pairCommNode.first->getX();
                                vCoord[0]+= pairCommNode.second->getX();
                                vCoord[0] *= 0.5;
                                vCoord[1]= pairCommNode.first->getY();
                                vCoord[1]+= pairCommNode.second->getY();
                                vCoord[1] *= 0.5;
                                vCoord[2]= pairCommNode.first->getZ();
                                vCoord[2]+= pairCommNode.second->getZ();
                                vCoord[2] *= 0.5;
                                pEdgeCommNode->setCoord(vCoord);
                                pEdgeCommNode->setLevel(mMGLevel+1);
                                pEdgeCommNode->setID(countID);
                                countID++;
                                mEdgeCount++;
                                if(pProgCommMesh)
                                    pProgCommMesh->addCommNode(pEdgeCommNode);
                                if(pFace->getOrder()==ElementOrder::Second) {
                                    mvCommNode.push_back(pEdgeCommNode);
                                }
                                pFace->setEdgeCommFace(pNeibFace, iedge);
                                pFace->setEdgeCommNode(pEdgeCommNode, iedge);
                                pFace->markingEdgeNode(iedge);
                                pNeibFace->setEdgeCommFace(pFace, pairCommNode);
                                pNeibFace->setEdgeCommNode(pEdgeCommNode, pairCommNode);
                                pNeibFace->markingEdgeNode(pairCommNode);
                                bfind=true;
                                break;
                            }
                        };
                        if(bfind) break;
                    };
                    if(!bfind) {
                        pEdgeCommNode= new CCommNode;//----------------------------------生成
                        vCoord[0]= pairCommNode.first->getX();
                        vCoord[0]+= pairCommNode.second->getX();
                        vCoord[0] *= 0.5;
                        vCoord[1]= pairCommNode.first->getY();
                        vCoord[1]+= pairCommNode.second->getY();
                        vCoord[1] *= 0.5;
                        vCoord[2]= pairCommNode.first->getZ();
                        vCoord[2]+= pairCommNode.second->getZ();
                        vCoord[2] *= 0.5;
                        pEdgeCommNode->setCoord(vCoord);
                        pEdgeCommNode->setLevel(mMGLevel+1);
                        pEdgeCommNode->setID(countID);
                        countID++;
                        mEdgeCount++;
                        if(pProgCommMesh)
                            pProgCommMesh->addCommNode(pEdgeCommNode);
                        if(pFace->getOrder()==ElementOrder::Second) {
                            mvCommNode.push_back(pEdgeCommNode);
                        }
                        pFace->setEdgeCommNode(pEdgeCommNode, iedge);
                        pFace->markingEdgeNode(iedge);
                    }
                }
            };
            pFace->replaceEdgeCommNode();
        }
    };
}
void CCommMesh2::setupFaceCommNode(CCommMesh2 *pProgCommMesh)
{
    uiint level;
    CCommNode *pFaceCommNode;
    CCommFace *pFace;
    uiint numOfFace= mvCommFace.size();
    uiint iface;
    uiint countID= mvCommNode.size() + mEdgeCount;
    for(iface=0; iface< numOfFace; iface++) {
        pFace= mvCommFace[iface];
        if(pFace->getNumOfEdge() > 2) {
            pFaceCommNode= new CCommNode;//-------------------------生成
            pFace->setFaceCommNode(pFaceCommNode);
            vdouble vCoord;
            vCoord.resize(3);
            uiint numOfVert= pFace->getNumOfVert();
            CCommNode *pVertCommNode;
            vCoord[0]=0.0;
            vCoord[1]=0.0;
            vCoord[2]=0.0;
            uiint ivert;
            level=0;
            for(ivert=0; ivert< numOfVert; ivert++) {
                pVertCommNode= pFace->getCommNode(ivert);
                vCoord[0] += pVertCommNode->getX();
                vCoord[1] += pVertCommNode->getY();
                vCoord[2] += pVertCommNode->getZ();
            };
            vCoord[0] /= numOfVert;
            vCoord[1] /= numOfVert;
            vCoord[2] /= numOfVert;
            pFaceCommNode->setLevel(mMGLevel+1);
            pFaceCommNode->setCoord(vCoord);
            pFaceCommNode->setID(countID);
            countID++;
            pProgCommMesh->addCommNode(pFaceCommNode);
        }
    };
}
void CCommMesh2::deleteProgData()
{
    uiint iface, nNumOfFace = mvCommFace.size();
    for(iface=0; iface < nNumOfFace; iface++) mvCommFace[iface]->deleteProgData();
}
