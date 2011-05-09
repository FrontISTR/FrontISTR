/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   CommMesh2.cxx
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
#include "CommMesh2.h"
using namespace pmw;
CCommMesh2::CCommMesh2()
{
    ;
}
CCommMesh2::~CCommMesh2()
{
    uint numOfNode= mvCommNode.size();
    CCommNode *pCommNode;
    uint inode;
    for(inode=0; inode< numOfNode; inode++){
        pCommNode= mvCommNode[inode];
        if(pCommNode->getLevel()==mMGLevel){
            delete pCommNode;
        }
    };
    uint numOfFace= mvCommFace.size();
    CCommFace *pCommFace;
    uint iface;
    for(iface=0; iface< numOfFace; iface++){
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
CCommNode* CCommMesh2::getCommVertNode(const uint& id)
{
    uint index;
    index= mmCommNodeID2Index[id];
    return mvCommNode[index];
}
CCommFace* CCommMesh2::getCommFace(const uint& id)
{
    uint index;
    index= mmCommFaceID2Index[id];
    return mvCommFace[index];
}
void CCommMesh2::setupVertCommNode(CCommMesh2* pProgCommMesh)
{
    uint numOfNode;
    CCommNode *pCommNode;
    numOfNode= mvCommNode.size();
    pProgCommMesh->reserveCommNode(numOfNode);
    uint inode;
    for(inode=0; inode< numOfNode; inode++){
        pCommNode= mvCommNode[inode];
        pProgCommMesh->addCommNode(pCommNode);
    };
}
void CCommMesh2::setupAggFace()
{
    CCommFace *pCommFace;
    CCommNode *pCommNode;
    uint numOfCommFace(mvCommFace.size());
    uint numOfCommNode;
    uint iface,inode;
    numOfCommNode= mvCommNode.size();
    for(inode=0; inode< numOfCommNode; inode++){
        pCommNode= mvCommNode[inode];
        pCommNode->clearAggElemID();
        pCommNode->clearNeibElemVert();
    };
    for(iface=0; iface< numOfCommFace; iface++){
        pCommFace= mvCommFace[iface];
        numOfCommNode= pCommFace->getVertCommNodeSize();
        for(inode=0; inode< numOfCommNode; inode++){
            pCommNode= pCommFace->getVertCommNode(inode);
            pCommNode->setAggElemID(pCommFace->getID());
            pCommNode->setNeibElemVert(pCommFace->getID(), inode);
        };
    };
}
void CCommMesh2::setupEdgeCommNode(CCommMesh2 *pProgCommMesh)
{
    PairCommNode pairCommNode;
    CCommFace *pFace, *pNeibFace;
    CCommNode* pEdgeCommNode;
    vdouble vCoord; vCoord.resize(3);
    uint numOfFace, numOfEdge;
    uint level;
    uint iface,iedge;
    uint countID= mvCommNode.size();
    numOfFace= mvCommFace.size();
    for(iface=0; iface< numOfFace; iface++){
        pFace= mvCommFace[iface];
        numOfEdge= pFace->getNumOfEdge();
        for(iedge=0; iedge< numOfEdge; iedge++){
            if(!pFace->isEdgeNodeMarking(iedge)){
                pairCommNode= pFace->getEdgePairCommNode(iedge);
                uint numOfAggFaceA= pairCommNode.first->getNumOfAggElem();
                uint numOfAggFaceB= pairCommNode.second->getNumOfAggElem();
                uint faceIDa,faceIDb;
                uint iagg,jagg;
                bool bfind(false);
                for(iagg=0; iagg< numOfAggFaceA; iagg++){
                    faceIDa= pairCommNode.first->getAggElemID(iagg);
                    for(jagg=0; jagg< numOfAggFaceB; jagg++){
                        faceIDb= pairCommNode.second->getAggElemID(jagg);
                        if(faceIDa==faceIDb && faceIDa!=pFace->getID()){
                            uint index= mmCommFaceID2Index[faceIDa];
                            pNeibFace= mvCommFace[index];
                            pEdgeCommNode= new CCommNode;
                            vCoord[0]= pairCommNode.first->getX(); vCoord[0]+= pairCommNode.second->getX(); vCoord[0] *= 0.5;
                            vCoord[1]= pairCommNode.first->getY(); vCoord[1]+= pairCommNode.second->getY(); vCoord[1] *= 0.5;
                            vCoord[2]= pairCommNode.first->getZ(); vCoord[2]+= pairCommNode.second->getZ(); vCoord[2] *= 0.5;
                            pEdgeCommNode->setCoord(vCoord);
                            pEdgeCommNode->setLevel(mMGLevel+1);
                            pEdgeCommNode->setID(countID);
                            countID++;
                            pProgCommMesh->addCommNode(pEdgeCommNode);
                            mvEdgeCommNode.push_back(pEdgeCommNode);  
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
                if(!bfind){
                    pEdgeCommNode= new CCommNode;
                    vCoord[0]= pairCommNode.first->getX(); vCoord[0]+= pairCommNode.second->getX(); vCoord[0] *= 0.5;
                    vCoord[1]= pairCommNode.first->getY(); vCoord[1]+= pairCommNode.second->getY(); vCoord[1] *= 0.5;
                    vCoord[2]= pairCommNode.first->getZ(); vCoord[2]+= pairCommNode.second->getZ(); vCoord[2] *= 0.5;
                    pEdgeCommNode->setCoord(vCoord);
                    pEdgeCommNode->setLevel(mMGLevel+1);
                    pEdgeCommNode->setID(countID);
                    countID++;
                    pProgCommMesh->addCommNode(pEdgeCommNode);
                    mvEdgeCommNode.push_back(pEdgeCommNode); 
                    pFace->setEdgeCommNode(pEdgeCommNode, iedge);
                    pFace->markingEdgeNode(iedge);
                }
            }
        };
    };
}
void CCommMesh2::setupFaceCommNode(CCommMesh2 *pProgCommMesh)
{
    uint level;
    CCommNode *pFaceCommNode;
    CCommFace *pFace;
    uint numOfFace= mvCommFace.size();
    uint iface;
    uint countID= mvCommNode.size()+mvEdgeCommNode.size();
    for(iface=0; iface< numOfFace; iface++){
        pFace= mvCommFace[iface];
        if(pFace->getNumOfEdge() > 2){
            pFaceCommNode= new CCommNode;
            pFace->setFaceCommNode(pFaceCommNode);
            vdouble vCoord; vCoord.resize(3);
            uint numOfVert= pFace->getVertCommNodeSize();
            CCommNode *pVertCommNode;
            vCoord[0]=0.0; vCoord[1]=0.0; vCoord[2]=0.0;
            uint ivert;
            for(ivert=0; ivert< numOfVert; ivert++){
                pVertCommNode= pFace->getVertCommNode(ivert);
                vCoord[0] += pVertCommNode->getX();  vCoord[1] += pVertCommNode->getY();  vCoord[2] += pVertCommNode->getZ();
            };
            vCoord[0] /= numOfVert; vCoord[1] /= numOfVert; vCoord[2] /= numOfVert;
            pFaceCommNode->setLevel(mMGLevel+1);
            pFaceCommNode->setCoord(vCoord);
            pFaceCommNode->setID(countID);
            countID++;
            pProgCommMesh->addCommNode(pFaceCommNode);
            mvFaceCommNode.push_back(pFaceCommNode);
        }
    };
}
