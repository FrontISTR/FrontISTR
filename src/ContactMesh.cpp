/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   ContactMesh.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "SkinFace.h"
#include <vector>
#include "Vertex.h"
#include "ContactNode.h"
#include "OctreeKnot.h"
#include "ContactMesh.h"
using namespace pmw;
CContactMesh::CContactMesh()
{
    ;
}
CContactMesh::~CContactMesh()
{
	vector<CContactNode*>::iterator itConNode;
    CContactNode *pConNode;
	uint i;
	for(i=0; i < mvConNode.size() ; i++){
        pConNode = mvConNode[i];
        if(mLevel==pConNode->getLevel()) delete pConNode;
		mvConNode.erase(mvConNode.begin()+i);
    };
    for_each(mvFace.begin(), mvFace.end(), DeleteObject());
    for_each(mvSlaveFace.begin(), mvSlaveFace.end(), DeleteObject());
    vector<COctreeKnot*> vKnot;
    uint numOfLayer= mvKnot.size();
    uint ilayer;
    for(ilayer=1; ilayer < numOfLayer; ilayer++){
        vKnot= mvKnot[ilayer];
        for_each(vKnot.begin(), vKnot.end(), DeleteObject());
    };
    std::cout  << "~CContactMesh" << ", mgLevel= " << mLevel << std::endl;
}
void CContactMesh::addConNode(CContactNode* pConNode, const uint& id)
{
    mvConNode.push_back(pConNode);
    mmConNodeID2Index[id]= mvConNode.size()-1;
}
void CContactMesh::addMasterConNode(CContactNode* pConNode, const uint& id)
{
    mvMasterConNode.push_back(pConNode);
    mmMasterConNodeID2Index[id]= mvMasterConNode.size()-1;
}
void CContactMesh::addSlaveConNode(CContactNode* pConNode, const uint& id)
{
    mvSlaveConNode.push_back(pConNode);
    mmSlaveConNodeID2Index[id]= mvSlaveConNode.size()-1;
}
void CContactMesh::addMasterMeshID(const uint& id)
{
    mvMasterMeshID.push_back(id);
    uint index= mvMasterMeshID.size()-1;
    mmMasterMeshID2Index[id]= index;
}
void CContactMesh::addSlaveMeshID(const uint& id)
{
    mvSlaveMeshID.push_back(id);
    uint index= mvSlaveMeshID.size()-1;
    mmSlaveMeshID2Index[id]= index;
}
void CContactMesh::setupCoarseConNode(CContactMesh* pProgConMesh)
{
    uint numOfNode;
    numOfNode= mvConNode.size();
    CContactNode *pConNode;
    uint inode, id;
    for(inode=0; inode< numOfNode; inode++){
        pConNode= mvConNode[inode];
		pConNode->pushLevelMarking();
        id= pConNode->getID();
        pProgConMesh->addConNode(pConNode, id);
    };
    numOfNode= mvMasterConNode.size();
    for(inode=0; inode< numOfNode; inode++){
        pConNode= mvMasterConNode[inode];
        id= pConNode->getID();
        pProgConMesh->addMasterConNode(pConNode, id);
    };
    numOfNode= mvSlaveConNode.size();
    for(inode=0; inode< numOfNode; inode++){
        pConNode= mvSlaveConNode[inode];
        id= pConNode->getID();
        pProgConMesh->addSlaveConNode(pConNode, id);
    };
}
void CContactMesh::setupAggSkinFace()
{
    CSkinFace  *pMFace;
    CSkinFace  *pSFace;
    CContactNode *pConNode;
    uint numOfMFace(mvFace.size()), numOfSFace(mvSlaveFace.size()), numOfConNode;
    uint iface,icnode;
    numOfConNode= mvMasterConNode.size();
    for(icnode=0; icnode< numOfConNode; icnode++){
        pConNode= mvMasterConNode[icnode];
        pConNode->clearAggElemID();
        pConNode->clearNeibElemVert();
    };
    numOfConNode= mvSlaveConNode.size();
    for(icnode=0; icnode< numOfConNode; icnode++){
        pConNode= mvSlaveConNode[icnode];
        pConNode->clearAggElemID();
        pConNode->clearNeibElemVert();
    };
    for(iface=0; iface< numOfMFace; iface++){
        pMFace= mvFace[iface];
        numOfConNode= pMFace->getNumOfNode();
        for(icnode=0; icnode< numOfConNode; icnode++){
            pConNode= pMFace->getNode(icnode);
            pConNode->setAggElemID(pMFace->getID());
            pConNode->setNeibElemVert(pMFace->getID(), icnode);
        };
    };
    for(iface=0; iface< numOfSFace; iface++){
        pSFace= mvSlaveFace[iface];
        numOfConNode= pSFace->getNumOfNode();
        for(icnode=0; icnode< numOfConNode; icnode++){
            pConNode= pSFace->getNode(icnode);
            pConNode->setAggElemID(pSFace->getID());
            pConNode->setNeibElemVert(pSFace->getID(), icnode);
        };
    };
}
void CContactMesh::setupEdgeConNode(CContactMesh *pProgConMesh)
{
    PairConNode pairConNode;
    CSkinFace *pFace, *pNeibFace;
    CContactNode* pEdgeConNode;
    vdouble vCoord; vCoord.resize(3);
    uint numOfFace, numOfEdge;
    uint level; 
    uint numOfScalar,numOfDisp;
    uint face_rank,rank;
    uint iface,iedge;
    uint countID=mvConNode.size();
    uint maslave;
    for(maslave=0; maslave< 2; maslave++){
        if(maslave==0) numOfFace= mvFace.size();     
        if(maslave==1) numOfFace= mvSlaveFace.size();
        for(iface=0; iface< numOfFace; iface++){
            if(maslave==0) pFace= mvFace[iface];     
            if(maslave==1) pFace= mvSlaveFace[iface];
            numOfEdge= pFace->getNumOfEdge();
            face_rank= pFace->getRank();
            for(iedge=0; iedge< numOfEdge; iedge++){
                if(!pFace->isEdgeNodeMarking(iedge)){
                    pairConNode= pFace->getEdgePairNode(iedge);
                    uint numOfAggFaceA= pairConNode.first->getNumOfAggElem();
                    uint numOfAggFaceB= pairConNode.second->getNumOfAggElem();
                    uint faceIDa,faceIDb;
                    uint iagg,jagg;
                    bool bfind(false);
                    for(iagg=0; iagg< numOfAggFaceA; iagg++){
                        faceIDa= pairConNode.first->getAggElemID(iagg);
                        for(jagg=0; jagg< numOfAggFaceB; jagg++){
                            faceIDb= pairConNode.second->getAggElemID(jagg);
                            if(faceIDa==faceIDb && faceIDa!=pFace->getID()){
                                uint index;
                                if(maslave==0){
                                    index= mmMasterFaceID2Index[faceIDa];
                                    pNeibFace= mvFace[index];
                                }
                                if(maslave==1){
                                    index= mmSlaveFaceID2Index[faceIDa];
                                    pNeibFace= mvSlaveFace[index];
                                }
                                pEdgeConNode= new CContactNode;
                                vCoord[0]= pairConNode.first->getX(); vCoord[0]+= pairConNode.second->getX(); vCoord[0] *= 0.5;
                                vCoord[1]= pairConNode.first->getY(); vCoord[1]+= pairConNode.second->getY(); vCoord[1] *= 0.5;
                                vCoord[2]= pairConNode.first->getZ(); vCoord[2]+= pairConNode.second->getZ(); vCoord[2] *= 0.5;
                                pEdgeConNode->setCoord(vCoord);
								level= mLevel+1;
                                pEdgeConNode->setLevel(level);
								pEdgeConNode->pushLevelMarking();
                                pEdgeConNode->setID(countID);
                                countID++;
                                if(pairConNode.first->getRank()==myRank){
                                    uint meshID= pairConNode.first->getMeshID();
                                    pEdgeConNode->setMeshID(meshID);
                                }
                                if((pairConNode.first->getRank()!=face_rank) && (pairConNode.second->getRank()!=face_rank)){
                                    rank= pairConNode.first->getRank();
                                }else{
                                    rank= face_rank;
                                }
                                pEdgeConNode->setRank(rank);
                                numOfScalar= pairConNode.first->getNumOfScalar();
                                numOfDisp= pairConNode.first->getNumOfDisp();
                                pEdgeConNode->resizeDisp(numOfDisp); pEdgeConNode->initDisp();
                                pEdgeConNode->resizeScalar(numOfScalar); pEdgeConNode->initScalar();
                                pProgConMesh->addConNode(pEdgeConNode, countID);                     
                                if(maslave==0) pProgConMesh->addMasterConNode(pEdgeConNode, countID);
                                if(maslave==1) pProgConMesh->addSlaveConNode(pEdgeConNode, countID); 
                                pFace->setEdgeFace(pNeibFace, iedge);
                                pFace->setEdgeConNode(pEdgeConNode, iedge);
                                pFace->markingEdgeNode(iedge);
                                pNeibFace->setEdgeFace(pFace, pairConNode);
                                pNeibFace->setEdgeConNode(pEdgeConNode, pairConNode);
                                pNeibFace->markingEdgeNode(pairConNode);
                                bfind=true;
                                break;
                            }
                        };
                        if(bfind) break;
                    };
                    if(!bfind){
                        pEdgeConNode= new CContactNode;
                        vCoord[0]= pairConNode.first->getX(); vCoord[0]+= pairConNode.second->getX(); vCoord[0] *= 0.5;
                        vCoord[1]= pairConNode.first->getY(); vCoord[1]+= pairConNode.second->getY(); vCoord[1] *= 0.5;
                        vCoord[2]= pairConNode.first->getZ(); vCoord[2]+= pairConNode.second->getZ(); vCoord[2] *= 0.5;
                        pEdgeConNode->setCoord(vCoord);
						level= mLevel+1;
                        pEdgeConNode->setLevel(level);
						pEdgeConNode->pushLevelMarking();
                        pEdgeConNode->setID(countID);
                        countID++;
                        if(pairConNode.first->getRank()==myRank){
                            uint meshID= pairConNode.first->getMeshID();
                            pEdgeConNode->setMeshID(meshID);
                        }
                        if((pairConNode.first->getRank()!=face_rank) && (pairConNode.second->getRank()!=face_rank)){
                            rank= pairConNode.first->getRank();
                        }else{
                            rank= face_rank;
                        }
                        pEdgeConNode->setRank(rank);
                        numOfScalar= pairConNode.first->getNumOfScalar();
                        numOfDisp= pairConNode.first->getNumOfDisp();
                        pEdgeConNode->resizeDisp(numOfDisp); pEdgeConNode->initDisp();
                        pEdgeConNode->resizeScalar(numOfScalar); pEdgeConNode->initScalar();
                        pProgConMesh->addConNode(pEdgeConNode, countID);                     
                        if(maslave==0) pProgConMesh->addMasterConNode(pEdgeConNode, countID);
                        if(maslave==1) pProgConMesh->addSlaveConNode(pEdgeConNode, countID); 
                        pFace->setEdgeConNode(pEdgeConNode, iedge);
                        pFace->markingEdgeNode(iedge);
                    }
                }
            };
        };
    };
}
void CContactMesh::setupFaceConNode(CContactMesh *pProgConMesh)
{
    CSkinFace *pFace;
    CContactNode *pConNode, *pConFaceNode;
    uint numOfFace,numOfNode;
    uint numOfScalar,numOfDisp;
    uint iface,inode;
    uint countID= pProgConMesh->getNumOfConNode();
    uint face_rank;
    vdouble vCoord; vCoord.resize(3);
    uint maslave;
    for(maslave=0; maslave< 2; maslave++){
        if(maslave==0) numOfFace= mvFace.size();     
        if(maslave==1) numOfFace= mvSlaveFace.size();
        for(iface=0; iface< numOfFace; iface++){
            if(maslave==0) pFace= mvFace[iface];     
            if(maslave==1) pFace= mvSlaveFace[iface];
            face_rank= pFace->getRank();
            numOfNode= pFace->getNumOfNode();
            vCoord[0]=0.0; vCoord[1]=0.0; vCoord[2]=0.0;
            for(inode=0; inode< numOfNode; inode++){
                pConNode= pFace->getNode(inode);
                vCoord[0]+= pConNode->getX(); vCoord[1]+= pConNode->getY(); vCoord[2]+= pConNode->getZ();
            };
            vCoord[0] /= (double)numOfNode;  vCoord[1] /= (double)numOfNode;  vCoord[2] /= (double)numOfNode;
            pConFaceNode = new CContactNode;
            pConFaceNode->setCoord(vCoord);
            pConFaceNode->setID(countID);
            countID++;
            pConFaceNode->setRank(face_rank);
            pConFaceNode->setLevel(mLevel+1);
			pConFaceNode->pushLevelMarking();
            if(pConNode->getRank()==myRank){
                uint meshID= pConNode->getMeshID();
                pConFaceNode->setMeshID(meshID);
            }
            numOfScalar= pConNode->getNumOfScalar();
            numOfDisp= pConNode->getNumOfDisp();    
            pConFaceNode->resizeDisp(numOfDisp); pConFaceNode->initDisp();
            pConFaceNode->resizeScalar(numOfScalar); pConFaceNode->initScalar();
            pProgConMesh->addConNode(pConFaceNode,countID);                      
            if(maslave==0) pProgConMesh->addMasterConNode(pConFaceNode, countID);
            if(maslave==1) pProgConMesh->addSlaveConNode(pConFaceNode, countID); 
            pFace->setFaceConNode(pConFaceNode);
        };
    };
}
void CContactMesh::addMasterFace(CSkinFace* pFace)
{
    mvFace.push_back(pFace);
    mmMasterFaceID2Index[pFace->getID()]= mvFace.size()-1;
}
void CContactMesh::addMasterFace(vector<CSkinFace*>& vface)
{
    CSkinFace* pFace;
    uint i, numOfFace(vface.size());
    for(i=0; i< numOfFace; i++){
        pFace= vface[i];
        mvFace.push_back(pFace);
        mmMasterFaceID2Index[pFace->getID()]= mvFace.size()-1;
    };
}
void CContactMesh::addSlaveFace(CSkinFace* pSlaveFace)
{
    mvSlaveFace.push_back(pSlaveFace);
    mmSlaveFaceID2Index[pSlaveFace->getID()]= mvSlaveFace.size()-1;
}
void CContactMesh::addSlaveFace(vector<CSkinFace*>& vface)
{
    CSkinFace* pFace;
    uint i, numOfFace(vface.size());
    for(i=0; i< numOfFace; i++){
        pFace= vface[i];
        mvSlaveFace.push_back(pFace);
        mmSlaveFaceID2Index[pFace->getID()]= mvSlaveFace.size()-1;
    };
}
void CContactMesh::setupSPointOnMFace()
{
    uint maxLayer= mvKnot.size()-1;
    uint numOfMaster= mvFace.size();
    CSkinFace *pMFace;
    uint imface;
    for(imface=0; imface< numOfMaster; imface++){
        pMFace= mvFace[imface];
        CContactNode *pVertNode;
        uint numOfVert= pMFace->getNumOfNode();
        uint sum_knotID, knotID;
        uint finalLayer(0);
        uint iLayer, ivert;
        for(iLayer=maxLayer; iLayer >= 0; iLayer--){
            sum_knotID =0;
            for(ivert=0; ivert < numOfVert; ivert++){
                pVertNode= pMFace->getNode(ivert);
                if(ivert==0) knotID= pVertNode->getOctreeID(iLayer);
                sum_knotID += pVertNode->getOctreeID(iLayer);
            };
            if(sum_knotID == knotID*numOfVert){
                finalLayer= iLayer;
                break;
            }
        };
        COctreeKnot *pKnot= mvKnot[finalLayer][knotID];
        uint numOfSlaveNode= pKnot->getNumOfSlaveNode();
        pMFace->CalcNzNormalVector();  
        mBoundingBox.sizingOBB(pMFace);
        CContactNode* pConNode;
        uint isnode;
        for(isnode=0; isnode< numOfSlaveNode; isnode++){
            pConNode= pKnot->getSlaveNode(isnode);
            if(mBoundingBox.judgeOBB(pConNode)){
				pMFace->addSlaveNode(pConNode);
			}
        };
    };
/*
    uint numOfSlaveNode= mvSlaveConNode.size();
    CContactNode* pConNode;
    uint isnode;
    uint numOfMaster= mvFace.size();
    CSkinFace *pMFace;
    uint imface;
    for(imface=0; imface< numOfMaster; imface++){
        pMFace= mvFace[imface];
        pMFace->CalcNzNormalVector();
        mBoundingBox.sizingOBB(pMFace);
        for(isnode=0; isnode< numOfSlaveNode; isnode++){
            pConNode= mvSlaveConNode[isnode];
            if(mBoundingBox.judgeOBB(pConNode)){
                pMFace->addSlaveNode(pConNode);
                if(pConNode->getID()==22)
                    cout << "conID==22: pMFace ID= " << pMFace->getID() << endl;
            }
        };
    };
*/
}
void CContactMesh::setupMPC_Coef()
{
    uint numOfMaster= mvFace.size();
    CSkinFace *pMFace;
    uint imface;
    uint numOfSlave;
    uint islave;
    for(imface=0; imface < numOfMaster; imface++){
        pMFace= mvFace[imface];
        numOfSlave= pMFace->getNumOfSlaveNode();
        for(islave=0; islave < numOfSlave; islave++){
            pMFace->CalcSlave(islave, MPCValueType::Displacement);
        };
    };
}
CSkinFace* CContactMesh::getMasterFace_ID(const uint& id)
{
    uint index= mmMasterFaceID2Index[id];
    return mvFace[index];
}
void CContactMesh::generateOctree(const uint& maxLayer)
{
    double minX,maxX;
    double minY,maxY;
    double minZ,maxZ;
    CContactNode *pConNode;
    uint numOfConNode= mvConNode.size();
    uint inode;
    for(inode=0; inode< numOfConNode; inode++){
        pConNode= mvConNode[inode];
        pConNode->resizeOctreeID(maxLayer+1);
        if(inode==0){
            minX= pConNode->getX(); maxX= pConNode->getX();
            minY= pConNode->getY(); maxY= pConNode->getY();
            minZ= pConNode->getZ(); maxZ= pConNode->getZ();
        }else{
            if(minX > pConNode->getX()) minX= pConNode->getX();
            if(minY > pConNode->getY()) minY= pConNode->getY();
            if(minZ > pConNode->getZ()) minZ= pConNode->getZ();
            if(maxX < pConNode->getX()) maxX= pConNode->getX();
            if(maxY < pConNode->getY()) maxY= pConNode->getY();
            if(maxZ < pConNode->getZ()) maxZ= pConNode->getZ();
        }
    };
    double difX= minX - maxX;
    double difY= minY - maxY;
    double difZ= minZ - maxZ;
    double Alpha(0.1);
    if(abs(difX) <= 1.0e-4) maxX += Alpha;
    if(abs(difY) <= 1.0e-4) maxY += Alpha;
    if(abs(difZ) <= 1.0e-4) maxZ += Alpha;
    moOctreeKnot.setX(minX-Alpha, maxX+Alpha);
    moOctreeKnot.setY(minY-Alpha, maxY+Alpha);
    moOctreeKnot.setZ(minZ-Alpha, maxZ+Alpha);
    moOctreeKnot.setLayerID(0);
    moOctreeKnot.setID(0);
    uint numOfMasterNode, numOfSlaveNode;
    numOfMasterNode= mvMasterConNode.size();
    numOfSlaveNode= mvSlaveConNode.size();
    moOctreeKnot.reserveMasterNode(numOfMasterNode);
    moOctreeKnot.reserveSlaveNode(numOfSlaveNode);
    for(inode=0; inode < numOfMasterNode; inode++){
        pConNode= mvMasterConNode[inode];
        moOctreeKnot.addMasterNode(pConNode);
    };
    for(inode=0; inode < numOfSlaveNode; inode++){
        pConNode= mvSlaveConNode[inode];
        moOctreeKnot.addSlaveNode(pConNode);
    };
    moOctreeKnot.setItemProp();
    COctreeKnot *pPrevKnot,*pNextKnot;
    vector<COctreeKnot*> nextKnot;
    vector<COctreeKnot*> prevKnot;
    mvKnot.resize(maxLayer+1);         
    mvKnot[0].push_back(&moOctreeKnot);
    uint lastPos;
    uint ilayer, prev_pos, child_pos;
    for(ilayer=0; ilayer < maxLayer; ilayer++){
        uint knotID(0); 
        prevKnot= mvKnot[ilayer];      
        lastPos= mvKnot[ilayer].size();
        nextKnot.clear();
        for(prev_pos=0; prev_pos < lastPos; prev_pos++){
            pPrevKnot= prevKnot[prev_pos];
            pPrevKnot->createChildKnot();
            pPrevKnot->distItem();
            for(child_pos=0; child_pos < 8; child_pos++){
                pNextKnot= pPrevKnot->getChildKnot(child_pos);
                pNextKnot->setLayerID(ilayer+1);
                pNextKnot->setID(knotID);
                knotID++;
                pNextKnot->setItemProp();
                nextKnot.push_back(pNextKnot);
            };
        };
        mvKnot[ilayer+1]= nextKnot;
    };
}
