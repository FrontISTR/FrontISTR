/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ContactMesh.cpp
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
    CContactNode *pConNode;
    uiint i;
    for(i=0; i < mvConNode.size() ; i++) {
        pConNode = mvConNode[i];
        if(mLevel==pConNode->getLevel()) delete pConNode;
        mvConNode.erase(mvConNode.begin()+i);
    };
    for_each(mvFace.begin(), mvFace.end(), DeleteObject());
    for_each(mvSlaveFace.begin(), mvSlaveFace.end(), DeleteObject());

    vector<COctreeKnot*> vKnot;
    uiint numOfLayer= mvKnot.size();
    uiint ilayer;
    for(ilayer=1; ilayer < numOfLayer; ilayer++) {
        vKnot= mvKnot[ilayer];
        for_each(vKnot.begin(), vKnot.end(), DeleteObject());
    };
}
void CContactMesh::addConNode(CContactNode* pConNode, const uiint& id)
{
    mvConNode.push_back(pConNode);
    mmConNodeID2Index[id]= mvConNode.size()-1;
}
void CContactMesh::addMasterConNode(CContactNode* pConNode, const uiint& id)
{
    mvMasterConNode.push_back(pConNode);
    mmMasterConNodeID2Index[id]= mvMasterConNode.size()-1;
}
void CContactMesh::addSlaveConNode(CContactNode* pConNode, const uiint& id)
{
    mvSlaveConNode.push_back(pConNode);
    mmSlaveConNodeID2Index[id]= mvSlaveConNode.size()-1;
}
void CContactMesh::addMasterMeshID(const uiint& id)
{
    mvMasterMeshID.push_back(id);
    uiint index= mvMasterMeshID.size()-1;
    mmMasterMeshID2Index[id]= index;
}
void CContactMesh::addSlaveMeshID(const uiint& id)
{
    mvSlaveMeshID.push_back(id);
    uiint index= mvSlaveMeshID.size()-1;
    mmSlaveMeshID2Index[id]= index;
}
void CContactMesh::setupCoarseConNode(CContactMesh* pProgConMesh)
{
    uiint numOfNode;
    numOfNode= mvConNode.size();
    CContactNode *pConNode;
    uiint inode, id;
    for(inode=0; inode< numOfNode; inode++) {
        pConNode= mvConNode[inode];
        pConNode->pushLevelMarking();
        id= pConNode->getID();
        pProgConMesh->addConNode(pConNode, id);
    };
    numOfNode= mvMasterConNode.size();
    for(inode=0; inode< numOfNode; inode++) {
        pConNode= mvMasterConNode[inode];
        id= pConNode->getID();
        pProgConMesh->addMasterConNode(pConNode, id);
    };
    numOfNode= mvSlaveConNode.size();
    for(inode=0; inode< numOfNode; inode++) {
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
    uiint numOfMFace(mvFace.size()), numOfSFace(mvSlaveFace.size()), numOfConNode;
    uiint iface,icnode;
    numOfConNode= mvMasterConNode.size();
    for(icnode=0; icnode< numOfConNode; icnode++) {
        pConNode= mvMasterConNode[icnode];
        pConNode->clearAggElemID();
        pConNode->clearNeibElemVert();
    };
    numOfConNode= mvSlaveConNode.size();
    for(icnode=0; icnode< numOfConNode; icnode++) {
        pConNode= mvSlaveConNode[icnode];
        pConNode->clearAggElemID();
        pConNode->clearNeibElemVert();
    };

    uiint nNumOfVert;
    // マスター面のConNodeのAggregateFace
    for(iface=0; iface< numOfMFace; iface++) {
        pMFace= mvFace[iface];
        nNumOfVert= pMFace->getNumOfVert();
        for(icnode=0; icnode< nNumOfVert; icnode++) {
            pConNode= pMFace->getNode(icnode);
            pConNode->setAggElemID(pMFace->getID());
            pConNode->setNeibElemVert(pMFace->getID(), icnode);
        };
    };
    // スレーブ面のConNodeのAggregateFace
    for(iface=0; iface< numOfSFace; iface++) {
        pSFace= mvSlaveFace[iface];
        nNumOfVert= pSFace->getNumOfVert();
        for(icnode=0; icnode< nNumOfVert; icnode++) {
            pConNode= pSFace->getNode(icnode);
            pConNode->setAggElemID(pSFace->getID());
            pConNode->setNeibElemVert(pSFace->getID(), icnode);
        };
    };

    //--
    // Overlap ConNode : Master
    //--
    numOfConNode= mvMasterConNode.size();
    for(icnode=0; icnode< numOfConNode; icnode++) {
        pConNode= mvMasterConNode[icnode];

        uiint nNumOfAggFace= pConNode->getNumOfAggElem();
        for(uiint iagg=0; iagg < nNumOfAggFace; iagg++) {
            uiint nFaceID= pConNode->getAggElemID(iagg);
            uiint nFaceIX= mmMasterFaceID2Index[nFaceID];
            pMFace= mvFace[nFaceIX];

            // Overlap 点 : rank値のセットはパーティショナーで利用すべき値
            if(pConNode->getRank()!=pMFace->getRank()) {
                pConNode->markingOverlap();//---周囲の面と点のRank違い=>複数Rankに属している点
                pConNode->addOverlapRank(pMFace->getRank());//---自身以外のRankを,まずはセット
            }
            if(pConNode->isOverlap())
                pConNode->addOverlapRank(pConNode->getRank());//--自身のRankも追加
        };//-- agg_face loop end

    };//- cnode loop end

    //--
    // Overlap ConNode : Slave
    //--
    numOfConNode= mvSlaveConNode.size();
    for(icnode=0; icnode< numOfConNode; icnode++) {
        pConNode= mvSlaveConNode[icnode];

        uiint nNumOfAggFace= pConNode->getNumOfAggElem();
        for(uiint iagg=0; iagg < nNumOfAggFace; iagg++) {
            uiint nFaceID= pConNode->getAggElemID(iagg);
            uiint nFaceIX= mmSlaveFaceID2Index[nFaceID];
            pSFace= mvSlaveFace[nFaceIX];

            // Overlap 点 : rank値のセットはパーティショナーで利用すべき値
            if(pConNode->getRank()!=pSFace->getRank()) {
                pConNode->markingOverlap();//---周囲の面と点のRank違い=>複数Rankに属している点
                pConNode->addOverlapRank(pSFace->getRank());//---自身以外のRankを,まずはセット
            }
        };//-- agg_face loop end

        if(pConNode->isOverlap())
            pConNode->addOverlapRank(pConNode->getRank());//--自身のRankも追加

    };//- cnode loop end

}
void CContactMesh::setupEdgeConNode(CContactMesh *pProgConMesh, const uiint& iLevel)
{
    PairConNode pairConNode;
    CSkinFace *pFace, *pNeibFace;
    CContactNode* pEdgeConNode;
    vdouble vCoord;
    vCoord.resize(3);
    uiint numOfFace, numOfEdge;
    uiint level;
    uiint numOfScalar,numOfDisp;
    uiint face_rank,rank;
    uiint iface,iedge;
    uiint countID=mvConNode.size();
    uiint maslave;
    for(maslave=0; maslave< 2; maslave++) {
        if(maslave==0) numOfFace= mvFace.size();
        if(maslave==1) numOfFace= mvSlaveFace.size();
        for(iface=0; iface< numOfFace; iface++) {
            if(maslave==0) pFace= mvFace[iface];
            if(maslave==1) pFace= mvSlaveFace[iface];
            bool bExecFlag(false);
            if(pFace->getOrder()==ElementOrder::Second && iLevel > 0 ) bExecFlag=true;
            if(pFace->getOrder()==ElementOrder::First )   bExecFlag=true;
            if(bExecFlag) {
                numOfEdge= pFace->getNumOfEdge();
                face_rank= pFace->getRank();
                for(iedge=0; iedge< numOfEdge; iedge++) {
                    if(!pFace->isEdgeNodeMarking(iedge)) {
                        pairConNode= pFace->getEdgePairNode(iedge);
                        uiint numOfAggFaceA= pairConNode.first->getNumOfAggElem();
                        uiint numOfAggFaceB= pairConNode.second->getNumOfAggElem();
                        uiint faceIDa,faceIDb;
                        uiint iagg,jagg;
                        bool bfind(false);
                        for(iagg=0; iagg< numOfAggFaceA; iagg++) {
                            faceIDa= pairConNode.first->getAggElemID(iagg);
                            for(jagg=0; jagg< numOfAggFaceB; jagg++) {
                                faceIDb= pairConNode.second->getAggElemID(jagg);
                                //--
                                // 隣接Face
                                //--
                                if(faceIDa==faceIDb && faceIDa!=pFace->getID()) {
                                    uiint index;
                                    if(maslave==0) {
                                        index= mmMasterFaceID2Index[faceIDa];
                                        pNeibFace= mvFace[index];
                                    }
                                    if(maslave==1) {
                                        index= mmSlaveFaceID2Index[faceIDa];
                                        pNeibFace= mvSlaveFace[index];
                                    }
                                    pEdgeConNode= new CContactNode;
                                    vCoord[0]= pairConNode.first->getX();
                                    vCoord[0]+= pairConNode.second->getX();
                                    vCoord[0] *= 0.5;
                                    vCoord[1]= pairConNode.first->getY();
                                    vCoord[1]+= pairConNode.second->getY();
                                    vCoord[1] *= 0.5;
                                    vCoord[2]= pairConNode.first->getZ();
                                    vCoord[2]+= pairConNode.second->getZ();
                                    vCoord[2] *= 0.5;
                                    pEdgeConNode->setCoord(vCoord);
                                    if(pFace->getOrder()==ElementOrder::First) {
                                        level= mLevel+1;
                                        pEdgeConNode->setLevel(level);
                                        pEdgeConNode->pushLevelMarking();
                                    }
                                    if(pFace->getOrder()==ElementOrder::Second) {
                                        pEdgeConNode->setLevel(mLevel);
                                        pEdgeConNode->pushLevelMarking();
                                    }
                                    pEdgeConNode->setID(countID);
                                    countID++;
                                    if(pairConNode.first->getRank()==myRank) {
                                        uiint meshID= pairConNode.first->getMeshID();
                                        pEdgeConNode->setMeshID(meshID);
                                    }

                                    // 節点ランクの設定 '12.03.19
                                    if( pairConNode.first->getRank() == pairConNode.second->getRank() ) {
                                        rank= pairConNode.first->getRank();//---両端のランクが等しい：点ランクに合わせる.
                                    } else {
                                        if( pFace->getRank() >  pNeibFace->getRank() ) rank= pNeibFace->getRank();
                                        if( pFace->getRank() <= pNeibFace->getRank() ) rank= pFace->getRank();
                                        ;//----------------------両端のランクが異なる：隣接する面ランクの小ランクを選択.
                                    }
                                    pEdgeConNode->setRank(rank);//------------------- rank セット(1/2)

                                    // Selfマーキング
                                    if( pairConNode.first->isSelfNode() && pairConNode.second->isSelfNode() ) {
                                        pEdgeConNode->markingSelfNode();
                                        pEdgeConNode->markingSelfMesh();
                                    }

                                    numOfScalar= pairConNode.first->getNumOfScalar();
                                    numOfDisp= pairConNode.first->getNumOfDisp();
                                    pEdgeConNode->resizeDisp(numOfDisp);
                                    pEdgeConNode->initDisp();
                                    pEdgeConNode->resizeScalar(numOfScalar);
                                    pEdgeConNode->initScalar();
                                    if(pProgConMesh) {
                                        pProgConMesh->addConNode(pEdgeConNode, countID);
                                        if(maslave==0) pProgConMesh->addMasterConNode(pEdgeConNode, countID);
                                        if(maslave==1) pProgConMesh->addSlaveConNode(pEdgeConNode, countID);
                                    }
                                    pFace->setEdgeFace(pNeibFace, iedge);
                                    pFace->setEdgeConNode(pEdgeConNode, iedge);
                                    pFace->markingEdgeNode(iedge);
                                    if(pFace->getOrder()==ElementOrder::Second) {
                                        this->addConNode(pEdgeConNode, countID);
                                        if(maslave==0) this->addMasterConNode(pEdgeConNode, countID);
                                        if(maslave==1) this->addSlaveConNode(pEdgeConNode,countID);
                                    }
                                    pNeibFace->setEdgeFace(pFace, pairConNode);
                                    pNeibFace->setEdgeConNode(pEdgeConNode, pairConNode);
                                    pNeibFace->markingEdgeNode(pairConNode);
                                    bfind=true;//-- 発見フラグ
                                    break;
                                }
                            };
                            if(bfind) break;
                        };
                        if(!bfind) {
                            //--
                            // 隣接する面が存在しない.
                            //--
                            pEdgeConNode= new CContactNode;
                            vCoord[0]= pairConNode.first->getX();
                            vCoord[0]+= pairConNode.second->getX();
                            vCoord[0] *= 0.5;
                            vCoord[1]= pairConNode.first->getY();
                            vCoord[1]+= pairConNode.second->getY();
                            vCoord[1] *= 0.5;
                            vCoord[2]= pairConNode.first->getZ();
                            vCoord[2]+= pairConNode.second->getZ();
                            vCoord[2] *= 0.5;
                            pEdgeConNode->setCoord(vCoord);
                            if(pFace->getOrder()==ElementOrder::First) {
                                level= mLevel+1;
                                pEdgeConNode->setLevel(level);
                                pEdgeConNode->pushLevelMarking();
                            }
                            if(pFace->getOrder()==ElementOrder::Second) {
                                pEdgeConNode->setLevel(mLevel);
                                pEdgeConNode->pushLevelMarking();
                            }
                            pEdgeConNode->setID(countID);
                            countID++;
                            if(pairConNode.first->getRank()==myRank) {
                                uiint meshID= pairConNode.first->getMeshID();
                                pEdgeConNode->setMeshID(meshID);
                            }

                            // 節点ランクの設定 '12.03.19
                            if( pairConNode.first->getRank() == pairConNode.second->getRank() ) {
                                rank= pairConNode.first->getRank();//---両端のランクが等しい：点ランクに合わせる.
                            } else {
                                rank= pFace->getRank();
                                ;//----------------------両端のランクが異なる：面ランクを選択.(隣接する面は無い)
                            }
                            pEdgeConNode->setRank(rank);//-------------------- rank セット(2/2)

                            // Selfマーキング
                            if( pairConNode.first->isSelfNode() && pairConNode.second->isSelfNode() ) {
                                pEdgeConNode->markingSelfNode();
                                pEdgeConNode->markingSelfMesh();
                            }

                            numOfScalar= pairConNode.first->getNumOfScalar();
                            numOfDisp= pairConNode.first->getNumOfDisp();
                            pEdgeConNode->resizeDisp(numOfDisp);
                            pEdgeConNode->initDisp();
                            pEdgeConNode->resizeScalar(numOfScalar);
                            pEdgeConNode->initScalar();
                            if(pProgConMesh) {
                                pProgConMesh->addConNode(pEdgeConNode, countID);
                                if(maslave==0) pProgConMesh->addMasterConNode(pEdgeConNode, countID);
                                if(maslave==1) pProgConMesh->addSlaveConNode(pEdgeConNode, countID);
                            }
                            pFace->setEdgeConNode(pEdgeConNode, iedge);
                            pFace->markingEdgeNode(iedge);
                            if(pFace->getOrder()==ElementOrder::Second) {
                                this->addConNode(pEdgeConNode, countID);
                                if(maslave==0) this->addMasterConNode(pEdgeConNode, countID);
                                if(maslave==1) this->addSlaveConNode(pEdgeConNode,countID);
                            }
                        }
                    }
                };
                pFace->replaceEdgeNode();
            }
        };
    };
}
void CContactMesh::setupFaceConNode(CContactMesh *pProgConMesh)
{
    CSkinFace *pFace;
    CContactNode *pConNode, *pConFaceNode;
    uiint numOfFace,numOfNode;
    uiint numOfScalar,numOfDisp;
    uiint iface,inode;
    uiint countID= pProgConMesh->getNumOfConNode();
    uiint face_rank;
    vdouble vCoord;
    vCoord.resize(3);
    uiint maslave;
    for(maslave=0; maslave< 2; maslave++) {
        if(maslave==0) numOfFace= mvFace.size();
        if(maslave==1) numOfFace= mvSlaveFace.size();
        for(iface=0; iface< numOfFace; iface++) {
            if(maslave==0) pFace= mvFace[iface];
            if(maslave==1) pFace= mvSlaveFace[iface];

            face_rank= pFace->getRank();
            numOfNode= pFace->getNumOfNode();
            vCoord[0]=0.0;
            vCoord[1]=0.0;
            vCoord[2]=0.0;
            for(inode=0; inode< numOfNode; inode++) {
                pConNode= pFace->getNode(inode);
                vCoord[0]+= pConNode->getX();
                vCoord[1]+= pConNode->getY();
                vCoord[2]+= pConNode->getZ();
            };
            vCoord[0] /= (double)numOfNode;
            vCoord[1] /= (double)numOfNode;
            vCoord[2] /= (double)numOfNode;
            pConFaceNode = new CContactNode;
            pConFaceNode->setCoord(vCoord);
            pConFaceNode->setID(countID);
            countID++;
            pConFaceNode->setRank(face_rank);//------------- rank セット
            pConFaceNode->setLevel(mLevel+1);
            pConFaceNode->pushLevelMarking();

            if(pConNode->getRank()==myRank) {
                uiint meshID= pConNode->getMeshID();
                pConFaceNode->setMeshID(meshID);
            }

            // Selfマーキング
            if( pFace->isSelf() ) {
                pConFaceNode->markingSelfNode();
                pConFaceNode->markingSelfMesh();
            }

            numOfScalar= pConNode->getNumOfScalar();
            numOfDisp= pConNode->getNumOfDisp();
            pConFaceNode->resizeDisp(numOfDisp);
            pConFaceNode->initDisp();
            pConFaceNode->resizeScalar(numOfScalar);
            pConFaceNode->initScalar();
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
    uiint i, numOfFace(vface.size());
    for(i=0; i< numOfFace; i++) {
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
    uiint i, numOfFace(vface.size());
    for(i=0; i< numOfFace; i++) {
        pFace= vface[i];
        mvSlaveFace.push_back(pFace);
        mmSlaveFaceID2Index[pFace->getID()]= mvSlaveFace.size()-1;
    };
}
void CContactMesh::setupSPointOnMFace()
{
    uiint maxLayer= mvKnot.size()-1;
    uiint numOfMaster= mvFace.size();
    CSkinFace *pMFace;
    uiint imface;
    for(imface=0; imface< numOfMaster; imface++) {
        pMFace= mvFace[imface];
        CContactNode *pVertNode;
        uiint numOfVert= pMFace->getNumOfNode();
        uiint sum_knotID, knotID;
        uiint finalLayer(0);
        uiint iLayer, ivert;
        for(iLayer=maxLayer; iLayer >= 0; iLayer--) {
            sum_knotID =0;
            for(ivert=0; ivert < numOfVert; ivert++) {
                pVertNode= pMFace->getNode(ivert);
                if(ivert==0) knotID= pVertNode->getOctreeID(iLayer);
                sum_knotID += pVertNode->getOctreeID(iLayer);
            };
            if(sum_knotID == knotID*numOfVert) {
                finalLayer= iLayer;
                break;
            }
        };
        COctreeKnot *pKnot= mvKnot[finalLayer][knotID];
        uiint numOfSlaveNode= pKnot->getNumOfSlaveNode();
        pMFace->CalcNzNormalVector();
        mBoundingBox.sizingOBB(pMFace);
        CContactNode* pConNode;
        uiint isnode;
        for(isnode=0; isnode< numOfSlaveNode; isnode++) {
            pConNode= pKnot->getSlaveNode(isnode);
            if(mBoundingBox.judgeOBB(pConNode)) {
                pMFace->addSlaveNode(pConNode);
            }
        };
    };
}
void CContactMesh::setupMPC_Coef()
{
    uiint numOfMaster= mvFace.size();
    CSkinFace *pMFace;
    uiint imface;
    uiint numOfSlave;
    uiint islave;
    for(imface=0; imface < numOfMaster; imface++) {
        pMFace= mvFace[imface];
        numOfSlave= pMFace->getNumOfSlaveNode();
        for(islave=0; islave < numOfSlave; islave++) {
            pMFace->CalcSlave(islave, MPCValueType::Displacement);
        };
    };
}
CSkinFace* CContactMesh::getMasterFace_ID(const uiint& id)
{
    uiint index= mmMasterFaceID2Index[id];
    return mvFace[index];
}
void CContactMesh::generateOctree(const uiint& maxLayer)
{
    double minX,maxX;
    double minY,maxY;
    double minZ,maxZ;
    CContactNode *pConNode;
    uiint numOfConNode= mvConNode.size();
    uiint inode;
    for(inode=0; inode< numOfConNode; inode++) {
        pConNode= mvConNode[inode];
        pConNode->resizeOctreeID(maxLayer+1);
        if(inode==0) {
            minX= pConNode->getX();
            maxX= pConNode->getX();
            minY= pConNode->getY();
            maxY= pConNode->getY();
            minZ= pConNode->getZ();
            maxZ= pConNode->getZ();
        } else {
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
    uiint numOfMasterNode, numOfSlaveNode;
    numOfMasterNode= mvMasterConNode.size();
    numOfSlaveNode= mvSlaveConNode.size();
    moOctreeKnot.reserveMasterNode(numOfMasterNode);
    moOctreeKnot.reserveSlaveNode(numOfSlaveNode);
    for(inode=0; inode < numOfMasterNode; inode++) {
        pConNode= mvMasterConNode[inode];
        moOctreeKnot.addMasterNode(pConNode);
    };
    for(inode=0; inode < numOfSlaveNode; inode++) {
        pConNode= mvSlaveConNode[inode];
        moOctreeKnot.addSlaveNode(pConNode);
    };
    moOctreeKnot.setItemProp();

    COctreeKnot *pPrevKnot,*pNextKnot;
    vector<COctreeKnot*> nextKnot;
    vector<COctreeKnot*> prevKnot;
    mvKnot.resize(maxLayer+1);
    mvKnot[0].push_back(&moOctreeKnot);
    uiint lastPos;
    uiint ilayer, prev_pos, child_pos;
    for(ilayer=0; ilayer < maxLayer; ilayer++) {
        uiint knotID(0);
        prevKnot= mvKnot[ilayer];
        lastPos= mvKnot[ilayer].size();
        nextKnot.clear();
        for(prev_pos=0; prev_pos < lastPos; prev_pos++) {
            pPrevKnot= prevKnot[prev_pos];
            pPrevKnot->createChildKnot();
            pPrevKnot->distItem();
            for(child_pos=0; child_pos < 8; child_pos++) {
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
void CContactMesh::deleteProgData()
{
    uiint nNumOfFace;
    uiint iFace;
    nNumOfFace = mvFace.size();
    for(iFace=0; iFace < nNumOfFace; iFace++) mvFace[iFace]->deleteProgData();
    nNumOfFace = mvSlaveFace.size();
    for(iFace=0; iFace < nNumOfFace; iFace++) mvSlaveFace[iFace]->deleteProgData();
}

void CContactMesh::clear()
{
    mmMasterMeshID2Index.clear();
    mmSlaveMeshID2Index.clear();

    uiint nNumOfCNode= mvConNode.size();
    for(uiint i=0; i < nNumOfCNode; i++) {
        delete mvConNode[i];
    };
    mvConNode.clear();
    mmConNodeID2Index.clear();

    uiint nNumOfMCNode= mvMasterConNode.size();
    for(uiint i=0; i < nNumOfMCNode; i++) {
        delete mvMasterConNode[i];
    }
    mvMasterConNode.clear();
    uiint nNumOfMFace= mvFace.size();
    for(uiint i=0; i < nNumOfMFace; i++) {
        delete mvFace[i];
    }
    mvFace.clear();

    mmMasterConNodeID2Index.clear();
    mmMasterFaceID2Index.clear();

    uiint nNumOfSCNode= mvSlaveConNode.size();
    for(uiint i=0; i < nNumOfSCNode; i++) {
        delete mvSlaveConNode[i];
    }
    mvSlaveConNode.clear();
    uiint nNumOfSFace= mvSlaveFace.size();
    for(uiint i=0; i < nNumOfSFace; i++) {
        delete mvSlaveFace[i];
    }
    mvSlaveFace.clear();

    mmSlaveConNodeID2Index.clear();
    mmSlaveFaceID2Index.clear();
}

//--
// 熱伝達率
//--
double& CContactMesh::getTransCoeff(const uiint& ieq)
{
    return mpFilm->getTransCoeff(ieq);
}


