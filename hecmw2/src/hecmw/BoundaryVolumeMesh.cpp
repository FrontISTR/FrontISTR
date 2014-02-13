/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryVolumeMesh.cpp
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
#include "BndVertex.h"
#include "BoundaryMesh.h"
#include "BoundaryVolume.h"
#include "BoundaryParts.h"
#include "BoundaryVolumeMesh.h"
using namespace pmw;
CBoundaryVolumeMesh::CBoundaryVolumeMesh()
{
    ;
}
CBoundaryVolumeMesh::~CBoundaryVolumeMesh()
{
    for_each(mvBVolume.begin(), mvBVolume.end(), DeleteObject());
}
void CBoundaryVolumeMesh::resizeVolume(const uiint& res_size)
{
    mvBVolume.resize(res_size);
}
void CBoundaryVolumeMesh::setBVolume(const uiint& index, CBoundaryVolume* pBVolume)
{
    uiint id;
    id= pBVolume->getID();
    mmBVolumeID2Index[id]= index;
    mvBVolume[index]= pBVolume;
}
void CBoundaryVolumeMesh::addBVolume(CBoundaryVolume* pBVolume)
{
    mvBVolume.push_back(pBVolume);
    uiint id;
    id= pBVolume->getID();
    mmBVolumeID2Index[id]= mvBVolume.size()-1;
}
void CBoundaryVolumeMesh::resizeAggVol()
{
    uiint numOfAgg= mvBNode.size();
    mvAggregateVol.resize(numOfAgg);
}
void CBoundaryVolumeMesh::setupAggVol()
{
    uiint numOfBNode= mvBNode.size();
    CBoundaryNode  *pBNode;
    uiint ibnode;
    for(ibnode=0; ibnode < numOfBNode; ibnode++) {
        pBNode= mvBNode[ibnode];
        pBNode->clearAggElemID();
        pBNode->clearNeibElemVert();
    };
    uiint numOfVol= mvBVolume.size();
    CBoundaryVolume *pBVol;
    uiint ivol;
    for(ivol=0; ivol < numOfVol; ivol++) {
        pBVol= mvBVolume[ivol];
        uiint numOfVert= pBVol->getNumOfVert();
        uiint ivert;
        for(ivert=0; ivert < numOfVert; ivert++) {
            pBNode= pBVol->getBNode(ivert);
            pBNode->setAggElemID(pBVol->getID());
        };
    };
    for(ibnode=0; ibnode < numOfBNode; ibnode++) {
        pBNode= mvBNode[ibnode];
        uiint numOfAgg= pBNode->getNumOfAggElem();
        uiint iAgg;
        for(iAgg=0; iAgg < numOfAgg; iAgg++) {
            mvAggregateVol[ibnode].push_back(pBNode->getAggElemID(iAgg));
        };
    };
}
void CBoundaryVolumeMesh::GeneEdgeBNode()
{
    uiint numOfVol= mvBVolume.size();
    uiint ivol;
    CBoundaryVolume *pBVol;
    CElement *pElem;
    CNode *pNode0,*pNode1;
    CNode *pEdgeNode;
    uiint countID= mvBNode.size();
    for(ivol=0; ivol < numOfVol; ivol++) {
        pBVol= mvBVolume[ivol];
        uiint numOfEdge= pBVol->getNumOfEdge();
        uiint iedge;
        PairBNode pairBNode;
        for(iedge=0; iedge < numOfEdge; iedge++) {
            if(!pBVol->isMarkingEdge(iedge)) {
                bool bfind(false);
                pairBNode= pBVol->getPairBNode(iedge);
                pElem = pBVol->getElement();
                pNode0 = pairBNode.first->getNode();
                pNode1 = pairBNode.second->getNode();
                uiint nElemEdge= pElem->getEdgeIndex(pNode0, pNode1);
                uiint numOfiAgg= pairBNode.first->getNumOfAggElem();
                uiint numOfjAgg= pairBNode.second->getNumOfAggElem();
                uiint iAgg,jAgg, iVolID, jVolID;
                for(iAgg=0; iAgg < numOfiAgg; iAgg++) {
                    iVolID= pairBNode.first->getAggElemID(iAgg);
                    for(jAgg=0; jAgg < numOfjAgg; jAgg++) {
                        jVolID= pairBNode.second->getAggElemID(jAgg);
                        if(iVolID==jVolID) {
                            if(iVolID != pBVol->getID()) {
                                bfind= true;
                                CBoundaryNode *pEdgeBNode= new CBoundaryNode;//--------------------
                                pEdgeBNode->setID(countID);
                                countID++;
                                mvBEdgeBNode.push_back(pEdgeBNode);
                                pBVol->setEdgeNeibVol(iedge, iVolID);
                                pBVol->markingEdge(iedge);
                                pBVol->setEdgeBNode(iedge, pEdgeBNode);
                                uiint neibIndex= mmBVolumeID2Index[iVolID];
                                CBoundaryVolume *pNeibBVol= mvBVolume[neibIndex];
                                uiint jedge= pNeibBVol->getEdgeID(pairBNode);
                                pNeibBVol->setEdgeNeibVol(jedge, pBVol->getID());
                                pNeibBVol->markingEdge(jedge);
                                pNeibBVol->setEdgeBNode(jedge, pEdgeBNode);
                                pEdgeNode = pElem->getEdgeInterNode(nElemEdge);
                                pEdgeBNode->setNode(pEdgeNode);
                                if(pBVol->getOrder()==ElementOrder::Second) {
                                    uiint renum = mvBNode.size();
                                    mvBNode.resize(renum+1);
                                    mvBNode[renum]=pEdgeBNode;
                                    mmBNodeID2Index[pEdgeBNode->getID()] = renum;
                                    mvAggregateVol.resize(renum+1);
                                    mvAggregateVol[renum].push_back(pBVol->getID());
                                    mvAggregateVol[renum].push_back(pNeibBVol->getID());
                                    pBVol->replaceEdgeBNode(iedge);
                                    pNeibBVol->replaceEdgeBNode(jedge);
                                    pEdgeBNode->setMGLevel(mMGLevel);
                                    pEdgeBNode->resizeValue(mMaxMGLevel-mMGLevel + 1);
                                } else {
                                    pEdgeBNode->setMGLevel(mMGLevel+1);
                                    pEdgeBNode->resizeValue(mMaxMGLevel-mMGLevel);
                                    mnEdgeNodeCount++;
                                }
                            }
                        }
                        if(bfind) break;
                    };
                    if(bfind) break;
                };
                if(!bfind) {
                    CBoundaryNode *pEdgeBNode= new CBoundaryNode;//---------------------
                    pEdgeBNode->setID(countID);
                    countID++;
                    mvBEdgeBNode.push_back(pEdgeBNode);
                    pBVol->markingEdge(iedge);
                    pBVol->setEdgeBNode(iedge, pEdgeBNode);
                    pEdgeNode = pElem->getEdgeInterNode(nElemEdge);
                    pEdgeBNode->setNode(pEdgeNode);
                    if(pBVol->getOrder()==ElementOrder::Second) {
                        uiint renum = mvBNode.size();
                        mvBNode.resize(renum+1);
                        mvBNode[renum]=pEdgeBNode;
                        mmBNodeID2Index[pEdgeBNode->getID()] = renum;
                        mvAggregateVol.resize(renum+1);
                        mvAggregateVol[renum].push_back(pBVol->getID());
                        pBVol->replaceEdgeBNode(iedge);
                        pEdgeBNode->setMGLevel(mMGLevel);
                        pEdgeBNode->resizeValue(mMaxMGLevel-mMGLevel + 1);
                    } else {
                        pEdgeBNode->setMGLevel(mMGLevel+1);
                        pEdgeBNode->resizeValue(mMaxMGLevel-mMGLevel);
                        mnEdgeNodeCount++;
                    }
                }
            }
        };
    };
}
void CBoundaryVolumeMesh::GeneFaceBNode()
{
    uiint numOfVol= mvBVolume.size();
    CBoundaryVolume *pBVol;
    CBoundaryVolume *pNeibBVol;
    CElement *pElem;
    vector<CNode*> vNode;
    vNode.resize(3);
    CNode *pFaceNode;
    uiint countID= mvBNode.size() + mvBEdgeBNode.size();
    uiint ivol;
    for(ivol=0; ivol < numOfVol; ivol++) {
        pBVol= mvBVolume[ivol];
        uiint numOfFace= pBVol->getNumOfFace();
        uiint iface;
        vector<CBoundaryNode*> vBNode;
        for(iface=0; iface < numOfFace; iface++) {
            if(!pBVol->isMarkingFace(iface)) {
                vBNode= pBVol->getFaceCnvNodes(iface);
                uiint iVolID, jVolID, kVolID;
                uiint iAgg, jAgg, kAgg;
                uiint numOfiAgg= vBNode[0]->getNumOfAggElem();
                uiint numOfjAgg= vBNode[1]->getNumOfAggElem();
                uiint numOfkAgg= vBNode[2]->getNumOfAggElem();
                bool bfind(false);
                for(iAgg=0; iAgg < numOfiAgg; iAgg++) {
                    iVolID= vBNode[0]->getAggElemID(iAgg);
                    for(jAgg=0; jAgg < numOfjAgg; jAgg++) {
                        jVolID= vBNode[1]->getAggElemID(jAgg);
                        for(kAgg=0; kAgg < numOfkAgg; kAgg++) {
                            kVolID= vBNode[2]->getAggElemID(kAgg);
                            if(iVolID==jVolID && jVolID==kVolID && kVolID!=pBVol->getID()) {
                                bfind= true;
                                CBoundaryNode *pFaceBNode= new CBoundaryNode;
                                pFaceBNode->setMGLevel(mMGLevel+1);
                                pFaceBNode->resizeValue(mMaxMGLevel-mMGLevel);
                                pFaceBNode->setID(countID);
                                countID++;
                                mvBFaceBNode.push_back(pFaceBNode);
                                pBVol->setFaceNeibVol(iface, iVolID);
                                pBVol->markingFace(iface);
                                pBVol->setFaceBNode(iface, pFaceBNode);
                                uiint neibIndex= mmBVolumeID2Index[iVolID];
                                pNeibBVol= mvBVolume[neibIndex];
                                uiint neibFace= pNeibBVol->getFaceID(vBNode);
                                pNeibBVol->setFaceNeibVol(neibFace, pBVol->getID());
                                pNeibBVol->markingFace(neibFace);
                                pNeibBVol->setFaceBNode(neibFace, pFaceBNode);
                                pElem = pBVol->getElement();
                                uiint ivert;
                                for(ivert=0; ivert < 3; ivert++) {
                                    vNode[ivert] = vBNode[ivert]->getNode();
                                };
                                uiint nFaceIndex = pElem->getFaceIndex(vNode[0],vNode[1],vNode[2]);
                                pFaceNode = pElem->getFaceNode(nFaceIndex);
                                pFaceBNode->setNode(pFaceNode);
                            }
                            if(bfind) break;
                        };
                        if(bfind) break;
                    };
                    if(bfind) break;
                };
                if(!bfind) {
                    CBoundaryNode *pFaceBNode= new CBoundaryNode;
                    pFaceBNode->setMGLevel(mMGLevel+1);
                    pFaceBNode->resizeValue(mMaxMGLevel-mMGLevel);
                    pFaceBNode->setID(countID);
                    countID++;
                    mvBFaceBNode.push_back(pFaceBNode);
                    pBVol->markingFace(iface);
                    pBVol->setFaceBNode(iface, pFaceBNode);
                    pElem = pBVol->getElement();
                    uiint ivert;
                    for(ivert=0; ivert < 3; ivert++) {
                        vNode[ivert] = vBNode[ivert]->getNode();
                    };
                    uiint nFaceIndex = pElem->getFaceIndex(vNode[0],vNode[1],vNode[2]);
                    pFaceNode = pElem->getFaceNode(nFaceIndex);
                    pFaceBNode->setNode(pFaceNode);
                }
            }
        };
    };
}
void CBoundaryVolumeMesh::GeneVolBNode()
{
    uiint numOfVol= mvBVolume.size();
    CBoundaryVolume *pBVol;
    CElement *pElem;
    CNode *pNode;
    uiint countID= mvBNode.size() + mvBEdgeBNode.size() + mvBFaceBNode.size();
    uiint ivol;
    for(ivol=0; ivol < numOfVol; ivol++) {
        pBVol= mvBVolume[ivol];
        CBoundaryNode *pBVolBNode= new CBoundaryNode;
        pBVolBNode->setMGLevel(mMGLevel+1);
        pBVolBNode->resizeValue(mMaxMGLevel-mMGLevel);
        pBVolBNode->setID(countID);
        countID++;
        pBVol->setVolBNode(pBVolBNode);
        pElem = pBVol->getElement();
        pNode = pElem->getVolumeNode();
        pBVolBNode->setNode(pNode);
        mvBVolBNode.push_back(pBVolBNode);
    };
}
void CBoundaryVolumeMesh::refine(CBoundaryVolumeMesh* pProgBVolMesh)
{
    CBoundaryVolume *pBVol;
    uiint ivol, numOfVol= mvBVolume.size();
    vector<CBoundaryVolume*> vProgVol;
    uiint countID(0);
    for(ivol=0; ivol < numOfVol; ivol++) {
        pBVol= mvBVolume[ivol];
        pBVol->refine(countID, mvDOF);
        vProgVol.clear();
        vProgVol= pBVol->getProgParts();
        uiint i;
        for(i=0; i < vProgVol.size(); i++) {
            pProgBVolMesh->addBVolume(vProgVol[i]);
        };
    };
    uiint numOfBNode    = mvBNode.size();
    uiint numOfFaceBNode= mvBFaceBNode.size();
    uiint numOfVolBNode = mvBVolBNode.size();
    uiint numOfProgBNode= numOfBNode + mnEdgeNodeCount + numOfFaceBNode + numOfVolBNode;
    pProgBVolMesh->resizeBNode(numOfProgBNode);
    uiint ibnode;
    uiint init = 0;
    uiint end  = numOfBNode;
    for(ibnode=init; ibnode < end; ibnode++) {
        pProgBVolMesh->setBNode(ibnode, mvBNode[ibnode]);
    };
    init = numOfBNode;
    end  = numOfBNode + mnEdgeNodeCount;
    for(ibnode=init; ibnode < end; ibnode++) {
        pProgBVolMesh->setBNode(ibnode, mvBEdgeBNode[ibnode-init]);
    };
    init = end;
    end += numOfFaceBNode;
    for(ibnode=init; ibnode < end; ibnode++) {
        pProgBVolMesh->setBNode(ibnode, mvBFaceBNode[ibnode-init]);
    };
    init = end;
    end += numOfVolBNode;
    for(ibnode=init; ibnode < end; ibnode++) {
        pProgBVolMesh->setBNode(ibnode, mvBVolBNode[ibnode-init]);
    }
}
void CBoundaryVolumeMesh::distNeumannValue()
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    if(mnBndType != BoundaryType::Neumann) {
        pLogger->Info(Utility::LoggerMode::Error, "BoundaryType Error,  CBoundaryVolumeMesh::distNeumannValue");
        return;
    }
    CShapeHexa  *pShHexa  = CShapeHexa::Instance();
    CShapeTetra *pShTetra = CShapeTetra::Instance();
    CShapePrism *pShPrism = CShapePrism::Instance();
    uiint inode, numOfBNode=mvBNode.size();
    CBoundaryNode *pBNode;
    uiint dof, idof;
    for(inode=0; inode < numOfBNode; inode++) {
        pBNode= mvBNode[inode];
        uiint numOfDOF = getNumOfDOF();
        for(idof=0; idof < numOfDOF; idof++) {
            dof = getDOF(idof);
            pBNode->initValue(dof, mMGLevel);//--- ノイマン値を加算するため初期化(0.0)
        };
    };
    double integVal, nodalVal, entVal;
    uiint ivol, numOfVol=mvBVolume.size();
    uiint numOfDOF;

    CCalc* pCalc=CCalc::Instance();
    double cx,cy,cz;//-------------面中心座標
    double calcVal;

    for(ivol=0; ivol < numOfVol; ivol++) {
        CBoundaryVolume *pBVol = mvBVolume[ivol];
        numOfDOF = getNumOfDOF();
        for(idof=0; idof < numOfDOF; idof++) {
            dof = getDOF(idof);
            entVal = pBVol->getBndValue(dof);

            if( existPoland(dof) ) { //---数式処理の有無
                cx= pBVol->getCenterX();
                cy= pBVol->getCenterY();
                cz= pBVol->getCenterZ();

                pCalc->setElementParam(entVal, cx,cy,cz);
                calcVal= pCalc->Exec( mmPoland[dof] );//------数式処理(idof)
            } else {
                calcVal=entVal;
            }

            uiint ivert;
            switch(pBVol->getElemType()) {
            case(ElementType::Hexa):
                for(ivert=0; ivert < 8; ivert++) {
                    integVal= pShHexa->getIntegralValue8(ivert);
                    //nodalVal= entVal * integVal;
                    nodalVal= integVal * calcVal;//-------------数式結果を形状関数分配
                    pBNode= pBVol->getBNode(ivert);
                    pBNode->addValue(dof, mMGLevel, nodalVal);
                };
                break;
            case(ElementType::Hexa2):
                for(ivert=0; ivert < 20; ivert++) {
                    integVal= pShHexa->getIntegralValue20(ivert);
                    //nodalVal= entVal * integVal;
                    nodalVal= integVal * calcVal;//-------------数式結果を形状関数分配
                    pBNode= pBVol->getBNode(ivert);
                    pBNode->addValue(dof, mMGLevel, nodalVal);
                };
                break;
            case(ElementType::Tetra):
                for(ivert=0; ivert < 4; ivert++) {
                    integVal= pShTetra->getIntegValue4(ivert);
                    //nodalVal= entVal * integVal;
                    nodalVal= integVal * calcVal;//-------------数式結果を形状関数分配
                    pBNode= pBVol->getBNode(ivert);
                    pBNode->addValue(dof, mMGLevel, nodalVal);
                }
                break;
            case(ElementType::Tetra2):
                for(ivert=0; ivert < 10; ivert++) {
                    integVal= pShTetra->getIntegValue10(ivert);
                    //nodalVal= entVal * integVal;
                    nodalVal= integVal * calcVal;//-------------数式結果を形状関数分配
                    pBNode= pBVol->getBNode(ivert);
                    pBNode->addValue(dof, mMGLevel, nodalVal);
                }
                break;
            case(ElementType::Prism):
                for(ivert=0; ivert < 6; ivert++) {
                    integVal= pShPrism->getIntegValue6(ivert);
                    //nodalVal= entVal * integVal;
                    nodalVal= integVal * calcVal;//-------------数式結果を形状関数分配
                    pBNode= pBVol->getBNode(ivert);
                    pBNode->addValue(dof, mMGLevel, nodalVal);
                }
                break;
            case(ElementType::Prism2):
                for(ivert=0; ivert < 15; ivert++) {
                    integVal= pShPrism->getIntegValue15(ivert);
                    //nodalVal= entVal * integVal;
                    nodalVal= integVal * calcVal;//-------------数式結果を形状関数分配
                    pBNode= pBVol->getBNode(ivert);
                    pBNode->addValue(dof, mMGLevel, nodalVal);
                }
                break;
            default:
                pLogger->Info(Utility::LoggerMode::Error, "CBoundaryVolumeMesh::distNeumannValue, invalid ElementType");
                break;
            }
        };
    };
}
void CBoundaryVolumeMesh::distDirichletValue()
{
    uiint ivol, numOfVol=mvBVolume.size();
    CBoundaryVolume *pBVol;
    for(ivol=0; ivol < numOfVol; ivol++) {
        pBVol = mvBVolume[ivol];
        uiint idof, dof;
        for(idof=0; idof < getNumOfDOF(); idof++) {
            dof = getDOF(idof);
            pBVol->distDirichletVal(dof, mMGLevel, mMaxMGLevel, mvPolandDOF, mmPoland);
        };
    };
}
void CBoundaryVolumeMesh::deleteProgData()
{
    uiint ivol, nNumOfBVol = mvBVolume.size();
    for(ivol=0; ivol < nNumOfBVol; ivol++) {
        mvBVolume[ivol]->deleteProgData();
    };
    vector<CBoundaryNode*>().swap(mvBEdgeBNode);
    vector<CBoundaryNode*>().swap(mvBFaceBNode);
    vector<CBoundaryNode*>().swap(mvBVolBNode);
}
