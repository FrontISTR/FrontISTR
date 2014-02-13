/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryFaceMesh.cpp
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
#include "BoundaryNode.h"
#include "BoundaryParts.h"
#include "BoundaryFace.h"
#include "BoundaryFaceMesh.h"
using namespace pmw;
CBoundaryFaceMesh::CBoundaryFaceMesh()
{
    ;
}
CBoundaryFaceMesh::~CBoundaryFaceMesh()
{
    for_each(mvBFace.begin(), mvBFace.end(), DeleteObject());
}
void CBoundaryFaceMesh::setBFace(const uiint& index, CBoundaryFace* pBFace)
{
    uiint id;
    id= pBFace->getID();
    mmBFaceID2Index[id]= index;
    mvBFace[index]= pBFace;
}
void CBoundaryFaceMesh::addBFace(CBoundaryFace* pBFace)
{
    mvBFace.push_back(pBFace);
    uiint id;
    id= pBFace->getID();
    mmBFaceID2Index[id]= mvBFace.size()-1;
}
void CBoundaryFaceMesh::resizeAggFace()
{
    uiint numOfAgg= mvBNode.size();
    mvAggregateFace.resize(numOfAgg);
}
void CBoundaryFaceMesh::setupAggFace()
{
    int numOfBNode= mvBNode.size();
    CBoundaryNode *pBNode;
    uiint ibnode;
    for(ibnode=0; ibnode < numOfBNode; ibnode++) {
        pBNode= mvBNode[ibnode];
        pBNode->clearAggElemID();
        pBNode->clearNeibElemVert();
    };
    uiint numOfFace= mvBFace.size();
    CBoundaryFace *pBFace;
    uiint iface;
    for(iface=0; iface < numOfFace; iface++) {
        pBFace= mvBFace[iface];
        uiint numOfVert= pBFace->getNumOfVert();
        uiint ivert, faceID;
        for(ivert=0; ivert < numOfVert; ivert++) {
            pBNode= pBFace->getBNode(ivert);
            faceID= pBFace->getID();
            pBNode->setAggElemID(faceID);
        };
    };
    for(ibnode=0; ibnode < numOfBNode; ibnode++) {
        pBNode= mvBNode[ibnode];
        uiint numOfAgg= pBNode->getNumOfAggElem();
        uiint iagg;
        for(iagg=0; iagg < numOfAgg; iagg++) {
            mvAggregateFace[ibnode].push_back(pBNode->getAggElemID(iagg));
        };
    };
}
void CBoundaryFaceMesh::GeneEdgeBNode()
{
    uiint countID= mvBNode.size();
    uiint numOfBFace= mvBFace.size();
    uiint iface;
    CBoundaryFace *pBFace;
    bool bfind(false);
    for(iface=0; iface < numOfBFace; iface++) {
        pBFace= mvBFace[iface];
        uiint numOfEdge= pBFace->getNumOfEdge();
        PairBNode pairBNode;
        uiint iedge;
        for(iedge=0; iedge < numOfEdge; iedge++) {
            if(!pBFace->isMarkingEdge(iedge)) {
                pairBNode= pBFace->getPairBNode(iedge);
                uiint numOfiAgg= pairBNode.first->getNumOfAggElem();
                uiint numOfjAgg= pairBNode.second->getNumOfAggElem();
                uiint iAgg,jAgg;
                uiint iFaceID,jFaceID;
                bfind= false;
                for(iAgg=0; iAgg < numOfiAgg; iAgg++) {
                    iFaceID= pairBNode.first->getAggElemID(iAgg);
                    for(jAgg=0; jAgg < numOfjAgg; jAgg++) {
                        jFaceID= pairBNode.second->getAggElemID(jAgg);
                        if(iFaceID == jFaceID) {
                            if(iFaceID != pBFace->getID()) {
                                bfind= true;
                                CBoundaryNode *pEdgeBNode= new CBoundaryNode;//--------------------生成
                                pEdgeBNode->setID(countID);
                                countID++;
                                mvBEdgeBNode.push_back(pEdgeBNode);
                                pBFace->setEdgeNeibFace(iedge,iFaceID);
                                pBFace->markingEdge(iedge);
                                pBFace->setEdgeBNode(iedge, pEdgeBNode);
                                uiint neib_index= mmBFaceID2Index[iFaceID];
                                CBoundaryFace* pNeibBFace= mvBFace[neib_index];
                                uiint jedge= pNeibBFace->getEdgeID(pairBNode);
                                pNeibBFace->setEdgeNeibFace(jedge, pBFace->getID());
                                pNeibBFace->markingEdge(jedge);
                                pNeibBFace->setEdgeBNode(jedge, pEdgeBNode);
                                if(pBFace->getOrder()==ElementOrder::Second) {
                                    mvBNode.push_back(pEdgeBNode);
                                    mmBNodeID2Index[pEdgeBNode->getID()] = mvBNode.size()-1;
                                    uiint index = mvBNode.size()-1;
                                    mvAggregateFace.resize(mvBNode.size());
                                    mvAggregateFace[index].push_back(pBFace->getID());
                                    mvAggregateFace[index].push_back(pNeibBFace->getID());
                                    pBFace->replaceEdgeBNode();
                                    pNeibBFace->replaceEdgeBNode();
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
                    CBoundaryNode *pEdgeBNode= new CBoundaryNode;//----------------------生成
                    pEdgeBNode->setID(countID);
                    countID++;
                    mvBEdgeBNode.push_back(pEdgeBNode);
                    pBFace->setEdgeBNode(iedge, pEdgeBNode);
                    pBFace->markingEdge(iedge);
                    if(pBFace->getOrder()==ElementOrder::Second) {
                        mvBNode.push_back(pEdgeBNode);
                        mmBNodeID2Index[pEdgeBNode->getID()] = mvBNode.size()-1;
                        uiint index = mvBNode.size()-1;
                        mvAggregateFace.resize(mvBNode.size());
                        mvAggregateFace[index].push_back(pBFace->getID());
                        pBFace->replaceEdgeBNode();
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
        pBFace->setupNode_Edge();
    };
}
void CBoundaryFaceMesh::GeneFaceBNode()
{
    uiint numOfFace= mvBFace.size();
    mvBFaceBNode.reserve(numOfFace);
    CBoundaryFace *pBFace;
    uiint iface;
    uiint countID= mvBNode.size() + mvBEdgeBNode.size();
    for(iface=0; iface < numOfFace; iface++) {
        CBoundaryNode *pFaceBNode= new CBoundaryNode;//------------------生成
        mvBFaceBNode.push_back(pFaceBNode);
        pFaceBNode->setMGLevel(mMGLevel+1);
        pFaceBNode->resizeValue(mMaxMGLevel-mMGLevel);
        pFaceBNode->setID(countID);
        countID++;
        pBFace= mvBFace[iface];
        pBFace->setFaceBNode(pFaceBNode);
        pBFace->setupNode_Face();
    };
}
void CBoundaryFaceMesh::refine(CBoundaryFaceMesh *pProgBFaceMesh)
{
    CBoundaryFace *pBFace;
    uiint iface, numOfFace= mvBFace.size();
    uiint countID(0);
    vector<CBoundaryFace*> vprogBFace;
    for(iface=0; iface < numOfFace; iface++) {
        pBFace= mvBFace[iface];
        pBFace->refine(countID, mvDOF);
        vprogBFace.clear();
        vprogBFace= pBFace->getProgParts();
        uiint i;
        for(i=0; i < vprogBFace.size(); i++) {
            pProgBFaceMesh->addBFace(vprogBFace[i]);
        };
    };
    uiint numOfBNode    = mvBNode.size();
    uiint numOfFaceBNode= mvBFaceBNode.size();
    uiint numOfProgBNode= numOfBNode + mnEdgeNodeCount + numOfFaceBNode;
    pProgBFaceMesh->resizeBNode(numOfProgBNode);
    uiint ibnode;
    uiint init= 0;
    uiint end = numOfBNode;
    for(ibnode=init; ibnode < end; ibnode++) {
        pProgBFaceMesh->setBNode(ibnode, mvBNode[ibnode]);
    };
    init= end;
    end = init + mnEdgeNodeCount;
    for(ibnode=init; ibnode < end; ibnode++) {
        pProgBFaceMesh->setBNode(ibnode, mvBEdgeBNode[ibnode - init]);
    };
    init= end;
    end = init + numOfFaceBNode;
    for(ibnode=init; ibnode < end; ibnode++) {
        pProgBFaceMesh->setBNode(ibnode, mvBFaceBNode[ibnode - init]);
    }
}
void CBoundaryFaceMesh::distNeumannValue()
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    if(mnBndType != BoundaryType::Neumann) {
        pLogger->Info(Utility::LoggerMode::Error, "BoundaryType Error,  CBoundaryFaceMesh::distNeumannValue");
        return;
    }
    CShapeQuad     *pQuadShape= CShapeQuad::Instance();
    CShapeTriangle *pTriShape= CShapeTriangle::Instance();
    CBoundaryNode *pBndNode;
    CBoundaryFace *pBndFace;
    double entVal;
    double integVal;
    double nodalVal;
    uiint ivert, inode, numOfBNode= mvBNode.size();
    uiint iface, numOfFace= mvBFace.size();
    uiint idof, numOfDOF, dof;
    for(inode=0; inode < numOfBNode; inode++) {
        pBndNode= mvBNode[inode];
        numOfDOF = getNumOfDOF();
        for(idof=0; idof < numOfDOF; idof++) {
            dof = getDOF(idof);
            pBndNode->initValue(dof, mMGLevel);//--- ノイマン値を加算するため初期化(0.0)
        };
    };

    CCalc* pCalc=CCalc::Instance();
    double cx,cy,cz;//-------------面中心座標
    double calcVal;
    // Face loop
    for(iface=0; iface < numOfFace; iface++) {
        pBndFace= mvBFace[iface];
        numOfDOF = getNumOfDOF();
        // DOF loop
        for(idof=0; idof < numOfDOF; idof++) {
            dof = getDOF(idof);
            entVal= pBndFace->getBndValue(dof);//---------------面の境界値

            if( existPoland(dof) ) { //---数式処理の有無
                cx= pBndFace->getCenterX();
                cy= pBndFace->getCenterY();
                cz= pBndFace->getCenterZ();

                pCalc->setElementParam(entVal, cx,cy,cz);
                calcVal= pCalc->Exec( mmPoland[dof] );//------数式処理(dof)
            } else {
                calcVal=entVal;
            }

            ////cout << "BoundaryFaceMesh::distNeumannValue  dof:" << dof << " calcVal:" << calcVal
            ////       << " x:" << cx << " y:" << cy << " z:" << cz << endl;


            switch(pBndFace->getBFaceShape()) {
            case(ElementType::Quad):
                for(ivert=0; ivert < 4; ivert++) {
                    integVal= pQuadShape->getIntegValue4(ivert);
                    //nodalVal= entVal * integVal;
                    nodalVal= calcVal * integVal;//-----------------数式結果を形状関数分配
                    pBndNode= pBndFace->getBNode(ivert);

                    ////cout << "BoundaryFaceMesh::distNeumannValue  nodalVal:" << nodalVal
                    ///        << " dof:" << dof << " NodeID:" << pBndNode->getNode()->getID() << endl;

                    pBndNode->addValue(dof, mMGLevel, nodalVal);
                };
                break;
            case(ElementType::Quad2):
                for(ivert=0; ivert < 8; ivert++) {
                    integVal= pQuadShape->getIntegValue8(ivert);
                    //nodalVal= entVal * integVal;
                    nodalVal= calcVal * integVal;//-----------------数式結果を形状関数分配
                    pBndNode= pBndFace->getBNode(ivert);
                    pBndNode->addValue(dof, mMGLevel, nodalVal);
                };
                break;
            case(ElementType::Triangle):
                for(ivert=0; ivert < 3; ivert++) {
                    integVal= pTriShape->getIntegValue3(ivert);
                    //nodalVal= entVal * integVal;
                    nodalVal= calcVal * integVal;//-----------------数式結果を形状関数分配
                    pBndNode= pBndFace->getBNode(ivert);
                    pBndNode->addValue(dof, mMGLevel, nodalVal);
                };
                break;
            case(ElementType::Triangle2):
                for(ivert=0; ivert < 6; ivert++) {
                    integVal= pTriShape->getIntegValue6(ivert);
                    //nodalVal= entVal * integVal;
                    nodalVal= calcVal * integVal;//-----------------数式結果を形状関数分配
                    pBndNode= pBndFace->getBNode(ivert);
                    pBndNode->addValue(dof, mMGLevel, nodalVal);
                };
                break;
            default:
                break;
            }
        };
    };
}
void CBoundaryFaceMesh::distDirichletValue()
{
    uiint iface, numOfFace=mvBFace.size();
    CBoundaryFace *pBFace;
    for(iface=0; iface < numOfFace; iface++) {
        pBFace = mvBFace[iface];
        uiint idof, dof;
        for(idof=0; idof < getNumOfDOF(); idof++) {
            dof = getDOF(idof);
            pBFace->distDirichletVal(dof, mMGLevel, mMaxMGLevel, mvPolandDOF, mmPoland);
        };
    };
}
void CBoundaryFaceMesh::deleteProgData()
{
    uiint iface, nNumOfBFace=mvBFace.size();
    for(iface=0; iface < nNumOfBFace; iface++) mvBFace[iface]->deleteProgData();
    vector<CBoundaryNode*>().swap(mvBEdgeBNode);
    vector<CBoundaryNode*>().swap(mvBFaceBNode);
}
