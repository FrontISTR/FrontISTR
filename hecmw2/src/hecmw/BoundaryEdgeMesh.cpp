/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryEdgeMesh.cpp
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
#include "BoundaryEdgeMesh.h"
#include "BoundaryNode.h"
using namespace pmw;
CBoundaryEdgeMesh::CBoundaryEdgeMesh()
{
    ;
}
CBoundaryEdgeMesh::~CBoundaryEdgeMesh()
{
    for_each(mvBEdge.begin(), mvBEdge.end(), DeleteObject());
}
void CBoundaryEdgeMesh::resizeEdge(const uiint& res_size)
{
    mvBEdge.resize(res_size);
}
void CBoundaryEdgeMesh::setBEdge(const uiint& index, CBoundaryEdge* pBEdge)
{
    uiint id;
    id= pBEdge->getID();
    mmBEdgeID2Index[id]= index;
    mvBEdge[index]= pBEdge;
}
void CBoundaryEdgeMesh::addBEdge(CBoundaryEdge* pBEdge)
{
    mvBEdge.push_back(pBEdge);
    uiint id;
    id= pBEdge->getID();
    mmBEdgeID2Index[id]= mvBEdge.size()-1;
}
void CBoundaryEdgeMesh::resizeAggEdge()
{
    uiint res_size = mvBNode.size();
    mvAggregateEdge.resize(res_size);
}
void CBoundaryEdgeMesh::setupAggEdge()
{
    uiint numOfBNode= mvBNode.size();
    CBoundaryNode *pBNode;
    uiint ibnode;
    for(ibnode=0; ibnode < numOfBNode; ibnode++) {
        pBNode= mvBNode[ibnode];
        pBNode->clearAggElemID();
        pBNode->clearNeibElemVert();
    };
    uiint numOfEdge= mvBEdge.size();
    CBoundaryEdge *pBEdge;
    uiint iedge;
    for(iedge=0; iedge < numOfEdge; iedge++) {
        pBEdge= mvBEdge[iedge];
        pBEdge->setupVertexElemID();
    };
    for(ibnode=0; ibnode < numOfBNode; ibnode++) {
        pBNode= mvBNode[ibnode];
        uiint numOfAgg= pBNode->getNumOfAggElem();
        uiint iagg;
        for(iagg=0; iagg < numOfAgg; iagg++) {
            mvAggregateEdge[ibnode].push_back(pBNode->getAggElemID(iagg));
        };
    };
}
void CBoundaryEdgeMesh::GeneEdgeBNode()
{
    uiint countID= mvBNode.size();
    uiint numOfEdge= mvBEdge.size();
    uiint iedge;
    CBoundaryEdge *pBEdge;
    mvBEdgeBNode.reserve(numOfEdge);
    CBoundaryNode *pBNode;
    for(iedge=0; iedge < numOfEdge; iedge++) {
        pBEdge= mvBEdge[iedge];
        if(!pBEdge->getEdgeBNode()) {
            pBNode = new CBoundaryNode;
            countID += iedge;
            pBNode->setID(countID);
            pBEdge->setEdgeBNode(pBNode);
            mvBEdgeBNode.push_back(pBNode);
            if(pBEdge->getOrder()==ElementOrder::Second) {
                mvBNode.push_back(pBNode);
                mmBNodeID2Index[pBNode->getID()] = mvBNode.size()-1;
                uiint index = mvBNode.size()-1;
                mvAggregateEdge.resize(mvBNode.size());
                mvAggregateEdge[index].push_back(pBEdge->getID());
                pBNode->setMGLevel(mMGLevel);
                pBNode->resizeValue(mMaxMGLevel-mMGLevel + 1);
            } else {
                pBNode->setMGLevel(mMGLevel+1);
                pBNode->resizeValue(mMaxMGLevel-mMGLevel);
                mnEdgeNodeCount++;
            }
            pBEdge->setupNode();
            if(pBEdge->getOrder()==ElementOrder::Second)  pBEdge->replaceEdgeBNode();
        }
    };
}
void CBoundaryEdgeMesh::refine(CBoundaryEdgeMesh* pProgEdgeMesh)
{
    CBoundaryEdge *pBEdge;
    vector<CBoundaryEdge*> vBEdge;
    uiint countID(0);
    uiint iedge, numOfBEdge= mvBEdge.size();
    for(iedge=0; iedge < numOfBEdge; iedge++) {
        pBEdge= mvBEdge[iedge];
        pBEdge->refine(countID, mvDOF);
        vBEdge.clear();
        vBEdge= pBEdge->getProgParts();
        uiint iprog;
        for(iprog=0; iprog < vBEdge.size(); iprog++) {
            pProgEdgeMesh->addBEdge(vBEdge[iprog]);
        };
    };
    uiint numOfBNode= mvBNode.size();
    uiint numOfProgBNode= numOfBNode + mnEdgeNodeCount;
    pProgEdgeMesh->resizeBNode(numOfProgBNode);
    uiint ibnode;
    CBoundaryNode *pBNode;
    for(ibnode=0; ibnode < numOfBNode; ibnode++) {
        pBNode= mvBNode[ibnode];
        pProgEdgeMesh->setBNode(ibnode, pBNode);
    };
    for(ibnode=numOfBNode; ibnode < numOfProgBNode; ibnode++) {
        pBNode= mvBEdgeBNode[ibnode-numOfBNode];
        pProgEdgeMesh->setBNode(ibnode, pBNode);
    };
}
void CBoundaryEdgeMesh::distNeumannValue()
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    if(mnBndType != BoundaryType::Neumann) {
        pLogger->Info(Utility::LoggerMode::Error, "BoundaryType Error,  CBoundaryEdgeMesh::distNeumannValue");
        return;
    }
    uiint idof, dof, numOfDOF;
    uiint inode,numOfBNode=mvBNode.size();
    CBoundaryNode *pBNode;
    for(inode=0; inode < numOfBNode; inode++) {
        pBNode= mvBNode[inode];
        numOfDOF= getNumOfDOF();
        for(idof=0; idof < numOfDOF; idof++) {
            dof= getDOF(idof);
            pBNode->initValue(dof, mMGLevel);//--- ノイマン値を加算するため初期化(0.0)
        };
    };
    CShapeLine *pShLine = CShapeLine::Instance();
    CBoundaryEdge *pBEdge;
    uiint  iedge, ivert, numOfEdge = mvBEdge.size();
    double entVal,integVal,nodalVal;

    CCalc* pCalc=CCalc::Instance();
    double cx,cy,cz;//-------------面中心座標
    double calcVal;

    for(iedge=0; iedge < numOfEdge; iedge++) {
        pBEdge = mvBEdge[iedge];
        uiint numOfDOF = getNumOfDOF();
        for(idof=0; idof < numOfDOF; idof++) {
            dof = getDOF(idof);
            entVal = pBEdge->getBndValue(dof);

            if( existPoland(dof) ) { //---数式処理の有無
                cx= pBEdge->getCenterX();
                cy= pBEdge->getCenterY();
                cz= pBEdge->getCenterZ();

                pCalc->setElementParam(entVal, cx,cy,cz);
                calcVal= pCalc->Exec( mmPoland[dof] );//------数式処理(idof)
            } else {
                calcVal=entVal;
            }

            switch(pBEdge->getBEdgeShape()) {
            case(ElementType::Beam):
                for(ivert=0; ivert < 2; ivert++) {
                    integVal= pShLine->getIntegValue2(ivert);
                    //nodalVal= integVal * entVal;
                    nodalVal= integVal * calcVal;//-------------数式結果を形状関数分配
                    pBNode= pBEdge->getBNode(ivert);
                    pBNode->addValue(dof, mMGLevel, nodalVal);
                };
                break;
            case(ElementType::Beam2):
                for(ivert=0; ivert < 3; ivert++) {
                    integVal= pShLine->getIntegValue3(ivert);
                    //nodalVal= integVal * entVal;
                    nodalVal= integVal * calcVal;//-------------数式結果を形状関数分配
                    pBNode= pBEdge->getBNode(ivert);
                    pBNode->addValue(dof, mMGLevel, nodalVal);
                };
                break;
            default:
                break;
            }
        };
    };
}
void CBoundaryEdgeMesh::distDirichletValue()
{
    uiint iedge, numOfEdge=mvBEdge.size();
    CBoundaryEdge *pBEdge;
    for(iedge=0; iedge < numOfEdge; iedge++) {
        pBEdge = mvBEdge[iedge];
        uiint idof, dof;
        for(idof=0; idof < getNumOfDOF(); idof++) {
            dof = getDOF(idof);
            pBEdge->distDirichletVal(dof, mMGLevel, mMaxMGLevel, mvPolandDOF, mmPoland);
        };
    };
}
void CBoundaryEdgeMesh::deleteProgData()
{
    vector<CBoundaryNode*>().swap(mvBEdgeBNode);
}
