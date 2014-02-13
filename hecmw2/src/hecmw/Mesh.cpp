/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Mesh.cpp
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
#include "Vertex.h"
#include "LoggerType.h"
#include "IndexBucket.h"
#include "BoundaryEdgeMesh.h"
#include "BoundaryNodeMesh.h"
#include "BoundaryMesh.h"
#include "CommunicationMesh.h"
#include "Logger.h"
#include "Node.h"
#include <vector>
#include "VectorNode.h"
#include "Element.h"
#include "AggregateElement.h"
#include "Mesh.h"
#include "ScalarNode.h"
#include "ScalarVectorNode.h"
using namespace pmw;
CMesh::CMesh(void)
{
    mpLogger = Utility::CLogger::Instance();
}
CMesh::CMesh(const uiint& numofnode, const uiint& numofelem)
{
    mpLogger = Utility::CLogger::Instance();
    mNumOfNode = numofnode;
    mNumOfElement = numofelem;
    mvNode.resize(mNumOfNode);
    mvElement.resize(mNumOfElement);

    mvMarkingBNode= NULL;
}
CMesh::~CMesh(void)
{
    if(mnSolutionType==SolutionType::FEM)
        for_each(mvAggElement.begin(), mvAggElement.end(), DeleteObject());

    if(mnSolutionType==SolutionType::FVM)
        for_each(mvAggNode.begin(), mvAggNode.end(), DeleteObject());

    for_each(mvCommMesh.begin(), mvCommMesh.end(), DeleteObject());

    if(mMGLevel==mMaxMGLevel)
        for_each(mvNode.begin(), mvNode.end(), DeleteObject());

    for_each(mvElement.begin(), mvElement.end(), DeleteObject());

    if(mvMarkingBNode)
        delete[] mvMarkingBNode;
    if(mvMarkingDirichlet)
        delete[] mvMarkingDirichlet;
    if(mvMarkingNeumann)
        delete[] mvMarkingNeumann;

    if(mvMarkingLargeRankNode)
        delete[] mvMarkingLargeRankNode;

}
void CMesh::initBucketNode(const uiint& max_id, const uiint& min_id)
{
    moBucket.clearBucketNode();
    moBucket.resizeBucketNode(max_id, min_id);
}
void CMesh::setupBucketNodeIndex(const uiint& id, const uiint& index)
{
    moBucket.setIndexNode(id, index);
}
void CMesh::setupBucketNode()
{
    CNode *pNode;
    uiint maxID,minID, i;
    pNode = mvNode[0];
    maxID = pNode->getID();
    minID = pNode->getID();
    for(i=0; i < mvNode.size(); i++) {
        pNode = mvNode[i];
        if(pNode->getID() > maxID) maxID = pNode->getID();
        if(pNode->getID() < minID) minID = pNode->getID();
    };
    moBucket.resizeBucketNode(maxID, minID);
    for(i=0; i < mvNode.size(); i++) {
        pNode = mvNode[i];
        moBucket.setIndexNode(pNode->getID(), i);
    };
}
void CMesh::initBucketElement(const uiint& max_id, const uiint& min_id)
{
    moBucket.clearBucketElement();
    moBucket.resizeBucketElement(max_id, min_id);
}
void CMesh::setupBucketElementIndex(const uiint& id, const uiint& index)
{
    moBucket.setIndexElement(id, index);
}
void CMesh::setupBucketElement()
{
    CElement *pElement;
    uiint maxID,minID, i;
    pElement = mvElement[0];
    maxID = pElement->getID();
    minID = pElement->getID();
    for(i=0; i < mvElement.size(); i++) {
        pElement = mvElement[i];
        if(pElement->getID() > maxID) maxID = pElement->getID();
        if(pElement->getID() < minID) minID = pElement->getID();
    };
    moBucket.resizeBucketElement(maxID,minID);
    for(i=0; i < mvElement.size(); i++) {
        pElement = mvElement[i];
        moBucket.setIndexElement(pElement->getID(), i);
    };
}
void CMesh::reserveNode(const uiint& num_of_node)
{
    mvNode.reserve(num_of_node);
}
void CMesh::setNode(CNode *pNode)
{
    mvNode.push_back(pNode);
}
CNode* CMesh::getNode(const uiint &id)
{
    uiint index;
    index = moBucket.getIndexNode(id);
    return mvNode[index];
}
CNode* CMesh::getNodeIX(const uiint& index)
{
    return mvNode[index];
}
void CMesh::reserveElement(const uiint& num_of_elem)
{
    mvElement.reserve(num_of_elem);
}
void CMesh::setElement(CElement *pElement)
{
    mvElement.push_back(pElement);
}
CElement* CMesh::getElement(const uiint& id)
{
    uiint index;
    index = moBucket.getIndexElement(id);
    return mvElement[index];
}
CElement* CMesh::getElementIX(const uiint& index)
{
    return mvElement[index];
}
void CMesh::resizeAggregate(const uiint& res_size)
{
    if(mnSolutionType==SolutionType::FEM) mvAggElement.resize(res_size);
    if(mnSolutionType==SolutionType::FVM) mvAggNode.resize(res_size);
}
void CMesh::setAggElement(CAggregateElement* pAggElem, const uiint& inode)
{
    mvAggElement[inode] = pAggElem;
}
void CMesh::setAggNode(CAggregateNode* pAggNode, const uiint& inode)
{
    mvAggNode[inode] = pAggNode;
}
CAggregateElement* CMesh::getAggElem(const uiint& node_id)
{
    uiint index;
    index = moBucket.getIndexNode(node_id);
    return mvAggElement[index];
}
CAggregateElement* CMesh::getAggElemIX(const uiint& inode)
{
    return mvAggElement[inode];
}
CAggregateNode* CMesh::getAggNode(const uiint& node_id)
{
    uiint index;
    index = moBucket.getIndexNode(node_id);
    return mvAggNode[index];
}
CAggregateNode* CMesh::getAggNodeIX(const uiint& inode)
{
    return mvAggNode[inode];
}
void CMesh::setupAggregate(const uiint& nLevel)
{
    CElement *pElem;
    CNode *pNode;
    uiint ielem, local_id, inode;
    mNumOfNode= mvNode.size();
    mNumOfElement= mvElement.size();
    for(inode=0; inode< mNumOfNode; inode++) {
        pNode= mvNode[inode];
        pNode->clearAggElemID();
        pNode->clearNeibElemVert();
    };
    uiint numOfLocalNode;
    uiint elemID;
    for(ielem=0; ielem< mNumOfElement; ielem++) {
        pElem = mvElement[ielem];
        if(nLevel==0 && pElem->getOrder()==ElementOrder::Second) {
            numOfLocalNode= pElem->getNumOfNode();
            elemID = pElem->getID();
            for(local_id=0; local_id < numOfLocalNode; local_id++) {
                pNode = pElem->getNode(local_id);
                pNode->setAggElemID(elemID);
                pNode->setNeibElemVert(elemID, local_id);
            };
        } else {
            numOfLocalNode= pElem->getNumOfVert();
            elemID = pElem->getID();
            for(local_id=0; local_id< numOfLocalNode; local_id++) {
                pNode = pElem->getNode(local_id);
                pNode->setAggElemID(elemID);
                pNode->setNeibElemVert(elemID,local_id);
            };
        }
    };

    if(mnSolutionType==SolutionType::FEM) {
        CAggregateElement *pAggElem;
        uiint iagg_elem, numOfAggElem;
        for(inode=0; inode< mNumOfNode; inode++) {
            pNode = mvNode[inode];
            numOfAggElem = pNode->getNumOfAggElem();
            pAggElem = mvAggElement[inode];
            pAggElem->reserve(numOfAggElem);
            for(iagg_elem=0; iagg_elem< numOfAggElem; iagg_elem++) {
                elemID = pNode->getAggElemID(iagg_elem);
                pElem= getElement(elemID);
                pAggElem->push(pElem);
            };
        };
    }
    if(mnSolutionType==SolutionType::FVM) {
        CAggregateNode *pAggNode;
        vector<CNode*> vConnNode;
        uiint ncon, nNumOfConnNode;
        uiint iagg_elem, numOfAggElem;
        for(inode=0; inode< mNumOfNode; inode++) {
            pNode= mvNode[inode];
            numOfAggElem= pNode->getNumOfAggElem();
            pAggNode = mvAggNode[inode];
            pAggNode->reserveNode(numOfAggElem);
            for(iagg_elem=0; iagg_elem< numOfAggElem; iagg_elem++) {
                elemID= pNode->getAggElemID(iagg_elem);
                pElem= getElement(elemID);
                vConnNode= pElem->getConnectNode(pNode);
                nNumOfConnNode= vConnNode.size();
                for(ncon=0; ncon < nNumOfConnNode; ncon++) {
                    pAggNode->setNode(vConnNode[ncon]);
                };
            };
        };
    }
}
void CMesh::presetProgMesh(CMesh* pProgMesh)
{
    uiint progNumOfNode= mNumOfNode;
    uiint progNumOfElem= mNumOfElement;
    CElement *pElem;
    uiint ielem;
    for(ielem=0; ielem< mNumOfElement; ielem++) {
        pElem= mvElement[ielem];
        switch(pElem->getType()) {
        case(ElementType::Hexa):
        case(ElementType::Hexa2):
            progNumOfElem += 7;
            break;
        case(ElementType::Tetra):
        case(ElementType::Tetra2):
            progNumOfElem += 3;
            break;
        case(ElementType::Prism):
        case(ElementType::Prism2):
            progNumOfElem += 5;
            break;
        case(ElementType::Quad):
        case(ElementType::Quad2):
            progNumOfElem += 3;
            break;
        case(ElementType::Triangle):
        case(ElementType::Triangle2):
            progNumOfElem += 2;
            break;
        case(ElementType::Beam):
        case(ElementType::Beam2):
            progNumOfElem += 1;
            break;
        }
    };
    progNumOfNode *= 8;
    pProgMesh->reserveNode(progNumOfNode);
    pProgMesh->reserveElement(progNumOfElem);
    uiint i;
    for(i=0; i< mNumOfNode; i++) {
        pProgMesh->setNode(mvNode[i]);
    };
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::MWDebug,"prolongation numOfElem",progNumOfElem);
    pLogger->Info(Utility::LoggerMode::MWDebug,"prolongation numOfNode",progNumOfNode);
}
void CMesh::setupEdgeElement(CMesh *pProgMesh, const uiint& nLevel)
{
    CElement *pElem;
    PairNode edgeNode;
    CNode *pNode0,*pNode1;
    CElement *pOtherElem;
    vector<CElement*> vEdgeElem;
    uiint nOtherEdge;
    uiint numOfEdge;
    uiint ielem,iedge;
    uiint nInitNodeSize;
    uiint indexCount;   
    CNode *inNode;  
////    if(pProgMesh){
////        indexCount = moBucket.getMaxNodeID() + 1;
////    }else{
////        indexCount = moBucket.getMaxNodeID() + 1;
////    }
    
    indexCount = moBucket.getMaxNodeID() + 1;


    for(ielem=0; ielem< mNumOfElement; ielem++) {

        //cout << "Mesh::setupEdgeElement ielem:" << ielem << endl;

        pElem = mvElement[ielem];
        numOfEdge = pElem->getNumOfEdge();
        bool bFlag=false;
        if(pElem->getOrder()==ElementOrder::Second && nLevel > 0 ) bFlag=true;
        if(pElem->getOrder()==ElementOrder::First)  bFlag=true;

        if(bFlag){
            for(iedge=0; iedge< numOfEdge; iedge++){

                //cout << "Mesh::setupEdgeElement iedge:" << iedge << endl;

                edgeNode = pElem->getPairNode(iedge);
                pNode0= edgeNode.first;  pNode1= edgeNode.second;

                if(!pElem->isEdgeElem(pNode0,pNode1)){

                    //cout << "Mesh::setupEdgeElement !isEdgeElem" << endl;

                    uiint iagg,jagg, elemIndex0, elemIndex1, elemid0, elemid1;
                    uiint numOfAggElem0=pNode0->getNumOfAggElem();
                    uiint numOfAggElem1=pNode1->getNumOfAggElem();
                    pElem->reserveEdgeElement(iedge,numOfAggElem0);

                    //cout << "Mesh::setupEdgeElement --- NumAggElem0:" << numOfAggElem0
                    //                                << " NumAggElem1:" << numOfAggElem1 << endl;

                    for(iagg=0; iagg< numOfAggElem0; iagg++){
                        elemid0= pNode0->getAggElemID(iagg);
                        elemIndex0= moBucket.getIndexElement(elemid0);

                        //cout << "Mesh::setupEdgeElement --- iagg:" << iagg << endl;

                        for(jagg=0; jagg< numOfAggElem1; jagg++){
                            elemid1= pNode1->getAggElemID(jagg);
                            elemIndex1= moBucket.getIndexElement(elemid1);
                            if(elemIndex0 == elemIndex1)
                                pElem->setEdgeElement(iedge, mvElement[elemIndex0]);
                        };
                    };
                    
                    //cout << "Mesh::setupEdgeElement --- A" << endl;

                    inNode= GeneInterNode(pNode0);
                    vdouble P0_coord= pNode0->getCoord();
                    vdouble P1_coord= pNode1->getCoord();
                    vdouble In_coord;
                    In_coord.reserve(3);
                    for(uiint i=0; i< 3; i++) {
                        In_coord.push_back( (P0_coord[i] + P1_coord[i])*0.5 );
                    }
                    inNode->setCoord(In_coord);
                    inNode->setID(indexCount);
                    if(pProgMesh) {
                        setupParentNode(pNode0,pNode1,inNode);//-------辺ノード
                    } else {
                        ////cout << "Mesh setupEdgeElement 2nd" << endl;
                        setupParentNode2nd(pNode0,pNode1,inNode);//-----最終Level 2次要素 辺ノード
                    }

                    //cout << "Mesh::setupEdgeElement --- B" << endl;

                    pElem->setEdgeInterNode(inNode,iedge);
                    pElem->setBoolEdgeElem(pNode0, pNode1);
                    if(pProgMesh) {
                        pProgMesh->setNode(inNode);
                    } else {
                        ////cout << "Mesh setupEdgeElement NULL" << endl;
                    }

                    //cout << "Mesh::setupEdgeElement --- C" << endl;

                    if(pElem->getOrder()==ElementOrder::Second){
                        mvNode.push_back(inNode);
                        mNumOfNode += 1;
                        CAggregateElement *pTAggElem = new CAggregateElement;
                        mvAggElement.push_back( pTAggElem );
                        moBucket.re_resizeBucketNode(indexCount);
                        moBucket.setIndexNode( indexCount, mvNode.size()-1 );
                    }

                    //cout << "Mesh::setupEdgeElement --- D" << endl;

                    indexCount++;
                    vEdgeElem.clear();
                    vEdgeElem= pElem->getEdgeElement(iedge);
                    uiint n;
                    for(n=0; n < vEdgeElem.size(); n++) {
                        pOtherElem=vEdgeElem[n];
                        if(pOtherElem->getID() != pElem->getID() && !pOtherElem->isEdgeElem(pNode0,pNode1)) {
                            nOtherEdge =pOtherElem->getEdgeIndex(pNode0, pNode1);
                            pOtherElem->setEdgeAggElement(nOtherEdge,vEdgeElem);
                            pOtherElem->setBoolEdgeElem(pNode0, pNode1);
                            pOtherElem->setEdgeInterNode(inNode,nOtherEdge);
                        }
                    };

                    //cout << "Mesh::setupEdgeElement --- E" << endl;

                    if(pElem->getOrder()==ElementOrder::Second && mnSolutionType==SolutionType::FEM){
                        uiint index = mvNode.size() - 1;
                        CAggregateElement *pAggElem = mvAggElement[index];
                        pAggElem->reserve(vEdgeElem.size());
                        for(uiint i=0; i < vEdgeElem.size(); i++) {
                            CElement *pEdgeElem = vEdgeElem[i];
                            pAggElem->push(pEdgeElem);
                        };
                    }

                    //cout << "Mesh::setupEdgeElement --- F" << endl;
                    
                }//if(!isEdgeElem)end
            };
        }
    };
}
void CMesh::replaceEdgeNode()
{
    uiint ielem;
    for(ielem=0; ielem< mNumOfElement; ielem++) {
        mvElement[ielem]->replaseEdgeNode();
    };
}
CNode* CMesh::GeneInterNode(CNode* pNode)
{
    CNode *inNode;
    uiint nSDOF, nVDOF;
    switch(pNode->getType()) {
    case(NodeType::Scalar):
        inNode= new CScalarNode;
        nSDOF = pNode->getScalarDOF();
        inNode->setScalarDOF(nSDOF);
        //inNode->resizeGridLevel(mMaxMGLevel-mMGLevel);
        inNode->resizeGridLevel(mMaxMGLevel+1);
        break;
    case(NodeType::Vector):
        inNode= new CVectorNode;
        nVDOF = pNode->getVectorDOF();
        inNode->setVectorDOF(nVDOF);
        //inNode->resizeGridLevel(mMaxMGLevel-mMGLevel);
        inNode->resizeGridLevel(mMaxMGLevel+1);
        break;
    case(NodeType::ScalarVector):
        inNode= new CScalarVectorNode;
        nSDOF = pNode->getScalarDOF();
        nVDOF = pNode->getVectorDOF();
        inNode->setScalarDOF(nSDOF);
        inNode->setVectorDOF(nVDOF);
        //inNode->resizeGridLevel(mMaxMGLevel-mMGLevel);
        inNode->resizeGridLevel(mMaxMGLevel+1);
        break;
    default:
        mpLogger->Info(Utility::LoggerMode::Error,"Node Generation Error, CMesh::GeneInterNode");
        break;
    }
    return inNode;
}
void CMesh::setupFaceElement(CMesh* pProgMesh)
{
    CElement *pElem;
    vector<CNode*> vFaceCnvNode;
    uiint  elemIndex0, elemIndex1, elemIndex2, elemid0, elemid1, elemid2;
    CNode *pNode0,*pNode1,*pNode2;
    CNode *inNode;
    uiint ielem,isurf,jsurf, iagg,jagg, kagg;
    uiint numAggElems0, numAggElems1, numAggElems2;
    uiint numOfFace;
    vuint vnShareElems0, vnShareElems1;
    uiint nInitNodeSize;
    uiint indexCount;
    nInitNodeSize = pProgMesh->getNodeSize();
    CNode *pNode= pProgMesh->getNodeIX(nInitNodeSize-1);
    indexCount = pNode->getID()+1;
    for(ielem=0; ielem< mNumOfElement; ielem++) {
        pElem= mvElement[ielem];
        numOfFace= pElem->getNumOfFace();
        for(isurf=0; isurf< numOfFace; isurf++) {
            vFaceCnvNode = pElem->getFaceCnvNodes(isurf);
            pNode0= vFaceCnvNode[0];
            pNode1= vFaceCnvNode[1];
            pNode2= vFaceCnvNode[2];
            vnShareElems0.clear();
            vnShareElems1.clear();
            if(!pElem->isFaceElem(pNode0,pNode1,pNode2)) {
                numAggElems0= pNode0->getNumOfAggElem();
                vnShareElems0.reserve(numAggElems0);
                numAggElems1= pNode1->getNumOfAggElem();
                for(iagg=0; iagg< numAggElems0; iagg++) {
                    for(jagg=0; jagg< numAggElems1; jagg++) {
                        elemid0 = pNode0->getAggElemID(iagg);
                        elemIndex0= moBucket.getIndexElement(elemid0);
                        elemid1 = pNode1->getAggElemID(jagg);
                        elemIndex1= moBucket.getIndexElement(elemid1);
                        if(elemIndex0 == elemIndex1)
                            vnShareElems0.push_back(elemIndex0);
                    };
                };
                numAggElems2= pNode2->getNumOfAggElem();
                vnShareElems1.reserve(numAggElems2);
                for(jagg=0; jagg< vnShareElems0.size(); jagg++) {
                    for(kagg=0; kagg< numAggElems2; kagg++) {
                        elemid2= pNode2->getAggElemID(kagg);
                        elemIndex2= moBucket.getIndexElement(elemid2);
                        if(vnShareElems0[jagg] == elemIndex2)
                            vnShareElems1.push_back(vnShareElems0[jagg]);
                    };
                };
                CElement* pAdjElem;
                if(vnShareElems1.size() > 1) {
                    for(uiint i=0; i< vnShareElems1.size(); i++) {
                        uiint pElemIndex= moBucket.getIndexElement(pElem->getID());
                        if(pElemIndex != vnShareElems1[i]) {
                            pAdjElem= mvElement[vnShareElems1[i]];
                            pElem->setFaceElement(pAdjElem, isurf);
                            pElem->setBoolFaceElem(pNode0,pNode1,pNode2);
                            inNode= GeneInterNode(pNode0);
                            avgCoord(vFaceCnvNode, inNode);
                            inNode->setID(indexCount);
                            pProgMesh->setNode(inNode);
                            pElem->setFaceNode(inNode, isurf);
                            setupParentNode(vFaceCnvNode,inNode);//------------------------------
                            jsurf= pAdjElem->getFaceIndex(pNode0, pNode1, pNode2);
                            pAdjElem->setFaceElement(pElem, jsurf);
                            pAdjElem->setBoolFaceElem(pNode0,pNode1,pNode2);
                            pAdjElem->setFaceNode(inNode, jsurf);
                            indexCount++;
                        }
                    };
                } else {
                    inNode= GeneInterNode(pNode0);
                    avgCoord(vFaceCnvNode, inNode);
                    inNode->setID(indexCount);
                    pProgMesh->setNode(inNode);
                    setupParentNode(vFaceCnvNode,inNode);//------------------------------
                    indexCount++;
                    pElem->setFaceNode(inNode,isurf);
                    pElem->setBoolFaceElem(pNode0, pNode1, pNode2);
                }
            }
        };
    };
}
void CMesh::setupFaceElement2(CMesh* pProgMesh)
{
    CElement *pElem;
    vector<CNode*> vFaceCnvNode;
    uiint  elemIndex0, elemIndex1, elemIndex2, elemid0, elemid1, elemid2;
    CNode *pNode0,*pNode1,*pNode2;
    CNode *inNode;
    uiint ielem,isurf,jsurf, iagg,jagg, kagg;
    uiint numAggElems0, numAggElems1, numAggElems2;
    uiint numOfFace;
    vuint vnShareElems0, vnShareElems1;
    uiint nInitNodeSize;
    uiint indexCount;
    if(pProgMesh) {
        nInitNodeSize = pProgMesh->getNodeSize();
        CNode *pNode= pProgMesh->getNodeIX(nInitNodeSize-1);
        indexCount = pNode->getID()+1;
    } else {
    }
    for(ielem=0; ielem< mNumOfElement; ielem++) {
        pElem= mvElement[ielem];
        numOfFace= pElem->getNumOfFace();
        for(isurf=0; isurf< numOfFace; isurf++) {
            vFaceCnvNode = pElem->getFaceCnvNodes(isurf);
            pNode0= vFaceCnvNode[0];
            pNode1= vFaceCnvNode[1];
            pNode2= vFaceCnvNode[2];
            vnShareElems0.clear();
            vnShareElems1.clear();
            if(!pElem->isFaceElem(pNode0,pNode1,pNode2)) {
                numAggElems0= pNode0->getNumOfAggElem();
                vnShareElems0.reserve(numAggElems0);
                numAggElems1= pNode1->getNumOfAggElem();
                for(iagg=0; iagg< numAggElems0; iagg++) {
                    for(jagg=0; jagg< numAggElems1; jagg++) {
                        elemid0 = pNode0->getAggElemID(iagg);
                        elemIndex0= moBucket.getIndexElement(elemid0);
                        elemid1 = pNode1->getAggElemID(jagg);
                        elemIndex1= moBucket.getIndexElement(elemid1);
                        if(elemIndex0 == elemIndex1)
                            vnShareElems0.push_back(elemIndex0);
                    };
                };
                numAggElems2= pNode2->getNumOfAggElem();
                vnShareElems1.reserve(numAggElems2);
                for(jagg=0; jagg< vnShareElems0.size(); jagg++) {
                    for(kagg=0; kagg< numAggElems2; kagg++) {
                        elemid2= pNode2->getAggElemID(kagg);
                        elemIndex2= moBucket.getIndexElement(elemid2);
                        if(vnShareElems0[jagg] == elemIndex2)
                            vnShareElems1.push_back(vnShareElems0[jagg]);
                    };
                };
                CElement* pAdjElem;
                if(vnShareElems1.size() > 1) {
                    for(uiint i=0; i< vnShareElems1.size(); i++) {
                        uiint pElemIndex= moBucket.getIndexElement(pElem->getID());
                        if(pElemIndex != vnShareElems1[i]) {
                            pAdjElem= mvElement[vnShareElems1[i]];
                            pElem->setFaceElement(pAdjElem, isurf);
                            pElem->setBoolFaceElem(pNode0,pNode1,pNode2);
                            if(pProgMesh) {
                                inNode= GeneInterNode(pNode0);
                                avgCoord(vFaceCnvNode, inNode);
                                inNode->setID(indexCount);
                                pProgMesh->setNode(inNode);
                                pElem->setFaceNode(inNode, isurf);
                                setupParentNode(vFaceCnvNode,inNode);//------------------------------
                            }
                            jsurf= pAdjElem->getFaceIndex(pNode0, pNode1, pNode2);
                            pAdjElem->setFaceElement(pElem, jsurf);
                            pAdjElem->setBoolFaceElem(pNode0,pNode1,pNode2);
                            pAdjElem->setFaceNode(inNode, jsurf);
                            if(pProgMesh) {
                                indexCount++;
                            }
                        }
                    };
                } else {
                    if(pProgMesh) {
                        inNode= GeneInterNode(pNode0);
                        avgCoord(vFaceCnvNode, inNode);
                        inNode->setID(indexCount);
                        pProgMesh->setNode(inNode);
                        setupParentNode(vFaceCnvNode,inNode);//------------------------------
                        indexCount++;
                        pElem->setFaceNode(inNode,isurf);
                        pElem->setBoolFaceElem(pNode0, pNode1, pNode2);
                    }
                }
            }
        };
    };
}
void CMesh::setupVolumeNode(CMesh *pProgMesh)
{
    CElement* pElem;
    CNode *pNode,*cntNode;
    vector<CNode*> vLocalNode;
    uiint nNumOfVert;
    vdouble vCoord;
    uiint ielem, inode;
    uiint nInitNodeSize;
    uiint indexCount;
    nInitNodeSize = pProgMesh->getNodeSize();
    pNode= pProgMesh->getNodeIX(nInitNodeSize-1);
    indexCount = pNode->getID()+1;
    for(ielem=0; ielem< mNumOfElement; ielem++) {
        pElem= mvElement[ielem];
        if(pElem->getEntityType()==BaseElementType::Solid) {
            pNode= pElem->getNode(0);
            cntNode= GeneInterNode(pNode);
            cntNode->setID(indexCount);
            vCoord.clear();
            vCoord.resize(3);
            vCoord[0]=0.0;
            vCoord[1]=0.0;
            vCoord[2]=0.0;
            nNumOfVert= pElem->getNumOfVert();
            vLocalNode.resize(nNumOfVert);
            for(inode=0; inode< nNumOfVert; inode++) {
                pNode= pElem->getNode(inode);
                vLocalNode[inode] = pNode;
                vCoord[0] += pNode->getX();
                vCoord[1] += pNode->getY();
                vCoord[2] += pNode->getZ();
            };
            vCoord[0] /= (double)nNumOfVert;
            vCoord[1] /= (double)nNumOfVert;
            vCoord[2] /= (double)nNumOfVert;
            cntNode->setCoord(vCoord);
            setupParentNode(vLocalNode, cntNode);//-------------------------------------------
            pProgMesh->setNode(cntNode);
            ++indexCount;
            pElem->setVolumeNode(cntNode);
        }
    };
}
void CMesh::avgCoord(vector<CNode*> vCnvNode, CNode *pNode)
{
    vdouble vCoord;
    vCoord.resize(3);
    uiint i;
    for(i=0; i< 3; i++) {
        vCoord[i]=0.0;
    }
    uiint numOfFaceNode = vCnvNode.size();
    for(i=0; i< numOfFaceNode; i++) {
        vCoord[0] += vCnvNode[i]->getX();
        vCoord[1] += vCnvNode[i]->getY();
        vCoord[2] += vCnvNode[i]->getZ();
    };
    for(i=0; i< 3; i++) {
        vCoord[i] /= (double)numOfFaceNode;
    };
    pNode->setCoord(vCoord);
}
void CMesh::setupParentNode(CNode* pNode0, CNode* pNode1, CNode* inNode)//辺ノード
{
    uiint nProgMGLevel= mMGLevel+1; // progMeshのMGLevel

    inNode->reserveParentNode(nProgMGLevel, 2);
    inNode->addParentNode(nProgMGLevel, pNode0);
    inNode->addParentNode(nProgMGLevel, pNode1);
}
void CMesh::setupParentNode2nd(CNode* pNode0, CNode* pNode1, CNode* inNode)//-----最終Level 2次要素 辺ノード
{
    inNode->reserveParentNode(mMGLevel, 2);
    inNode->addParentNode(mMGLevel, pNode0);
    inNode->addParentNode(mMGLevel, pNode1);
}
void CMesh::setupParentNode(vector<CNode*>& vNode, CNode* inNode)//面、体積ノード
{
    uiint nProgMGLevel= mMGLevel+1; // progMeshのMGLevel

    uiint nNumOfParent= vNode.size();
    inNode->reserveParentNode(nProgMGLevel, nNumOfParent);
    uiint ipare;
    for(ipare=0; ipare< nNumOfParent; ipare++) {
        inNode->addParentNode(nProgMGLevel, vNode[ipare]);
    };
}
void CMesh::setCommMesh(CCommMesh* pCommMesh)
{
    mvCommMesh.push_back(pCommMesh);
    uiint comID= pCommMesh->getCommID();
    mmCommIndex[comID]= mvCommMesh.size()-1;
}
CCommMesh* CMesh::getCommMesh(const uiint& comID)
{
    uiint comIndex= mmCommIndex[comID];
    return mvCommMesh[comIndex];
}
void CMesh::sortMesh()
{
    mpLogger->Info(Utility::LoggerMode::MWDebug,"Mesh::sortMesh, initial mvNode.size    => ",(uiint)mvNode.size());
    mpLogger->Info(Utility::LoggerMode::MWDebug,"Mesh::sortMesh, initial mvElement.size => ",(uiint)mvElement.size());
    CNode *pNode, *pDNode;
    CElement *pElem, *pDElem;
    uiint numOfDNode, numOfDElement;
    CCommMesh *pCommMesh;
    uiint numOfComm= mvCommMesh.size();
    uiint icom;
    uiint idel;
    vector<CNode*>::iterator    itNode;
    vector<CElement*>::iterator itElement;
    for(icom=0; icom< numOfComm; icom++) {
        pCommMesh= mvCommMesh[icom];
        numOfDNode= pCommMesh->getNumOfDNode();
        numOfDElement= pCommMesh->getNumOfDCommElement();
        sortID<CNode*>(mvNode, mvNode.size());
        idel=0;
        for(itNode=mvNode.begin(); itNode< mvNode.end(); itNode++) {
            pNode= *itNode;
            if(idel < numOfDNode) {
                pDNode= pCommMesh->getDNode(idel);
                if(pNode->getID() == pDNode->getID()) {
                    mvNode.erase(itNode);
                    idel++;
                    if(itNode != mvNode.begin()) itNode--;
                }
            }
        };
        sortID<CElement*>(mvElement, mvElement.size());
        idel=0;
        for(itElement=mvElement.begin(); itElement< mvElement.end(); itElement++) {
            pElem= *itElement;
            if(idel< numOfDElement) {
                pDElem= pCommMesh->getDElement(idel);
                if(pDElem->getID() == pElem->getID()) {
                    mvElement.erase(itElement);
                    idel++;
                    if(itElement != mvElement.begin()) itElement--;
                }
            }
        };
    };
    mNodeEndIndex= mvNode.size();
    mElemEndIndex= mvElement.size();
    mpLogger->Info(Utility::LoggerMode::MWDebug,"Mesh::sortMesh, mNodeEndIndex => ",mNodeEndIndex);
    mpLogger->Info(Utility::LoggerMode::MWDebug,"Mesh::sortMesh, mElemEndIndex => ",mElemEndIndex);
    uiint idnode, idelem;
    for(icom=0; icom< numOfComm; icom++) {
        pCommMesh= mvCommMesh[icom];
        numOfDNode= pCommMesh->getNumOfDNode();
        numOfDElement= pCommMesh->getNumOfDCommElement();
        for(idnode=0; idnode< numOfDNode; idnode++) {
            pDNode= pCommMesh->getDNode(idnode);
            mvNode.push_back(pDNode);
        };
        for(idelem=0; idelem< numOfDElement; idelem++) {
            pDElem= pCommMesh->getDElement(idelem);
            mvElement.push_back(pDElem);
        };
    };
    mpLogger->Info(Utility::LoggerMode::MWDebug,"Mesh::sortMesh, mvNode.size    => ",(uiint)mvNode.size());
    mpLogger->Info(Utility::LoggerMode::MWDebug,"Mesh::sortMesh, mvElement.size => ",(uiint)mvElement.size());
}
void CMesh::setCommMesh2(CCommMesh2* pCommMesh2)
{
    mvCommMesh2.push_back(pCommMesh2);
    mmComm2Index[pCommMesh2->getID()]= mvCommMesh2.size()-1;
}
CCommMesh2* CMesh::getCommMesh2(const uiint& comID)
{
    uiint index;
    index= mmComm2Index[comID];
    return mvCommMesh2[index];
}
void CMesh::deleteProgData()
{
    uiint ielem;
    for(ielem=0; ielem < mNumOfElement; ielem++)
        mvElement[ielem]->deleteProgData();
}
void CMesh::deleteAggregate_on_Node()
{
    uiint inode;
    for(inode=0; inode < mNumOfNode; inode++)
        mvNode[inode]->deleteAggregate();
}
void CMesh::addElemGrp(CElementGroup* pElemGrp)
{
    mvElementGroup.push_back(pElemGrp);
    uiint nGrpID = pElemGrp->getID();
    mmElemGrpID2IX[nGrpID] = mvElementGroup.size() - 1;
}
uiint CMesh::getNumOfElemGrp()
{
    return mvElementGroup.size();
}
CElementGroup* CMesh::getElemGrpIX(const uiint& index)
{
    return mvElementGroup[index];
}
CElementGroup* CMesh::getElemGrpID(const uiint& nGrpID)
{
    uiint index = mmElemGrpID2IX[nGrpID];
    return mvElementGroup[index];
}


void CMesh::clear()
{
    uiint nNumOfNode= mvNode.size();
    for(uiint i=0; i < nNumOfNode; i++) {
        delete mvNode[i];
    }
    mvNode.clear();
    uiint nNumOfElem= mvElement.size();
    for(uiint i=0; i < nNumOfElem; i++) {
        delete mvElement[i];
    }
    mvElement.clear();


    uiint nNumOfAggElem= mvAggElement.size();
    for(uiint i=0; i < nNumOfAggElem; i++) {
        delete mvAggElement[i];
    }
    mvAggElement.clear();
    uiint nNumOfAggNode= mvAggNode.size();
    for(uiint i=0; i < nNumOfAggNode; i++) {
        delete mvAggNode[i];
    }
    mvAggNode.clear();

    moBucket.clearBucketElement();
    moBucket.clearBucketNode();


    uiint nNumOfCommMesh2= mvCommMesh2.size();
    for(uiint i=0; i < nNumOfCommMesh2; i++) {
        delete mvCommMesh2[i];
    }
    mvCommMesh2.clear();
    mmComm2Index.clear();

    mpBNodeMeshGrp->clear();
    mGrpBndFaceMesh.clear();
    mGrpBndEdgeMesh.clear();
    mGrpBndVolumeMesh.clear();

    if(mvMarkingBNode) delete[] mvMarkingBNode;
    if(mvMarkingDirichlet) delete[] mvMarkingDirichlet;
    if(mvMarkingNeumann) delete[] mvMarkingNeumann;
}

//--
// 境界Node判定 関数 : MG境界処理
//--
// 1.All Boundary
bool CMesh::isBNode(uiint& inode)
{
    return mvMarkingBNode[inode];
}
bool* CMesh::getBNodeMarkingArray()
{
    return mvMarkingBNode;
}
// 2.Dirichlet Boundary
bool CMesh::isDirichletBNode(uiint& inode)
{
    return mvMarkingDirichlet[inode];
}
bool* CMesh::getDirichletBNodeMarkingArray()
{
    return mvMarkingDirichlet;
}
// 3.Neaumann Boundary
bool CMesh::isNeumannBNode(uiint& inode)
{
    return mvMarkingNeumann[inode];
}
bool* CMesh::getNeumannBNodeArray()
{
    return mvMarkingNeumann;
}
//--
// 境界条件Node(All,Dirichlet,Neumann) マーキング配列の生成
//--
void CMesh::setupBNodeMarking()
{
    uiint nNumOfNode= mvNode.size();

    mvMarkingBNode = new bool[nNumOfNode];    // All Boundary
    mvMarkingDirichlet = new bool[nNumOfNode];// Dirichlet Boundary
    mvMarkingNeumann = new bool[nNumOfNode];  // Neumann Boundary

    //Marking Array 初期化
    for(uiint inode=0; inode < nNumOfNode; inode++) {
        mvMarkingBNode[inode]= false;
        mvMarkingDirichlet[inode]= false;
        mvMarkingNeumann[inode]= false;
    };

    //--
    //点Gr
    //--
    uiint nNumOfBNodeMesh = mpBNodeMeshGrp->getNumOfBoundaryNodeMesh();
    for(uiint ibmesh=0; ibmesh < nNumOfBNodeMesh; ibmesh++) {
        CBoundaryNodeMesh *pBNodeMesh = mpBNodeMeshGrp->getBndNodeMeshIX(ibmesh);

        uiint nNumOfBNode = pBNodeMesh->getNumOfBNode();
        //All Boundary
        for(uiint ibnode=0; ibnode < nNumOfBNode; ibnode++) {
            CBoundarySBNode *pBNode = pBNodeMesh->getBNodeIX(ibnode);
            CNode *pNode = pBNode->getNode();
            uiint id = pNode->getID();
            uiint index = moBucket.getIndexNode(id);

            mvMarkingBNode[index]= true;//----------------- BNode マーキング
        };
        //Dirichlet
        if(pBNodeMesh->getBndType()==BoundaryType::Dirichlet) {
            for(uiint ibnode=0; ibnode < nNumOfBNode; ibnode++) {
                CBoundarySBNode *pBNode = pBNodeMesh->getBNodeIX(ibnode);
                CNode *pNode = pBNode->getNode();
                uiint id = pNode->getID();
                uiint index = moBucket.getIndexNode(id);

                mvMarkingDirichlet[index]= true;//----------------- Dirichlet BNode マーキング
            };
        }
        //Neumann
        if(pBNodeMesh->getBndType()==BoundaryType::Neumann) {
            for(uiint ibnode=0; ibnode < nNumOfBNode; ibnode++) {
                CBoundarySBNode *pBNode = pBNodeMesh->getBNodeIX(ibnode);
                CNode *pNode = pBNode->getNode();
                uiint id = pNode->getID();
                uiint index = moBucket.getIndexNode(id);

                mvMarkingNeumann[index]= true;//----------------- Neumann BNode マーキング
            };
        }
    };//点Gr Meshループ

    //--
    //面Gr
    //--
    uiint nNumOfBFaceMesh = mGrpBndFaceMesh.NumOfBoundary();
    for(uiint ibmesh=0; ibmesh < nNumOfBFaceMesh; ibmesh++) {
        CBoundaryFaceMesh *pBFaceMesh = mGrpBndFaceMesh.get_withIndex(ibmesh);

        uiint nNumOfBNode = pBFaceMesh->getNumOfBNode();
        //All Boundary
        for(uiint ibnode=0; ibnode < nNumOfBNode; ibnode++) {
            CBoundaryNode *pBNode = pBFaceMesh->getBNodeIX(ibnode);
            CNode *pNode = pBNode->getNode();
            uiint id = pNode->getID();
            uiint index= moBucket.getIndexNode(id);

            mvMarkingBNode[index]= true;//----------------- BNode マーキング
        };
        //Dirichlet
        if(pBFaceMesh->getBndType()==BoundaryType::Dirichlet) {
            for(uiint ibnode=0; ibnode < nNumOfBNode; ibnode++) {
                CBoundaryNode *pBNode = pBFaceMesh->getBNodeIX(ibnode);
                CNode *pNode = pBNode->getNode();
                uiint id = pNode->getID();
                uiint index= moBucket.getIndexNode(id);

                mvMarkingDirichlet[index]= true;//----------------- Dirichlet BNode マーキング
            };
        }
        //Neumann
        if(pBFaceMesh->getBndType()==BoundaryType::Neumann) {
            for(uiint ibnode=0; ibnode < nNumOfBNode; ibnode++) {
                CBoundaryNode *pBNode = pBFaceMesh->getBNodeIX(ibnode);
                CNode *pNode = pBNode->getNode();
                uiint id = pNode->getID();
                uiint index= moBucket.getIndexNode(id);

                mvMarkingNeumann[index]= true;//----------------- Neumann BNode マーキング
            };
        }
    };//面Gr Meshループ

    //--
    //辺Gr
    //--
    uiint nNumOfBEdgeMesh = mGrpBndEdgeMesh.NumOfBoundary();
    for(uiint ibmesh=0; ibmesh < nNumOfBEdgeMesh; ibmesh++) {
        CBoundaryEdgeMesh *pBEdgeMesh = mGrpBndEdgeMesh.get_withIndex(ibmesh);

        uiint nNumOfBNode = pBEdgeMesh->getNumOfBNode();
        //All Boundary
        for(uiint ibnode=0; ibnode < nNumOfBNode; ibnode++) {
            CBoundaryNode *pBNode = pBEdgeMesh->getBNodeIX(ibnode);
            CNode *pNode = pBNode->getNode();
            uiint id = pNode->getID();
            uiint index= moBucket.getIndexNode(id);

            mvMarkingBNode[index]= true;//----------------- BNode マーキング
        }
        //Dirichlet
        if(pBEdgeMesh->getBndType()==BoundaryType::Dirichlet) {
            for(uiint ibnode=0; ibnode < nNumOfBNode; ibnode++) {
                CBoundaryNode *pBNode = pBEdgeMesh->getBNodeIX(ibnode);
                CNode *pNode = pBNode->getNode();
                uiint id = pNode->getID();
                uiint index= moBucket.getIndexNode(id);

                mvMarkingDirichlet[index]= true;//----------------- Dirichlet BNode マーキング
            }
        }
        //Neumann
        if(pBEdgeMesh->getBndType()==BoundaryType::Neumann) {
            for(uiint ibnode=0; ibnode < nNumOfBNode; ibnode++) {
                CBoundaryNode *pBNode = pBEdgeMesh->getBNodeIX(ibnode);
                CNode *pNode = pBNode->getNode();
                uiint id = pNode->getID();
                uiint index= moBucket.getIndexNode(id);

                mvMarkingNeumann[index]= true;//----------------- Neumann BNode マーキング
            }
        }

    };//辺Gr Meshループ

    //--
    //体積Gr
    //--
    uiint nNumOfBVolMesh = mGrpBndVolumeMesh.NumOfBoundary();
    for(uiint ibmesh=0; ibmesh < nNumOfBVolMesh; ibmesh++) {
        CBoundaryVolumeMesh *pBVolMesh = mGrpBndVolumeMesh.get_withIndex(ibmesh);

        uiint nNumOfBNode = pBVolMesh->getNumOfBNode();
        //All Boundary
        for(uiint ibnode=0; ibnode < nNumOfBNode; ibnode++) {
            CBoundaryNode *pBNode = pBVolMesh->getBNodeIX(ibnode);
            CNode *pNode = pBNode->getNode();
            uiint id= pNode->getID();
            uiint index= moBucket.getIndexNode(id);

            mvMarkingBNode[index]= true;//----------------- BNode マーキング
        }
        //Dirichlet
        if(pBVolMesh->getBndType()==BoundaryType::Dirichlet) {
            for(uiint ibnode=0; ibnode < nNumOfBNode; ibnode++) {
                CBoundaryNode *pBNode = pBVolMesh->getBNodeIX(ibnode);
                CNode *pNode = pBNode->getNode();
                uiint id= pNode->getID();
                uiint index= moBucket.getIndexNode(id);

                mvMarkingDirichlet[index]= true;//----------------- Dirichlet BNode マーキング
            }
        }
        //Neumann
        if(pBVolMesh->getBndType()==BoundaryType::Neumann) {
            for(uiint ibnode=0; ibnode < nNumOfBNode; ibnode++) {
                CBoundaryNode *pBNode = pBVolMesh->getBNodeIX(ibnode);
                CNode *pNode = pBNode->getNode();
                uiint id= pNode->getID();
                uiint index= moBucket.getIndexNode(id);

                mvMarkingNeumann[index]= true;//----------------- Neumann BNode マーキング
            }
        }

    };//体積Gr Meshループ
}
//--
// Node内のRank大の通信ノードをマーキング
//--
void CMesh::setupLargeRank_CommNodeMarking()
{
    // マーキング初期化
    uiint nNumOfNode= mvNode.size();
    mvMarkingLargeRankNode= new bool[nNumOfNode];

    for(uiint inode=0; inode < nNumOfNode; inode++) {
        mvMarkingLargeRankNode[inode]= false;
    };

    // 通信メッシュからデータ取得
    uiint nNumOfCommMesh = mvCommMesh2.size();

    for(uiint icom=0; icom < nNumOfCommMesh; icom++) {
        CCommMesh2 *pCommMesh2= mvCommMesh2[icom];

        uiint myRank= pCommMesh2->getRank();
        uiint transRank= pCommMesh2->getTrasmitRank();

        // 相手よりも自身Rankが大の場合：マーキング
        // # ひとつの通信ペアでもRankが大であれば、Rank大マーキング
        // # 最小Rank以外は全てマーキングされる.
        if(myRank > transRank) {
            uiint nNumOfCNode= pCommMesh2->getCommNodeSize();
            for(uiint icnode=0; icnode < nNumOfCNode; icnode++) {
                CCommNode *pCommNode= pCommMesh2->getCommNodeIX(icnode);
                CNode *pNode= pCommNode->getNode();
                uiint index = moBucket.getIndexNode(pNode->getID());

                mvMarkingLargeRankNode[index]= true;
            };
        }//if(myRank大)
    };
}
bool CMesh::isLargeRankCommNode(const uiint& inode)
{
    return mvMarkingLargeRankNode[inode];
}
////void CMesh::addPairRank()
////{
////    uiint nNumOfCommMesh = mvCommMesh2.size();
////
////    for(uiint icom=0; icom < nNumOfCommMesh; icom++){
////        CCommMesh2 *pCommMesh2= mvCommMesh2[icom];
////
////        uiint nNumOfCNode= pCommMesh2->getCommNodeSize();
////        for(uiint icnode=0; icnode < nNumOfCNode; icnode++){
////            CCommNode *pCommNode= pCommMesh2->getCommNodeIX(icnode);
////            CNode *pNode= pCommNode->getNode();
////
////            uiint myRank= pCommMesh2->getRank();
////            uiint transRank= pCommMesh2->getTrasmitRank();
////
////            pNode->addRank(myRank, transRank);//----------------- rankペア
////        };
////    };
////}




