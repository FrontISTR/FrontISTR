/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Mesh.cxx
|
|                     Written by T.Takeda,    2010/06/01
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
#include "CommunicationMesh.h"
#include "Logger.h"
#include "Node.h"
#include <vector>
#include "VectorNode.h"
#include "BoundaryGroup.h"
#include "Element.h"
#include "AggregateElement.h"
#include "Mesh.h"
#include "ScalarNode.h"
#include "ScalarVectorNode.h"
using namespace pmw;
CMesh::CMesh(void)
{
    mpLogger = Utility::CLogger::Instance();
    mnDummyCount= 0;
}
CMesh::CMesh(const uint& numofnode, const uint& numofelem)
{
    mpLogger = Utility::CLogger::Instance();
    mNumOfNode = numofnode;  mNumOfElement = numofelem;
    mvNode.resize(mNumOfNode);
    mvElement.resize(mNumOfElement);
}
CMesh::~CMesh(void)
{
    cout << "~CMesh start,  mMGLevel==" << mMGLevel << endl;
    uint i;
    CNode* pNode;
    for(i=0; i< mvNode.size(); i++){
        pNode= mvNode[i];
        if(pNode){
            if(pNode->getMGLevel()==mMGLevel) delete pNode;
        }
        mvNode.erase(mvNode.begin()+i);
    };
    CElement *pElem;
    for(i=0; i< mvElement.size(); i++){
        pElem= mvElement[i];
        if(pElem){
            if(pElem->getMGLevel()==mMGLevel) delete pElem;
        }
        mvElement.erase(mvElement.begin()+i);
    };
    for_each(mvAggElement.begin(), mvAggElement.end(), DeleteObject());
    for_each(mvAggNode.begin(), mvAggNode.end(), DeleteObject());
    for_each(mvCommMesh.begin(), mvCommMesh.end(), DeleteObject());
    CBoundaryFace *pBoundFace;
    for(i=0; i< mBoundaryFaces.NumOfBoundary(); i++){
        pBoundFace = mBoundaryFaces.get_withIndex(i);
        delete pBoundFace;
    };
    CBoundaryVolume *pBoundVol;
    for(i=0; i < mBoundaryVolumes.NumOfBoundary(); i++){
        pBoundVol = mBoundaryVolumes.get_withIndex(i);
        delete pBoundVol;
    };
    CBoundaryNode *pBoundNode;
    for(i=0; i < mBoundaryNodes.NumOfBoundary(); i++){
        pBoundNode = mBoundaryNodes.get_withIndex(i);
        delete pBoundNode;
    };
    cout << "~CMesh   end,  mMGLevel==" << mMGLevel << endl;
}
void CMesh::initBucketNode(const uint& max_id, const uint& min_id)
{
    moBucket.resizeBucketNode(max_id, min_id);
}
void CMesh::setupBucketNodeIndex(const uint& id, const uint& index)
{
    moBucket.setIndexNode(id, index);
}
void CMesh::setupBucketNode()
{
    CNode *pNode;
    uint maxID,minID, i;
    pNode = mvNode[0];
    maxID = pNode->getID(); minID = pNode->getID();
    for(i=0; i < mvNode.size(); i++){
        pNode = mvNode[i];
        if(pNode->getID() > maxID) maxID = pNode->getID();
        if(pNode->getID() < minID) minID = pNode->getID();
    };
    moBucket.resizeBucketNode(maxID, minID);
    for(i=0; i < mvNode.size(); i++){
        pNode = mvNode[i];
        moBucket.setIndexNode(pNode->getID(), i);
    };
}
void CMesh::initBucketElement(const uint& max_id, const uint& min_id)
{
    moBucket.resizeBucketElement(max_id, min_id);
}
void CMesh::setupBucketElementIndex(const uint& id, const uint& index)
{
    moBucket.setIndexElement(id, index);
}
void CMesh::setupBucketElement()
{
    CElement *pElement;
    uint maxID,minID, i;
    pElement = mvElement[0];
    maxID = pElement->getID(); minID = pElement->getID();
    for(i=0; i < mvElement.size(); i++){
        pElement = mvElement[i];
        if(pElement->getID() > maxID) maxID = pElement->getID();
        if(pElement->getID() < minID) minID = pElement->getID();
    };
    moBucket.resizeBucketElement(maxID,minID);
    for(i=0; i < mvElement.size(); i++){
        pElement = mvElement[i];
        moBucket.setIndexElement(pElement->getID(), i);
    };
}
void CMesh::reserveNode(const uint& num_of_node)
{
    mvNode.reserve(num_of_node);
}
void CMesh::setNode(CNode *pNode)
{
    mvNode.push_back(pNode);
}
CNode* CMesh::getNode(const uint &id)
{
    uint index;
    index = moBucket.getIndexNode(id);
    return mvNode[index];
}
CNode* CMesh::getNodeIX(const uint& index)
{
    return mvNode[index];
}
void CMesh::reserveElement(const uint& num_of_elem)
{
    mvElement.reserve(num_of_elem);
}
void CMesh::setElement(CElement *pElement)
{
    mvElement.push_back(pElement);
}
CElement* CMesh::getElement(const uint& id)
{
    uint index;
    index = moBucket.getIndexElement(id);
    return mvElement[index];
}
CElement* CMesh::getElementIX(const uint& index)
{
    return mvElement[index];
}
void CMesh::reserveAggregate(const uint& res_size)
{
    mvAggElement.reserve(res_size);
    mvAggNode.reserve(res_size);
}
void CMesh::setAggElement(CAggregateElement* pAggElem)
{
    mvAggElement.push_back(pAggElem);
}
void CMesh::setAggNode(CAggregateNode* pAggNode)
{
    mvAggNode.push_back(pAggNode);
}
CAggregateElement* CMesh::getAggElem(const uint& node_id)
{
    return mvAggElement[node_id];
}
CAggregateNode* CMesh::getAggNode(const uint& node_id)
{
    return mvAggNode[node_id];
}
void CMesh::setupAggregate()
{
    CElement *pElem;   CNode *pNode;
    uint ielem, local_id, inode;
    mNumOfNode= mvNode.size();
    mNumOfElement= mvElement.size();
    for(inode=0; inode< mNumOfNode; inode++){
        pNode= mvNode[inode];
        pNode->clearAggElemID();
        pNode->clearNeibElemVert();
    };
    uint numOfLocalNode; uint elemID;
    for(ielem=0; ielem< mNumOfElement; ielem++){
        pElem = mvElement[ielem];
        numOfLocalNode= pElem->getNumOfNode();
        elemID = pElem->getID();
        for(local_id=0; local_id< numOfLocalNode;local_id++){
            pNode = pElem->getNode(local_id);
            pNode->setAggElemID(elemID);      
            pNode->setNeibElemVert(elemID,local_id);
        };
    };
    CAggregateElement *pAggElem;
    uint elemIndex, iagg_elem, numOfAggElem;
    for(inode=0; inode< mNumOfNode; inode++){
        pNode = mvNode[inode];
        numOfAggElem = pNode->getNumOfAggElem();
        pAggElem = mvAggElement[pNode->getID()];
        pAggElem->reserve(numOfAggElem);
        for(iagg_elem=0; iagg_elem< numOfAggElem; iagg_elem++){
            elemID = pNode->getAggElemID(iagg_elem);
            pElem= getElement(elemID);
            pAggElem->push(pElem);
        };
    };
    CAggregateNode *pAggNode;
    vector<CNode*> vConnNode;
    uint ncon, nNumOfConnNode;
    for(inode=0; inode< mNumOfNode; inode++){
        pNode= mvNode[inode];
        numOfAggElem= pNode->getNumOfAggElem();
        pAggNode = mvAggNode[pNode->getID()];
        pAggNode->reserveNode(numOfAggElem);
        for(iagg_elem=0; iagg_elem< numOfAggElem; iagg_elem++){
            elemID= pNode->getAggElemID(iagg_elem);
            pElem= getElement(elemID);
            vConnNode= pElem->getConnectNode(pNode);
            nNumOfConnNode= vConnNode.size();
            for(ncon=0; ncon < nNumOfConnNode; ncon++){
                pAggNode->setNode(vConnNode[ncon]);
            };
        };
    };
}
void CMesh::presetProgMesh(CMesh* pProgMesh)
{
    uint progNumOfNode= mNumOfNode;
    uint progNumOfElem= mNumOfElement;
    CElement *pElem;
    uint ielem;
    for(ielem=0; ielem< mNumOfElement; ielem++){
        pElem= mvElement[ielem];
        switch(pElem->getType()){
            case(ElementType::Hexa):case(ElementType::Hexa2):
                progNumOfElem += 7;
                break;
            case(ElementType::Tetra):case(ElementType::Tetra2):
                progNumOfElem += 3;
                break;
            case(ElementType::Prism):case(ElementType::Prism2):
                progNumOfElem += 5;
                break;
            case(ElementType::Pyramid):case(ElementType::Pyramid2):
                progNumOfElem += 7;
                break;
            case(ElementType::Quad):case(ElementType::Quad2):
                progNumOfElem += 3;
                break;
            case(ElementType::Triangle):case(ElementType::Triangle2):
                progNumOfElem += 2;
                break;
            case(ElementType::Beam):case(ElementType::Beam2):
                progNumOfElem += 1;
                break;
        }
    };
    progNumOfNode *= 8;
    pProgMesh->reserveNode(progNumOfNode);   
    pProgMesh->reserveElement(progNumOfElem);
    uint i;
    for(i=0; i< mNumOfNode; i++){
        pProgMesh->setNode(mvNode[i]);
    };
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::MWDebug,"prolongation numOfElem",progNumOfElem);
    pLogger->Info(Utility::LoggerMode::MWDebug,"prolongation numOfNode",progNumOfNode);
}
void CMesh::setupEdgeElement(CMesh *pProgMesh)
{
    CElement *pElem;
    PairNode edgeNode;
    CNode *pNode0,*pNode1;
    CElement *pOtherElem;
    vector<CElement*> vEdgeElem;
    uint nOtherEdge;
    uint numOfEdge;
    uint ielem,iedge;
    uint indexCount;
    CNode *inNode;  
    if(pProgMesh){
        indexCount = pProgMesh->getNodeSize();
    }else{
        indexCount= mvNode.size();
    }
    for(ielem=0; ielem< mNumOfElement; ielem++){
        pElem = mvElement[ielem];
        numOfEdge = pElem->getNumOfEdge();
        for(iedge=0; iedge< numOfEdge; iedge++){
            edgeNode = pElem->getPairNode(iedge);
            pNode0= edgeNode.first;  pNode1= edgeNode.second;
            if(!pElem->isEdgeElem(pNode0,pNode1)){
                uint iagg,jagg, elemIndex0, elemIndex1, elemid0, elemid1;
                uint numOfAggElem0=pNode0->getNumOfAggElem();
                uint numOfAggElem1=pNode1->getNumOfAggElem();
                pElem->reserveEdgeElement(iedge,numOfAggElem0);
                for(iagg=0; iagg< numOfAggElem0; iagg++){
                    elemid0= pNode0->getAggElemID(iagg);
                    elemIndex0= moBucket.getIndexElement(elemid0);
                    for(jagg=0; jagg< numOfAggElem1; jagg++){
                        elemid1= pNode1->getAggElemID(jagg);
                        elemIndex1= moBucket.getIndexElement(elemid1);
                        if(elemIndex0 == elemIndex1)
                            pElem->setEdgeElement(iedge, mvElement[elemIndex0]);
                    };
                };
                inNode= GeneInterNode(pNode0);
                vdouble P0_coord= pNode0->getCoord();
                vdouble P1_coord= pNode1->getCoord();
                vdouble In_coord; In_coord.reserve(3);
                for(uint i=0; i< 3; i++){ In_coord.push_back( (P0_coord[i] + P1_coord[i])*0.5 );}
                inNode->setCoord(In_coord);
                inNode->setID(indexCount); 
                setupParentNode(pNode0,pNode1,inNode);
                setupChildNode(pNode0,pNode1, inNode);
                pElem->setEdgeInterNode(inNode,iedge);
                pElem->setBoolEdgeElem(pNode0, pNode1);
                if(pProgMesh){
                    pProgMesh->setNode(inNode);
                    setEdgeNode(inNode); 
                }else{
                    setEdgeNode(inNode);
                }
                indexCount++;
                vEdgeElem.clear();
                vEdgeElem= pElem->getEdgeElement(iedge);
                uint n;
                for(n=0; n< vEdgeElem.size(); n++){
                    pOtherElem=vEdgeElem[n];
                    if(pOtherElem->getID() != pElem->getID() && !pOtherElem->isEdgeElem(pNode0,pNode1)){
                        nOtherEdge =pOtherElem->getEdgeIndex(pNode0, pNode1);
                        pOtherElem->setEdgeAggElement(nOtherEdge,vEdgeElem);
                        pOtherElem->setBoolEdgeElem(pNode0, pNode1);        
                        pOtherElem->setEdgeInterNode(inNode,nOtherEdge);
                    }
                };
            }
        };
    };
}
CNode* CMesh::GeneInterNode(CNode* pNode)
{
    CNode *inNode;
    uint numOfScalar, numOfVector;
    switch(pNode->getType()){
        case(NodeType::Scalar):
            inNode= new CScalarNode;
            numOfScalar= pNode->numOfScalarParam();
            inNode->resizeScalar(numOfScalar);
            break;
        case(NodeType::Vector):
            inNode= new CVectorNode;
            numOfVector= pNode->numOfVectorParam();
            inNode->resizeVector(numOfVector);
            break;
        case(NodeType::ScalarVector):
            inNode= new CScalarVectorNode;
            numOfScalar= pNode->numOfScalarParam();
            numOfVector= pNode->numOfVectorParam();
            inNode->resizeScalar(numOfScalar);
            inNode->resizeVector(numOfVector);
            break;
        default:
            mpLogger->Info(Utility::LoggerMode::Error,"Node Generation Error, CMesh::GeneInterNode");
            break;
    }
    inNode->setMGLevel(mMGLevel+1);
    return inNode;
}
void CMesh::setupFaceElement(CMesh* pProgMesh)
{
    CElement *pElem;
    vector<CNode*> vFaceCnvNode;
    uint  elemIndex0, elemIndex1, elemIndex2, elemid0, elemid1, elemid2;
    CNode *pNode0,*pNode1,*pNode2;
    CNode *inNode;
    uint ielem,isurf,jsurf, iagg,jagg, kagg;
    uint numAggElems0, numAggElems1, numAggElems2;
    uint numOfFace;
    vuint vnShareElems0, vnShareElems1;
    uint indexCount;
    indexCount = pProgMesh->getNodeSize();
    for(ielem=0; ielem< mNumOfElement; ielem++){
        pElem= mvElement[ielem];
        numOfFace= pElem->getNumOfFace();
        for(isurf=0; isurf< numOfFace; isurf++){
            vFaceCnvNode = pElem->getFaceCnvNodes(isurf);
            pNode0= vFaceCnvNode[0]; pNode1= vFaceCnvNode[1]; pNode2= vFaceCnvNode[2];
            vnShareElems0.clear(); vnShareElems1.clear();
            if(!pElem->isFaceElem(pNode0,pNode1,pNode2)){
                numAggElems0= pNode0->getNumOfAggElem(); vnShareElems0.reserve(numAggElems0);
                numAggElems1= pNode1->getNumOfAggElem();
                for(iagg=0; iagg< numAggElems0; iagg++){
                for(jagg=0; jagg< numAggElems1; jagg++){
                    elemid0 = pNode0->getAggElemID(iagg);
                    elemIndex0= moBucket.getIndexElement(elemid0);
                    elemid1 = pNode1->getAggElemID(jagg);
                    elemIndex1= moBucket.getIndexElement(elemid1);
                    if(elemIndex0 == elemIndex1)
                        vnShareElems0.push_back(elemIndex0);
                };
                };
                numAggElems2= pNode2->getNumOfAggElem(); vnShareElems1.reserve(numAggElems2);
                for(jagg=0; jagg< vnShareElems0.size(); jagg++){
                for(kagg=0; kagg< numAggElems2; kagg++){
                    elemid2= pNode2->getAggElemID(kagg);
                    elemIndex2= moBucket.getIndexElement(elemid2);
                    if(vnShareElems0[jagg] == elemIndex2)
                        vnShareElems1.push_back(vnShareElems0[jagg]);
                };
                };
                CElement* pAdjElem;
                if(vnShareElems1.size()> 1){
                    for(uint i=0; i< vnShareElems1.size(); i++){
                        uint pElemIndex= moBucket.getIndexElement(pElem->getID());
                        if(pElemIndex != vnShareElems1[i]){
                            pAdjElem= mvElement[vnShareElems1[i]];
                            pElem->setFaceElement(pAdjElem, isurf);
                            pElem->setBoolFaceElem(pNode0,pNode1,pNode2);
                            inNode= GeneInterNode(pNode0);
                            avgCoord(vFaceCnvNode, inNode);
                            inNode->setID(indexCount);
                            pProgMesh->setNode(inNode);
                            pElem->setFaceNode(inNode, isurf);
                            setupParentNode(vFaceCnvNode,inNode);
                            setupChildNode(vFaceCnvNode, inNode);
                            jsurf= pAdjElem->getFaceIndex(pNode0, pNode1, pNode2);
                            pAdjElem->setFaceElement(pElem, jsurf);
                            pAdjElem->setBoolFaceElem(pNode0,pNode1,pNode2);
                            pAdjElem->setFaceNode(inNode, jsurf);           
                            indexCount++;
                        }
                    };
                }else{
                    inNode= GeneInterNode(pNode0);
                    avgCoord(vFaceCnvNode, inNode);
                    inNode->setID(indexCount);
                    pProgMesh->setNode(inNode);
                    setupParentNode(vFaceCnvNode,inNode);
                    setupChildNode(vFaceCnvNode, inNode);
                    indexCount++;
                    pElem->setFaceNode(inNode,isurf);
                    pElem->setBoolFaceElem(pNode0, pNode1, pNode2);
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
    uint numOfLocalNode;
    vdouble vCoord;
    uint ielem, inode;
    uint indexCount;
    indexCount = pProgMesh->getNodeSize();
    for(ielem=0; ielem< mNumOfElement; ielem++){
        pElem= mvElement[ielem];
        if(pElem->getEntityType()==BaseElementType::Solid){
            pNode= pElem->getNode(0);
            cntNode= GeneInterNode(pNode);
            cntNode->setID(indexCount);
            vCoord.clear();
            vCoord.resize(3);
            vCoord[0]=0.0; vCoord[1]=0.0; vCoord[2]=0.0;
            vLocalNode = pElem->getNode();
            numOfLocalNode= pElem->getNumOfNode();
            for(inode=0; inode< numOfLocalNode; inode++){
                pNode= vLocalNode[inode];
                vCoord[0] += pNode->getX();
                vCoord[1] += pNode->getY();
                vCoord[2] += pNode->getZ();
            };
            vCoord[0] /= (double)numOfLocalNode;  vCoord[1] /= (double)numOfLocalNode;  vCoord[2] /= (double)numOfLocalNode;
            cntNode->setCoord(vCoord);
            setupParentNode(vLocalNode, cntNode);
            setupChildNode(vLocalNode, cntNode); 
            pProgMesh->setNode(cntNode);
            ++indexCount;
            pElem->setVolumeNode(cntNode);
        }
    };
}
void CMesh::avgCoord(vector<CNode*> vCnvNode, CNode *pNode)
{
    vdouble vCoord; vCoord.resize(3);
    uint i;
    for(i=0; i< 3; i++){ vCoord[i]=0.0;}
    uint numOfFaceNode = vCnvNode.size();
    for(i=0; i< numOfFaceNode; i++){
        vCoord[0] += vCnvNode[i]->getX();
        vCoord[1] += vCnvNode[i]->getY();
        vCoord[2] += vCnvNode[i]->getZ();
    };
    for(i=0; i< 3; i++){
        vCoord[i] /= (double)numOfFaceNode;
    };
    pNode->setCoord(vCoord);
}
void CMesh::setupParentNode(CNode* pNode0, CNode* pNode1, CNode* inNode)
{
    inNode->reserveParentNode(2);
    inNode->addParentNode(pNode0);
    inNode->addParentNode(pNode1);
}
void CMesh::setupParentNode(vector<CNode*>& vNode, CNode* inNode)
{
    uint numOfParent= vNode.size();
    inNode->reserveParentNode(numOfParent);
    uint ipare;
    for(ipare=0; ipare< numOfParent; ipare++){
        inNode->addParentNode(vNode[ipare]);
    };
}
void CMesh::setupChildNode(CNode* pNode0, CNode* pNode1, CNode* inNode)
{
    pNode0->addChildNode(inNode);
    pNode1->addChildNode(inNode);
}
void CMesh::setupChildNode(vector<CNode*>& vNode, CNode* inNode)
{
    uint numOfParent= vNode.size();
    uint ipare;
    for(ipare=0; ipare < numOfParent; ipare++){
        vNode[ipare]->addChildNode(inNode);
    };
}
void CMesh::setCommMesh(CCommMesh* pCommMesh)
{
    mvCommMesh.push_back(pCommMesh);
    uint comID= pCommMesh->getCommID();
    mmCommIndex[comID]= mvCommMesh.size()-1;
}
CCommMesh* CMesh::getCommMesh(const uint& comID)
{
    uint comIndex= mmCommIndex[comID];
    return mvCommMesh[comIndex];
}
void CMesh::sortMesh()
{
    mpLogger->Info(Utility::LoggerMode::MWDebug,"Mesh::sortMesh, initial mvNode.size    => ",(uint)mvNode.size());
    mpLogger->Info(Utility::LoggerMode::MWDebug,"Mesh::sortMesh, initial mvElement.size => ",(uint)mvElement.size());
    CNode *pNode, *pDNode;
    CElement *pElem, *pDElem;
    uint numOfDNode, numOfDElement;
    CCommMesh *pCommMesh;
    uint numOfComm= mvCommMesh.size();
    uint icom;
    uint idel;
    vector<CNode*>::iterator    itNode;
    vector<CElement*>::iterator itElement;
    for(icom=0; icom< numOfComm; icom++){
        pCommMesh= mvCommMesh[icom];
        numOfDNode= pCommMesh->getNumOfDNode();
        numOfDElement= pCommMesh->getNumOfDCommElement();
        sortID<CNode*>(mvNode, mvNode.size());
        idel=0;
        for(itNode=mvNode.begin(); itNode< mvNode.end(); itNode++){
            pNode= *itNode;
            if(idel < numOfDNode){
                pDNode= pCommMesh->getDNode(idel);
                if(pNode->getID() == pDNode->getID()){
                    mvNode.erase(itNode);
                    idel++; 
                    if(itNode != mvNode.begin()) itNode--;
                }
            }
        };
        sortID<CElement*>(mvElement, mvElement.size());
        idel=0;
        for(itElement=mvElement.begin(); itElement< mvElement.end(); itElement++){
            pElem= *itElement;
            if(idel< numOfDElement){
                pDElem= pCommMesh->getDElement(idel);
                if(pDElem->getID() == pElem->getID()){
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
    uint idnode, idelem;
    for(icom=0; icom< numOfComm; icom++){
        pCommMesh= mvCommMesh[icom];
        numOfDNode= pCommMesh->getNumOfDNode();
        numOfDElement= pCommMesh->getNumOfDCommElement();
        for(idnode=0; idnode< numOfDNode; idnode++){
            pDNode= pCommMesh->getDNode(idnode);
            mvNode.push_back(pDNode);
        };
        for(idelem=0; idelem< numOfDElement; idelem++){
            pDElem= pCommMesh->getDElement(idelem);
            mvElement.push_back(pDElem);
        };
    };
    mpLogger->Info(Utility::LoggerMode::MWDebug,"Mesh::sortMesh, mvNode.size    => ",(uint)mvNode.size());
    mpLogger->Info(Utility::LoggerMode::MWDebug,"Mesh::sortMesh, mvElement.size => ",(uint)mvElement.size());
}
void CMesh::setCommMesh2(CCommMesh2* pCommMesh2)
{
    mvCommMesh2.push_back(pCommMesh2);
    mmComm2Index[pCommMesh2->getID()]= mvCommMesh2.size()-1;
}
CCommMesh2* CMesh::getCommMesh2(const uint& comID)
{
    uint index;
    index= mmComm2Index[comID];
    return mvCommMesh2[index];
}
void CMesh::setEdgeNode(CNode* pNode)
{
    uint id,index;
    id= pNode->getID();
    index= mvEdgeNode.size();
    mvEdgeNode.push_back(pNode);
    mmEdgeNodeID2IX[id]= index;
}
CNode* CMesh::getEdgeNode(const uint& id)
{
    uint index;
    index= mmEdgeNodeID2IX[id];
    return mvEdgeNode[index];
}
