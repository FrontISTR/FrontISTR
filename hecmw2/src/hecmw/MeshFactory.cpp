/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/MeshFactory.cpp
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "HEC_MPI.h"
#include "ElementType.h"
#include "BoundaryMesh.h"
#include "BoundaryParts.h"
#include "BoundaryFace.h"
#include "BndVertex.h"
#include "BoundaryNodeMesh.h"
#include "BoundaryFaceMesh.h"
#include "CommNode.h"
#include "CommMesh2.h"
#include "BoundaryNode.h"
#include <vector>
#include "FaceTree.h"
#include "EdgeTree.h"
#include "ContactNode.h"
#include "GMGModel.h"
#include "CommFace.h"
#include "ContactMesh.h"
#include "AssyModel.h"
#include "Element.h"
#include "Mesh.h"
#include "MeshFactory.h"
#include "BoundaryHexa.h"
#include "SkinFace.h"
using namespace pmw;
CMeshFactory::CMeshFactory(void)
{
    mpLogger = Utility::CLogger::Instance();
}
CMeshFactory::~CMeshFactory(void)
{
}
void CMeshFactory::refineMesh()
{
    if(mMGLevel > 0){
        MGMeshConstruct();
    }else{
        SGMeshConstruct();
    }
}
void CMeshFactory::SGMeshConstruct()
{
    CAssyModel *pAssy;
    CMesh      *pMesh;
    pAssy = mpGMGModel->getAssyModel(0);
    uiint nLevel=0;
    uiint numOfMesh = pAssy->getNumOfMesh();
    uiint imesh;
    for(imesh=0; imesh < numOfMesh; imesh++){
        pMesh = pAssy->getMesh(imesh);
        pMesh->setupAggregate(nLevel);
        string str = "SGMeshConstruct , pMesh->setupAggElement  at " + boost::lexical_cast<string>(imesh);
        mpLogger->Info(Utility::LoggerMode::MWDebug, str);
    };
}
void CMeshFactory::MGMeshConstruct()
{
    CAssyModel *pAssy, *pProgAssy;
    CMesh      *pMesh, *pProgMesh;
    CElement         *pElem=NULL;
    vector<CElement*> vProgElem; 
    CCommMesh  *pCommMesh, *pProgCommMesh;
    CCommElement      *pCommElem;
    vector<CCommElement*> vProgCommElem;
    uiint numOfCommMesh,numOfCommElemAll,numOfCommNode;
    uiint icommesh,icomelem,iprocom;
    pAssy= mpGMGModel->getAssyModel(0);
    uiint numOfMesh= pAssy->getNumOfMesh();
    uiint ilevel,imesh,ielem;
    for(ilevel=0; ilevel< mMGLevel; ilevel++){
        pAssy= mpGMGModel->getAssyModel(ilevel);
        pProgAssy= mpGMGModel->getAssyModel(ilevel+1);
        pProgAssy->resizeMesh(numOfMesh);
        pProgAssy->intializeBucket(pAssy->getMaxMeshID(),pAssy->getMinMeshID());
        pProgAssy->setMaxMeshID(pAssy->getMaxMeshID());
        pProgAssy->setMinMeshID(pAssy->getMinMeshID());
        for(imesh=0; imesh< numOfMesh; imesh++){
            pMesh= pAssy->getMesh(imesh);
            pProgAssy->setBucket(pMesh->getMeshID(), imesh);
            pProgMesh = new CMesh;          
            pProgMesh->setMGLevel(ilevel+1);
            pProgMesh->setMeshID(pMesh->getMeshID());
            if(ilevel==0){
                pMesh->setupAggregate(ilevel); 
                mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->setupAggElement finish at ilevel==0");
            }
            pMesh->presetProgMesh(pProgMesh);
            mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->presetProgMesh finish");
            pMesh->setupEdgeElement(pProgMesh, ilevel);
            mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->setupEdgeElement finish");
            pMesh->setupFaceElement(pProgMesh);
            mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->setupFaceElement finish");
            pMesh->setupVolumeNode(pProgMesh); 
            mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->setupVolumeNode finish");
            pMesh->replaceEdgeNode();
            uiint numOfElem = pMesh->getNumOfElement();
            uiint elementID= 0;
            for(ielem=0; ielem< numOfElem; ielem++){
                pElem= pMesh->getElementIX(ielem);
                vProgElem.clear();
                GeneProgElem(ilevel, pElem, vProgElem, elementID, pProgMesh);
            };
            pProgMesh->setupNumOfNode();
            pProgMesh->setupNumOfElement();
            pProgMesh->setSolutionType(mnSolutionType);
            uiint numOfNode= pProgMesh->getNodeSize();
            CAggregateElement *pAggElem;
            CAggregateNode    *pAggNode;
            pProgMesh->resizeAggregate(numOfNode);
            if(mnSolutionType==SolutionType::FEM){
                for(uiint iAgg=0; iAgg< numOfNode; iAgg++){
                    pAggElem = new CAggregateElement;
                    pProgMesh->setAggElement(pAggElem, iAgg);
                };
            }
            if(mnSolutionType==SolutionType::FVM){
                for(uiint iAgg=0; iAgg< numOfNode; iAgg++){
                    pAggNode = new CAggregateNode;
                    pProgMesh->setAggNode(pAggNode, iAgg);
                };
            }
            pProgAssy->setMesh(pProgMesh,pMesh->getMeshID());
            uiint progNodeSize = pProgMesh->getNodeSize();
            CNode *pNode = pProgMesh->getNodeIX(progNodeSize-1);
            uiint nProgMaxID = pNode->getID();
            pProgMesh->initBucketNode(nProgMaxID+1, pMesh->getMaxNodeID());
            pProgMesh->setupBucketNode();              
            uiint progElemSize = pProgMesh->getElementSize();
            pElem = pProgMesh->getElementIX(progElemSize-1);
            nProgMaxID = pElem->getID();
            pProgMesh->initBucketElement(nProgMaxID+1, 0);
            pProgMesh->setupBucketElement();
            uiint iGrp, nNumOfGrp = pMesh->getNumOfElemGrp();
            CElementGroup *pElemGrp, *pProgElemGrp;
            for(iGrp=0; iGrp < nNumOfGrp; iGrp++){
                pProgElemGrp = new CElementGroup;
                pElemGrp = pMesh->getElemGrpIX(iGrp);
                pElemGrp->refine(pProgElemGrp);
                pProgElemGrp->setMesh(pProgMesh);
                pProgElemGrp->setID(pElemGrp->getID());
                pProgElemGrp->setName(pElemGrp->getName());
                pProgMesh->addElemGrp(pProgElemGrp);
            };
            pProgMesh->setupAggregate(ilevel+1);
            numOfCommMesh= pMesh->getNumOfCommMesh();
            for(icommesh=0; icommesh< numOfCommMesh; icommesh++){
                pCommMesh= pMesh->getCommMesh(icommesh);
                CIndexBucket *pBucket= pProgMesh->getBucket();
                pProgCommMesh= new CCommMesh(pBucket);
                pProgCommMesh->setCommID( pCommMesh->getCommID());
                pProgCommMesh->setRankID( pCommMesh->getRankID());
                pProgCommMesh->setTransmitRankID( pCommMesh->getTransmitRankID());
                pProgMesh->setCommMesh(pProgCommMesh);
                numOfCommElemAll= pCommMesh->getNumOfCommElementAll();
                pProgCommMesh->reserveCommElementAll(numOfCommElemAll*8);
                numOfCommNode= pCommMesh->getNumOfNode();
                pProgCommMesh->reserveNode(numOfCommNode*8);
                for(icomelem=0; icomelem< numOfCommElemAll; icomelem++){
                    pCommElem= pCommMesh->getCommElementAll(icomelem);
                    pCommElem->setupProgNodeRank(ilevel+1);
                    vProgCommElem.clear();
                    GeneProgCommElem(pCommElem, vProgCommElem);
                    for(iprocom=0; iprocom< vProgCommElem.size(); iprocom++){
                        pProgCommMesh->setCommElementAll(vProgCommElem[iprocom]);
                    };
                };
                pProgCommMesh->AllocateCommElement();
                pProgCommMesh->setupAggCommElement(pProgMesh->getElements());
                pProgCommMesh->sortCommNodeIndex();
                pProgCommMesh->setupMapID2CommID();
            };
            pProgMesh->sortMesh();
            pProgMesh->setupBucketNode();
            pProgMesh->setupBucketElement();
        };
    };
    pAssy= mpGMGModel->getAssyModel(mMGLevel);
    numOfMesh= pAssy->getNumOfMesh();
    for(imesh=0; imesh< numOfMesh; imesh++){
        pMesh= pAssy->getMesh(imesh);
        pMesh->setupEdgeElement(NULL, mMGLevel);
        pMesh->replaceEdgeNode();
    };
}
void CMeshFactory::GeneProgElem(const uiint& ilevel,CElement* pElem, vector<CElement*>& vProgElem, uiint& elementID, CMesh* pProgMesh)
{
    CElement *pProgElem;
    uiint i;
    switch(pElem->getType()){
        case(ElementType::Hexa):
            vProgElem.reserve(8);
            for(i=0; i< 8; i++){
                pProgElem= new CHexa; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();
                vProgElem.push_back(pProgElem);
            };
            dividHexa(pElem,vProgElem, elementID, pProgMesh);
            break;
        case(ElementType::Hexa2):
            vProgElem.reserve(8);
            for(i=0; i< 8; i++){
                pProgElem= new CHexa2; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();
                vProgElem.push_back(pProgElem);
            };
            dividHexa(pElem,vProgElem, elementID, pProgMesh);
            break;
        case(ElementType::Tetra):
            vProgElem.reserve(4);
            for(i=0; i< 4; i++){
                pProgElem= new CHexa; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();
                vProgElem.push_back(pProgElem);
            };
            dividTetra(pElem,vProgElem, elementID, pProgMesh);
            break;
        case(ElementType::Tetra2):
            vProgElem.reserve(4);
            for(i=0; i< 4; i++){
                pProgElem= new CHexa2; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();
                vProgElem.push_back(pProgElem);
            };
            dividTetra(pElem,vProgElem, elementID, pProgMesh);
            break;
        case(ElementType::Prism):
            vProgElem.reserve(6);
            for(i=0; i< 6; i++){
                pProgElem= new CHexa; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();
                vProgElem.push_back(pProgElem);
            };
            dividPrism(pElem,vProgElem, elementID,pProgMesh);
            break;
        case(ElementType::Prism2):
            vProgElem.reserve(6);
            for(i=0; i< 6; i++){
                pProgElem= new CHexa2; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();
                vProgElem.push_back(pProgElem);
            };
            dividPrism(pElem,vProgElem, elementID,pProgMesh);
            break;
        case(ElementType::Quad):
            vProgElem.reserve(4);
            for(i=0; i< 4; i++){
                pProgElem= new CQuad; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();
                vProgElem.push_back(pProgElem);
            };
            dividQuad(pElem,vProgElem, elementID,pProgMesh);
            break;
        case(ElementType::Quad2):
            vProgElem.reserve(4);
            for(i=0; i< 4; i++){
                pProgElem= new CQuad2; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();
                vProgElem.push_back(pProgElem);
            };
            dividQuad(pElem,vProgElem, elementID,pProgMesh);
            break;
        case(ElementType::Triangle):
            vProgElem.reserve(3);
            for(i=0; i< 3; i++){
                pProgElem= new CQuad; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();
                vProgElem.push_back(pProgElem);
            };
            dividTriangle(pElem,vProgElem, elementID,pProgMesh);
            break;
        case(ElementType::Triangle2):
            vProgElem.reserve(3);
            for(i=0; i< 3; i++){
                pProgElem= new CQuad2; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();
                vProgElem.push_back(pProgElem);
            };
            dividTriangle(pElem,vProgElem, elementID,pProgMesh);
            break;
        case(ElementType::Beam):
            vProgElem.reserve(2);
            for(i=0; i< 2; i++){
                pProgElem= new CBeam; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();
                vProgElem.push_back(pProgElem);
            };
            dividBeam(pElem,vProgElem, elementID,pProgMesh);
            break;
        case(ElementType::Beam2):
            vProgElem.reserve(2);
            for(i=0; i< 2; i++){
                pProgElem= new CBeam2; pProgElem->setMGLevel(ilevel+1);
                pProgElem->initialize();
                vProgElem.push_back(pProgElem);
            };
            dividBeam(pElem,vProgElem, elementID,pProgMesh);
            break;
    }
}
void CMeshFactory::dividHexa(CElement* pElem, vector<CElement*>& vProgElem, uiint& elementID, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;
    vector<CNode*> vEdgeNode;
    vector<CNode*> vFaceNode;
    CNode          *pVolNode;
    uiint i;
    vVertNode.resize(8); for(i=0; i< 8; i++){ vVertNode[i] = pElem->getNode(i);}
    vEdgeNode.resize(12);for(i=0; i< 12; i++){ vEdgeNode[i] = pElem->getEdgeInterNode(i);}
    vFaceNode.resize(6); for(i=0; i< 6; i++){ vFaceNode[i] = pElem->getFaceNode(i);}
    pVolNode = pElem->getVolumeNode();
    vProgElem[0]->setNode(vVertNode[0],0); vProgElem[0]->setNode(vEdgeNode[0],1);
    vProgElem[0]->setNode(vFaceNode[0],2); vProgElem[0]->setNode(vEdgeNode[3],3);
    vProgElem[0]->setNode(vEdgeNode[8],4); vProgElem[0]->setNode(vFaceNode[4],5);
    vProgElem[0]->setNode(pVolNode,    6); vProgElem[0]->setNode(vFaceNode[3],7);
    pElem->setProgElem(vProgElem[0], 0);
    vProgElem[1]->setNode(vEdgeNode[0],0); vProgElem[1]->setNode(vVertNode[1],1);
    vProgElem[1]->setNode(vEdgeNode[1],2); vProgElem[1]->setNode(vFaceNode[0],3);
    vProgElem[1]->setNode(vFaceNode[4],4); vProgElem[1]->setNode(vEdgeNode[9],5);
    vProgElem[1]->setNode(vFaceNode[2],6); vProgElem[1]->setNode(pVolNode,    7);
    pElem->setProgElem(vProgElem[1], 1);
    vProgElem[2]->setNode(vEdgeNode[8],0); vProgElem[2]->setNode(vFaceNode[4],1);
    vProgElem[2]->setNode(pVolNode,    2); vProgElem[2]->setNode(vFaceNode[3],3);
    vProgElem[2]->setNode(vVertNode[4],4); vProgElem[2]->setNode(vEdgeNode[4],5);
    vProgElem[2]->setNode(vFaceNode[1],6); vProgElem[2]->setNode(vEdgeNode[7],7);
    pElem->setProgElem(vProgElem[2], 4);
    vProgElem[3]->setNode(vFaceNode[4],0); vProgElem[3]->setNode(vEdgeNode[9],1);
    vProgElem[3]->setNode(vFaceNode[2],2); vProgElem[3]->setNode(pVolNode,    3);
    vProgElem[3]->setNode(vEdgeNode[4],4); vProgElem[3]->setNode(vVertNode[5],5);
    vProgElem[3]->setNode(vEdgeNode[5],6); vProgElem[3]->setNode(vFaceNode[1],7);
    pElem->setProgElem(vProgElem[3], 5);
    vProgElem[4]->setNode(vEdgeNode[3],0); vProgElem[4]->setNode(vFaceNode[0],1);
    vProgElem[4]->setNode(vEdgeNode[2],2); vProgElem[4]->setNode(vVertNode[3],3);
    vProgElem[4]->setNode(vFaceNode[3],4); vProgElem[4]->setNode(pVolNode,    5);
    vProgElem[4]->setNode(vFaceNode[5],6); vProgElem[4]->setNode(vEdgeNode[11],7);
    pElem->setProgElem(vProgElem[4], 3);
    vProgElem[5]->setNode(vFaceNode[0],0); vProgElem[5]->setNode(vEdgeNode[1],1);
    vProgElem[5]->setNode(vVertNode[2],2); vProgElem[5]->setNode(vEdgeNode[2],3);
    vProgElem[5]->setNode(pVolNode,    4); vProgElem[5]->setNode(vFaceNode[2],5);
    vProgElem[5]->setNode(vEdgeNode[10],6);vProgElem[5]->setNode(vFaceNode[5],7);
    pElem->setProgElem(vProgElem[5], 2);
    vProgElem[6]->setNode(vFaceNode[3],0); vProgElem[6]->setNode(pVolNode,    1);
    vProgElem[6]->setNode(vFaceNode[5],2); vProgElem[6]->setNode(vEdgeNode[11],3);
    vProgElem[6]->setNode(vEdgeNode[7],4); vProgElem[6]->setNode(vFaceNode[1],5);
    vProgElem[6]->setNode(vEdgeNode[6],6); vProgElem[6]->setNode(vVertNode[7],7);
    pElem->setProgElem(vProgElem[6], 7);
    vProgElem[7]->setNode(pVolNode,    0); vProgElem[7]->setNode(vFaceNode[2],1);
    vProgElem[7]->setNode(vEdgeNode[10],2); vProgElem[7]->setNode(vFaceNode[5],3);
    vProgElem[7]->setNode(vFaceNode[1],4); vProgElem[7]->setNode(vEdgeNode[5],5);
    vProgElem[7]->setNode(vVertNode[6],6); vProgElem[7]->setNode(vEdgeNode[6],7);
    pElem->setProgElem(vProgElem[7], 6);
    for(i=0; i< 8; i++){
        vProgElem[i]->setID(elementID);          
        ++elementID;
        pProgMesh->setElement(vProgElem[i]);
    };
    uiint iface;
    if(pElem->isMPCMaster()){
        for(iface=0; iface< 6; iface++){
            if(pElem->isMPCFace(iface)){
                switch(iface){
                    case(0):
                        setProgHexaMPCMaster(vProgElem, iface, 0,1,4,5);
                        break;
                    case(1):
                        setProgHexaMPCMaster(vProgElem, iface, 2,3,6,7);
                        break;
                    case(2):
                        setProgHexaMPCMaster(vProgElem, iface, 1,3,5,7);
                        break;
                    case(3):
                        setProgHexaMPCMaster(vProgElem, iface, 0,2,4,6);
                        break;
                    case(4):
                        setProgHexaMPCMaster(vProgElem, iface, 0,1,2,3);
                        break;
                    case(5):
                        setProgHexaMPCMaster(vProgElem, iface, 4,5,6,7);
                        break;
                }
            }
        };
    }
    if(pElem->isMPCSlave()){
        for(iface=0; iface< 6; iface++){
            if(pElem->isMPCFace(iface)){
                switch(iface){
                    case(0):
                        setProgHexaMPCSlave(vProgElem, iface, 0,1,4,5);
                        break;
                    case(1):
                        setProgHexaMPCSlave(vProgElem, iface, 2,3,6,7);
                        break;
                    case(2):
                        setProgHexaMPCSlave(vProgElem, iface, 1,3,5,7);
                        break;
                    case(3):
                        setProgHexaMPCSlave(vProgElem, iface, 0,2,4,6);
                        break;
                    case(4):
                        setProgHexaMPCSlave(vProgElem, iface, 0,1,2,3);
                        break;
                    case(5):
                        setProgHexaMPCSlave(vProgElem, iface, 4,5,6,7);
                        break;
                }
            }
        };
    }
    if(pElem->isCommMesh2()){
        for(iface=0; iface< 6; iface++){
            if(pElem->isCommEntity(iface)){
                switch(iface){
                    case(0):
                        setProgHexaCommMesh2Entity(vProgElem, iface, 0,1,4,5);
                        break;
                    case(1):
                        setProgHexaCommMesh2Entity(vProgElem, iface, 2,3,6,7);
                        break;
                    case(2):
                        setProgHexaCommMesh2Entity(vProgElem, iface, 1,3,5,7);
                        break;
                    case(3):
                        setProgHexaCommMesh2Entity(vProgElem, iface, 0,2,4,6);
                        break;
                    case(4):
                        setProgHexaCommMesh2Entity(vProgElem, iface, 0,1,2,3);
                        break;
                    case(5):
                        setProgHexaCommMesh2Entity(vProgElem, iface, 4,5,6,7);
                        break;
                }
            }
        };
    }
}
void CMeshFactory::setProgHexaMPCMaster(vector<CElement*>& vProgElem, const uiint& iface, const uiint& i, const uiint& j, const uiint& k, const uiint& l)
{
    vProgElem[i]->markingMPCMaster(); vProgElem[i]->markingMPCFace(iface);
    vProgElem[j]->markingMPCMaster(); vProgElem[j]->markingMPCFace(iface);
    vProgElem[k]->markingMPCMaster(); vProgElem[k]->markingMPCFace(iface);
    vProgElem[l]->markingMPCMaster(); vProgElem[l]->markingMPCFace(iface);
}
void CMeshFactory::setProgHexaMPCSlave(vector<CElement*>& vProgElem, const uiint& iface, const uiint& i, const uiint& j, const uiint& k, const uiint& l)
{
    vProgElem[i]->markingMPCSlave(); vProgElem[i]->markingMPCFace(iface);
    vProgElem[j]->markingMPCSlave(); vProgElem[j]->markingMPCFace(iface);
    vProgElem[k]->markingMPCSlave(); vProgElem[k]->markingMPCFace(iface);
    vProgElem[l]->markingMPCSlave(); vProgElem[l]->markingMPCFace(iface);
}
void CMeshFactory::setProgHexaCommMesh2Entity(vector<CElement*>& vProgElem, const uiint& iface, const uiint& i, const uiint& j, const uiint& k, const uiint& l)
{
    vProgElem[i]->markingCommMesh2(); vProgElem[i]->markingCommEntity(iface);
    vProgElem[j]->markingCommMesh2(); vProgElem[j]->markingCommEntity(iface);
    vProgElem[k]->markingCommMesh2(); vProgElem[k]->markingCommEntity(iface);
    vProgElem[l]->markingCommMesh2(); vProgElem[l]->markingCommEntity(iface);
}
void CMeshFactory::dividTetra(CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;
    vector<CNode*> vEdgeNode;
    vector<CNode*> vFaceNode;
    CNode          *pVolNode;
    uiint numOfVert,numOfEdge,numOfFace;
    numOfVert= NumberOfVertex::Tetra(); numOfFace= NumberOfFace::Tetra(); numOfEdge= NumberOfEdge::Tetra();
    uiint i;
    vVertNode.resize(numOfVert);
    for(i=0; i< numOfVert; i++){
        vVertNode[i] = pElem->getNode(i);
    }
    vEdgeNode.resize(numOfEdge);
    for(i=0; i< numOfEdge; i++){
        vEdgeNode[i] = pElem->getEdgeInterNode(i);
    }
    vFaceNode.resize(numOfFace);
    for(i=0; i< numOfFace; i++){
        vFaceNode[i] = pElem->getFaceNode(i);
    }
    pVolNode = pElem->getVolumeNode();
    vProgElem[0]->setNode(vEdgeNode[2],0); vProgElem[0]->setNode(vVertNode[0],1);
    vProgElem[0]->setNode(vEdgeNode[0],2); vProgElem[0]->setNode(vFaceNode[0],3);
    vProgElem[0]->setNode(vFaceNode[3],4); vProgElem[0]->setNode(vEdgeNode[3],5);
    vProgElem[0]->setNode(vFaceNode[1],6); vProgElem[0]->setNode(pVolNode,    7);
    pElem->setProgElem(vProgElem[0], 0);
    vProgElem[1]->setNode(vVertNode[2],0); vProgElem[1]->setNode(vEdgeNode[2],1);
    vProgElem[1]->setNode(vFaceNode[0],2); vProgElem[1]->setNode(vEdgeNode[1],3);
    vProgElem[1]->setNode(vEdgeNode[5],4); vProgElem[1]->setNode(vFaceNode[3],5);
    vProgElem[1]->setNode(pVolNode,    6); vProgElem[1]->setNode(vFaceNode[2],7);
    pElem->setProgElem(vProgElem[1], 2);
    vProgElem[2]->setNode(vFaceNode[0],0); vProgElem[2]->setNode(vEdgeNode[0],1);
    vProgElem[2]->setNode(vVertNode[1],2); vProgElem[2]->setNode(vEdgeNode[1],3);
    vProgElem[2]->setNode(pVolNode,    4); vProgElem[2]->setNode(vFaceNode[1],5);
    vProgElem[2]->setNode(vEdgeNode[4],6); vProgElem[2]->setNode(vFaceNode[2],7);
    pElem->setProgElem(vProgElem[2], 1);
    vProgElem[3]->setNode(vFaceNode[3],0); vProgElem[3]->setNode(vEdgeNode[3],1);
    vProgElem[3]->setNode(vFaceNode[1],2); vProgElem[3]->setNode(pVolNode,    3);
    vProgElem[3]->setNode(vEdgeNode[5],4); vProgElem[3]->setNode(vVertNode[3],5);
    vProgElem[3]->setNode(vEdgeNode[4],6); vProgElem[3]->setNode(vFaceNode[2],7);
    pElem->setProgElem(vProgElem[3], 3);
    for(i=0; i< 4; i++){
        vProgElem[i]->setID(indexCount);
        ++indexCount;
        pProgMesh->setElement(vProgElem[i]);
    };
    uiint iface;
    if(pElem->isMPCMaster()){
        for(iface=0; iface< 4; iface++){
            if(pElem->isMPCFace(iface)){
                switch(iface){
                    case(0):
                        vProgElem[0]->markingMPCMaster(); vProgElem[0]->markingMPCFace(0);
                        vProgElem[1]->markingMPCMaster(); vProgElem[1]->markingMPCFace(0);
                        vProgElem[2]->markingMPCMaster(); vProgElem[2]->markingMPCFace(0);
                        break;
                    case(1):
                        vProgElem[0]->markingMPCMaster(); vProgElem[0]->markingMPCFace(2);
                        vProgElem[2]->markingMPCMaster(); vProgElem[2]->markingMPCFace(2);
                        vProgElem[3]->markingMPCMaster(); vProgElem[3]->markingMPCFace(2);
                        break;
                    case(2):
                        vProgElem[1]->markingMPCMaster(); vProgElem[1]->markingMPCFace(3);
                        vProgElem[2]->markingMPCMaster(); vProgElem[2]->markingMPCFace(5);
                        vProgElem[3]->markingMPCMaster(); vProgElem[3]->markingMPCFace(1);
                        break;
                    case(3):
                        vProgElem[0]->markingMPCMaster(); vProgElem[0]->markingMPCFace(4);
                        vProgElem[1]->markingMPCMaster(); vProgElem[1]->markingMPCFace(4);
                        vProgElem[3]->markingMPCMaster(); vProgElem[3]->markingMPCFace(4);
                        break;
                }
            }
        };
    }
    if(pElem->isMPCSlave()){
        for(iface=0; iface< 4; iface++){
            if(pElem->isMPCFace(iface)){
                switch(iface){
                    case(0):
                        vProgElem[0]->markingMPCSlave(); vProgElem[0]->markingMPCFace(0);
                        vProgElem[1]->markingMPCSlave(); vProgElem[1]->markingMPCFace(0);
                        vProgElem[2]->markingMPCSlave(); vProgElem[2]->markingMPCFace(0);
                        break;
                    case(1):
                        vProgElem[0]->markingMPCSlave(); vProgElem[0]->markingMPCFace(2);
                        vProgElem[2]->markingMPCSlave(); vProgElem[2]->markingMPCFace(2);
                        vProgElem[3]->markingMPCSlave(); vProgElem[3]->markingMPCFace(2);
                        break;
                    case(2):
                        vProgElem[1]->markingMPCSlave(); vProgElem[1]->markingMPCFace(3);
                        vProgElem[2]->markingMPCSlave(); vProgElem[2]->markingMPCFace(5);
                        vProgElem[3]->markingMPCSlave(); vProgElem[3]->markingMPCFace(1);
                        break;
                    case(3):
                        vProgElem[0]->markingMPCSlave(); vProgElem[0]->markingMPCFace(4);
                        vProgElem[1]->markingMPCSlave(); vProgElem[1]->markingMPCFace(4);
                        vProgElem[3]->markingMPCSlave(); vProgElem[3]->markingMPCFace(4);
                        break;
                }
            }
        };
    }
    if(pElem->isCommMesh2()){
        for(iface=0; iface< 4; iface++){
            if(pElem->isCommEntity(iface)){
                switch(iface){
                    case(0):
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(0);
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(0);
                        vProgElem[2]->markingCommMesh2(); vProgElem[2]->markingCommEntity(0);
                        break;
                    case(1):
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(2);
                        vProgElem[2]->markingCommMesh2(); vProgElem[2]->markingCommEntity(2);
                        vProgElem[3]->markingCommMesh2(); vProgElem[3]->markingCommEntity(2);
                        break;
                    case(2):
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(3);
                        vProgElem[2]->markingCommMesh2(); vProgElem[2]->markingCommEntity(5);
                        vProgElem[3]->markingCommMesh2(); vProgElem[3]->markingCommEntity(1);
                        break;
                    case(3):
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(4);
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(4);
                        vProgElem[3]->markingCommMesh2(); vProgElem[3]->markingCommEntity(4);
                        break;
                }
            }
        };
    }
}
void CMeshFactory::dividPrism(CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;
    vector<CNode*> vEdgeNode;
    vector<CNode*> vFaceNode;
    CNode          *pVolNode;
    uiint numOfVert,numOfEdge,numOfFace;
    numOfVert= NumberOfVertex::Prism(); numOfFace= NumberOfFace::Prism(); numOfEdge= NumberOfEdge::Prism();
    uiint i;
    vVertNode.resize(numOfVert); for(i=0; i< numOfVert; i++){ vVertNode[i] = pElem->getNode(i);}
    vEdgeNode.resize(numOfEdge);for(i=0; i< numOfEdge; i++){ vEdgeNode[i] = pElem->getEdgeInterNode(i);}
    vFaceNode.resize(numOfFace); for(i=0; i< numOfFace; i++){ vFaceNode[i] = pElem->getFaceNode(i);}
    pVolNode = pElem->getVolumeNode();
    vProgElem[0]->setNode(vVertNode[2],0); vProgElem[0]->setNode(vEdgeNode[1],1);
    vProgElem[0]->setNode(vFaceNode[0],2); vProgElem[0]->setNode(vEdgeNode[2],3);
    vProgElem[0]->setNode(vEdgeNode[5],4); vProgElem[0]->setNode(vFaceNode[4],5);
    vProgElem[0]->setNode(pVolNode,    6); vProgElem[0]->setNode(vFaceNode[3],7);
    pElem->setProgElem(vProgElem[0], 2);
    vProgElem[1]->setNode(vEdgeNode[1],0); vProgElem[1]->setNode(vVertNode[0],1);
    vProgElem[1]->setNode(vEdgeNode[0],2); vProgElem[1]->setNode(vFaceNode[0],3);
    vProgElem[1]->setNode(vFaceNode[4],4); vProgElem[1]->setNode(vEdgeNode[3],5);
    vProgElem[1]->setNode(vFaceNode[2],6); vProgElem[1]->setNode(pVolNode,    7);
    pElem->setProgElem(vProgElem[1], 0);
    vProgElem[2]->setNode(vFaceNode[0],0); vProgElem[2]->setNode(vEdgeNode[0],1);
    vProgElem[2]->setNode(vVertNode[1],2); vProgElem[2]->setNode(vEdgeNode[2],3);
    vProgElem[2]->setNode(pVolNode,    4); vProgElem[2]->setNode(vFaceNode[2],5);
    vProgElem[2]->setNode(vEdgeNode[4],6); vProgElem[2]->setNode(vFaceNode[3],7);
    pElem->setProgElem(vProgElem[2], 1);
    vProgElem[3]->setNode(vEdgeNode[5],0); vProgElem[3]->setNode(vFaceNode[4],1);
    vProgElem[3]->setNode(pVolNode,    2); vProgElem[3]->setNode(vFaceNode[3],3);
    vProgElem[3]->setNode(vVertNode[5],4); vProgElem[3]->setNode(vEdgeNode[8],5);
    vProgElem[3]->setNode(vFaceNode[1],6); vProgElem[3]->setNode(vEdgeNode[7],7);
    pElem->setProgElem(vProgElem[3], 5);
    vProgElem[4]->setNode(vFaceNode[4],0); vProgElem[4]->setNode(vEdgeNode[3],1);
    vProgElem[4]->setNode(vFaceNode[2],2); vProgElem[4]->setNode(pVolNode,    3);
    vProgElem[4]->setNode(vEdgeNode[8],4); vProgElem[4]->setNode(vVertNode[3],5);
    vProgElem[4]->setNode(vEdgeNode[6],6); vProgElem[4]->setNode(vFaceNode[1],7);
    pElem->setProgElem(vProgElem[4], 3);
    vProgElem[5]->setNode(pVolNode,    0); vProgElem[5]->setNode(vFaceNode[2],1);
    vProgElem[5]->setNode(vEdgeNode[4],2); vProgElem[5]->setNode(vFaceNode[3],3);
    vProgElem[5]->setNode(vFaceNode[1],4); vProgElem[5]->setNode(vEdgeNode[6],5);
    vProgElem[5]->setNode(vVertNode[4],6); vProgElem[5]->setNode(vEdgeNode[7],7);
    pElem->setProgElem(vProgElem[5], 4);
    for(i=0; i< 6; i++){
        vProgElem[i]->setID(indexCount);
        ++indexCount;
        pProgMesh->setElement(vProgElem[i]);
    };
    uiint iface;
    if(pElem->isMPCMaster()){
        for(iface=0; iface< 5; iface++){
            if(pElem->isMPCFace(iface)){
                switch(iface){
                    case(0):
                        vProgElem[0]->markingMPCMaster(); vProgElem[0]->markingMPCFace(0);
                        vProgElem[1]->markingMPCMaster(); vProgElem[1]->markingMPCFace(0);
                        vProgElem[2]->markingMPCMaster(); vProgElem[2]->markingMPCFace(0);
                        break;
                    case(1):
                        vProgElem[3]->markingMPCMaster(); vProgElem[3]->markingMPCFace(1);
                        vProgElem[4]->markingMPCMaster(); vProgElem[4]->markingMPCFace(1);
                        vProgElem[5]->markingMPCMaster(); vProgElem[5]->markingMPCFace(1);
                        break;
                    case(2):
                        vProgElem[1]->markingMPCMaster(); vProgElem[1]->markingMPCFace(2);
                        vProgElem[2]->markingMPCMaster(); vProgElem[2]->markingMPCFace(2);
                        vProgElem[4]->markingMPCMaster(); vProgElem[4]->markingMPCFace(2);
                        vProgElem[5]->markingMPCMaster(); vProgElem[5]->markingMPCFace(2);
                        break;
                    case(3):
                        vProgElem[0]->markingMPCMaster(); vProgElem[0]->markingMPCFace(3);
                        vProgElem[2]->markingMPCMaster(); vProgElem[2]->markingMPCFace(5);
                        vProgElem[3]->markingMPCMaster(); vProgElem[3]->markingMPCFace(3);
                        vProgElem[5]->markingMPCMaster(); vProgElem[5]->markingMPCFace(5);
                        break;
                    case(4):
                        vProgElem[0]->markingMPCMaster(); vProgElem[0]->markingMPCFace(4);
                        vProgElem[1]->markingMPCMaster(); vProgElem[1]->markingMPCFace(4);
                        vProgElem[3]->markingMPCMaster(); vProgElem[3]->markingMPCFace(4);
                        vProgElem[4]->markingMPCMaster(); vProgElem[4]->markingMPCFace(4);
                        break;
                }
            }
        };
    }
    if(pElem->isMPCSlave()){
        for(iface=0; iface< 5; iface++){
            if(pElem->isMPCFace(iface)){
                switch(iface){
                    case(0):
                        vProgElem[0]->markingMPCSlave(); vProgElem[0]->markingMPCFace(0);
                        vProgElem[1]->markingMPCSlave(); vProgElem[1]->markingMPCFace(0);
                        vProgElem[2]->markingMPCSlave(); vProgElem[2]->markingMPCFace(0);
                        break;
                    case(1):
                        vProgElem[3]->markingMPCSlave(); vProgElem[3]->markingMPCFace(1);
                        vProgElem[4]->markingMPCSlave(); vProgElem[4]->markingMPCFace(1);
                        vProgElem[5]->markingMPCSlave(); vProgElem[5]->markingMPCFace(1);
                        break;
                    case(2):
                        vProgElem[1]->markingMPCSlave(); vProgElem[1]->markingMPCFace(2);
                        vProgElem[2]->markingMPCSlave(); vProgElem[2]->markingMPCFace(2);
                        vProgElem[4]->markingMPCSlave(); vProgElem[4]->markingMPCFace(2);
                        vProgElem[5]->markingMPCSlave(); vProgElem[5]->markingMPCFace(2);
                        break;
                    case(3):
                        vProgElem[0]->markingMPCSlave(); vProgElem[0]->markingMPCFace(3);
                        vProgElem[2]->markingMPCSlave(); vProgElem[2]->markingMPCFace(5);
                        vProgElem[3]->markingMPCSlave(); vProgElem[3]->markingMPCFace(3);
                        vProgElem[5]->markingMPCSlave(); vProgElem[5]->markingMPCFace(5);
                        break;
                    case(4):
                        vProgElem[0]->markingMPCSlave(); vProgElem[0]->markingMPCFace(4);
                        vProgElem[1]->markingMPCSlave(); vProgElem[1]->markingMPCFace(4);
                        vProgElem[3]->markingMPCSlave(); vProgElem[3]->markingMPCFace(4);
                        vProgElem[4]->markingMPCSlave(); vProgElem[4]->markingMPCFace(4);
                        break;
                }
            }
        };
    }
    if(pElem->isCommMesh2()){
        for(iface=0; iface< 5; iface++){
            if(pElem->isCommEntity(iface)){
                switch(iface){
                    case(0):
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(0);
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(0);
                        vProgElem[2]->markingCommMesh2(); vProgElem[2]->markingCommEntity(0);
                        break;
                    case(1):
                        vProgElem[3]->markingCommMesh2(); vProgElem[3]->markingCommEntity(1);
                        vProgElem[4]->markingCommMesh2(); vProgElem[4]->markingCommEntity(1);
                        vProgElem[5]->markingCommMesh2(); vProgElem[5]->markingCommEntity(1);
                        break;
                    case(2):
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(2);
                        vProgElem[2]->markingCommMesh2(); vProgElem[2]->markingCommEntity(2);
                        vProgElem[4]->markingCommMesh2(); vProgElem[4]->markingCommEntity(2);
                        vProgElem[5]->markingCommMesh2(); vProgElem[5]->markingCommEntity(2);
                        break;
                    case(3):
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(3);
                        vProgElem[2]->markingCommMesh2(); vProgElem[2]->markingCommEntity(5);
                        vProgElem[3]->markingCommMesh2(); vProgElem[3]->markingCommEntity(3);
                        vProgElem[5]->markingCommMesh2(); vProgElem[5]->markingCommEntity(5);
                        break;
                    case(4):
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(4);
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(4);
                        vProgElem[3]->markingCommMesh2(); vProgElem[3]->markingCommEntity(4);
                        vProgElem[4]->markingCommMesh2(); vProgElem[4]->markingCommEntity(4);
                        break;
                }
            }
        };
    }
}
void CMeshFactory::dividPyramid(CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;
    vector<CNode*> vEdgeNode;
    vector<CNode*> vFaceNode;
    CNode          *pVolNode;
    uiint numOfVert,numOfEdge,numOfFace;
    numOfVert= NumberOfVertex::Pyramid(); numOfFace= NumberOfFace::Pyramid(); numOfEdge= NumberOfEdge::Pyramid();
    uiint i;
    vVertNode.resize(numOfVert); for(i=0; i< numOfVert; i++){ vVertNode[i] = pElem->getNode(i);}
    vEdgeNode.resize(numOfEdge);for(i=0; i< numOfEdge; i++){ vEdgeNode[i] = pElem->getEdgeInterNode(i);}
    vFaceNode.resize(numOfFace); for(i=0; i< numOfFace; i++){ vFaceNode[i] = pElem->getFaceNode(i);}
    pVolNode = pElem->getVolumeNode();
    vProgElem[0]->setNode(vVertNode[0],0); vProgElem[0]->setNode(vEdgeNode[0],1);
    vProgElem[0]->setNode(vFaceNode[0],2); vProgElem[0]->setNode(vEdgeNode[3],3);
    vProgElem[0]->setNode(vEdgeNode[7],4); vProgElem[0]->setNode(vFaceNode[4],5);
    vProgElem[0]->setNode(pVolNode,    6); vProgElem[0]->setNode(vFaceNode[3],7);
    pElem->setProgElem(vProgElem[0], 0);
    vProgElem[1]->setNode(vEdgeNode[0],0); vProgElem[1]->setNode(vVertNode[1],1);
    vProgElem[1]->setNode(vEdgeNode[1],2); vProgElem[1]->setNode(vFaceNode[0],3);
    vProgElem[1]->setNode(vFaceNode[4],4); vProgElem[1]->setNode(vEdgeNode[4],5);
    vProgElem[1]->setNode(vFaceNode[1],6); vProgElem[1]->setNode(pVolNode,    7);
    pElem->setProgElem(vProgElem[1], 1);
    vProgElem[2]->setNode(vFaceNode[0],0); vProgElem[2]->setNode(vEdgeNode[1],1);
    vProgElem[2]->setNode(vVertNode[2],2); vProgElem[2]->setNode(vEdgeNode[2],3);
    vProgElem[2]->setNode(pVolNode,    4); vProgElem[2]->setNode(vFaceNode[1],5);
    vProgElem[2]->setNode(vEdgeNode[5],6); vProgElem[2]->setNode(vFaceNode[2],7);
    pElem->setProgElem(vProgElem[2], 2);
    vProgElem[3]->setNode(vEdgeNode[3],0); vProgElem[3]->setNode(vFaceNode[0],1);
    vProgElem[3]->setNode(vEdgeNode[2],2); vProgElem[3]->setNode(vVertNode[3],3);
    vProgElem[3]->setNode(vFaceNode[3],4); vProgElem[3]->setNode(pVolNode,    5);
    vProgElem[3]->setNode(vFaceNode[2],6); vProgElem[3]->setNode(vEdgeNode[6],7);
    pElem->setProgElem(vProgElem[3], 3);
    vProgElem[4]->setNode(vEdgeNode[7],0); vProgElem[4]->setNode(vVertNode[4],1);
    vProgElem[4]->setNode(vEdgeNode[4],2); vProgElem[4]->setNode(vFaceNode[4],3);
    vProgElem[4]->setNode(pVolNode,    4);
    pElem->setProgElem(vProgElem[4], 7);
    vProgElem[5]->setNode(vEdgeNode[4],0); vProgElem[5]->setNode(vVertNode[4],1);
    vProgElem[5]->setNode(vEdgeNode[5],2); vProgElem[5]->setNode(vFaceNode[1],3);
    vProgElem[5]->setNode(pVolNode,    4);
    pElem->setProgElem(vProgElem[5], 4);
    vProgElem[6]->setNode(vEdgeNode[5],0); vProgElem[6]->setNode(vVertNode[4],1);
    vProgElem[6]->setNode(vEdgeNode[6],2); vProgElem[6]->setNode(vFaceNode[2],3);
    vProgElem[6]->setNode(pVolNode,    4);
    pElem->setProgElem(vProgElem[6], 5);
    vProgElem[7]->setNode(vEdgeNode[6],0); vProgElem[7]->setNode(vVertNode[4],1);
    vProgElem[7]->setNode(vEdgeNode[7],2); vProgElem[7]->setNode(vFaceNode[3],3);
    vProgElem[7]->setNode(pVolNode,    4);
    pElem->setProgElem(vProgElem[7], 6);
    for(i=0; i< 8; i++){
        vProgElem[i]->setID(indexCount);
        ++indexCount;
        pProgMesh->setElement(vProgElem[i]);
    };
}
void CMeshFactory::dividQuad(CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;
    vector<CNode*> vEdgeNode;
    vector<CNode*> vFaceNode;
    uiint numOfVert,numOfEdge,numOfFace;
    numOfVert= NumberOfVertex::Quad(); numOfFace= NumberOfFace::Quad(); numOfEdge= NumberOfEdge::Quad();
    uiint i;
    vVertNode.resize(numOfVert); for(i=0; i< numOfVert; i++){ vVertNode[i] = pElem->getNode(i);}
    vEdgeNode.resize(numOfEdge);for(i=0; i< numOfEdge; i++){ vEdgeNode[i] = pElem->getEdgeInterNode(i);}
    vFaceNode.resize(numOfFace); for(i=0; i< numOfFace; i++){ vFaceNode[i] = pElem->getFaceNode(i);}
    vProgElem[0]->setNode(vVertNode[0],0); vProgElem[0]->setNode(vEdgeNode[0],1);
    vProgElem[0]->setNode(vFaceNode[0],2); vProgElem[0]->setNode(vEdgeNode[3],3);
    pElem->setProgElem(vProgElem[0], 0);
    vProgElem[1]->setNode(vEdgeNode[0],0); vProgElem[1]->setNode(vVertNode[1],1);
    vProgElem[1]->setNode(vEdgeNode[1],2); vProgElem[1]->setNode(vFaceNode[0],3);
    pElem->setProgElem(vProgElem[1], 1);
    vProgElem[2]->setNode(vEdgeNode[1],0); vProgElem[2]->setNode(vVertNode[2],1);
    vProgElem[2]->setNode(vEdgeNode[2],2); vProgElem[2]->setNode(vFaceNode[0],3);
    pElem->setProgElem(vProgElem[2], 2);
    vProgElem[3]->setNode(vEdgeNode[2],0); vProgElem[3]->setNode(vVertNode[3],1);
    vProgElem[3]->setNode(vEdgeNode[3],2); vProgElem[3]->setNode(vFaceNode[0],3);
    pElem->setProgElem(vProgElem[3], 3);
    for(i=0; i< 4; i++){
        vProgElem[i]->setID(indexCount);
        ++indexCount;
        pProgMesh->setElement(vProgElem[i]);
    };
    uiint iprog;
    if(pElem->isMPCMaster()){
        for(iprog=0; iprog< 4; iprog++){ vProgElem[iprog]->markingMPCMaster(); vProgElem[iprog]->markingMPCFace(0);}
    }
    if(pElem->isMPCSlave()){
        for(iprog=0; iprog< 4; iprog++){ vProgElem[iprog]->markingMPCSlave(); vProgElem[iprog]->markingMPCFace(0);}
    }
    uiint iedge;
    if(pElem->isCommMesh2()){
        for(iedge=0; iedge< 4; iedge++){
            if(pElem->isCommEntity(iedge)){
                switch(iedge){
                    case(0):
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(0);
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(0);
                        break;
                    case(1):
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(1);
                        vProgElem[2]->markingCommMesh2(); vProgElem[2]->markingCommEntity(0);
                        break;
                    case(2):
                        vProgElem[2]->markingCommMesh2(); vProgElem[2]->markingCommEntity(1);
                        vProgElem[3]->markingCommMesh2(); vProgElem[3]->markingCommEntity(0);
                        break;
                    case(3):
                        vProgElem[3]->markingCommMesh2(); vProgElem[3]->markingCommEntity(1);
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(3);
                        break;
                }
            }
        };
    }
}
void CMeshFactory::dividTriangle(CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;
    vector<CNode*> vEdgeNode;
    vector<CNode*> vFaceNode;
    uiint numOfVert,numOfEdge,numOfFace;
    numOfVert= NumberOfVertex::Triangle(); numOfFace= NumberOfFace::Triangle(); numOfEdge= NumberOfEdge::Triangle();
    uiint i;
    vVertNode.resize(numOfVert); for(i=0; i< numOfVert; i++){ vVertNode[i] = pElem->getNode(i);}
    vEdgeNode.resize(numOfEdge);for(i=0; i< numOfEdge; i++){ vEdgeNode[i] = pElem->getEdgeInterNode(i);}
    vFaceNode.resize(numOfFace); for(i=0; i< numOfFace; i++){ vFaceNode[i] = pElem->getFaceNode(i);}
    vProgElem[0]->setNode(vEdgeNode[0],0); vProgElem[0]->setNode(vVertNode[1],1);
    vProgElem[0]->setNode(vEdgeNode[1],2); vProgElem[0]->setNode(vFaceNode[0],3);
    pElem->setProgElem(vProgElem[0], 1);
    vProgElem[1]->setNode(vEdgeNode[1],0); vProgElem[1]->setNode(vVertNode[2],1);
    vProgElem[1]->setNode(vEdgeNode[2],2); vProgElem[1]->setNode(vFaceNode[0],3);
    pElem->setProgElem(vProgElem[1], 2);
    vProgElem[2]->setNode(vEdgeNode[2],0); vProgElem[2]->setNode(vVertNode[0],1);
    vProgElem[2]->setNode(vEdgeNode[0],2); vProgElem[2]->setNode(vFaceNode[0],3);
    pElem->setProgElem(vProgElem[2], 0);
    for(i=0; i< 3; i++){
        vProgElem[i]->setID(indexCount);
        ++indexCount;
        pProgMesh->setElement(vProgElem[i]);
    };
    uiint iprog;
    if(pElem->isMPCMaster()){
        for(iprog=0; iprog< 3; iprog++){ vProgElem[iprog]->markingMPCMaster(); vProgElem[iprog]->markingMPCFace(0);}
    }
    if(pElem->isMPCSlave()){
        for(iprog=0; iprog< 3; iprog++){ vProgElem[iprog]->markingMPCSlave(); vProgElem[iprog]->markingMPCFace(0);}
    }
    uiint iedge;
    if(pElem->isCommMesh2()){
        for(iedge=0; iedge< 3; iedge++){
            if(pElem->isCommEntity(iedge)){
                switch(iedge){
                    case(0):
                        vProgElem[2]->markingCommMesh2(); vProgElem[2]->markingCommEntity(1);
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(0);
                        break;
                    case(1):
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(1);
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(0);
                        break;
                    case(2):
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(1);
                        vProgElem[2]->markingCommMesh2(); vProgElem[2]->markingCommEntity(0);
                        break;
                }
            }
        };
    }
}
void CMeshFactory::dividBeam(CElement* pElem, vector<CElement*>& vProgElem, uiint& indexCount, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;
    vector<CNode*> vEdgeNode;
    uiint numOfVert,numOfEdge;
    numOfVert= NumberOfVertex::Beam(); numOfEdge= NumberOfEdge::Beam();
    uiint i;
    vVertNode.resize(numOfVert); for(i=0; i< numOfVert; i++){ vVertNode[i] = pElem->getNode(i);}
    vEdgeNode.resize(numOfEdge); for(i=0; i< numOfEdge; i++){ vEdgeNode[i] = pElem->getEdgeInterNode(i);}
    vProgElem[0]->setNode(vVertNode[0],0); vProgElem[0]->setNode(vEdgeNode[0],1);
    pElem->setProgElem(vProgElem[0], 0);
    vProgElem[1]->setNode(vEdgeNode[0],0); vProgElem[1]->setNode(vVertNode[1],1);
    pElem->setProgElem(vProgElem[1], 1);
    for(i=0; i< 2; i++){
        vProgElem[i]->setID(indexCount);
        ++indexCount;
        pProgMesh->setElement(vProgElem[i]);
    };
    uiint iprog;
    if(pElem->isMPCMaster()){
        for(iprog=0; iprog< 2; iprog++) vProgElem[iprog]->markingMPCMaster();
    }
    if(pElem->isMPCSlave()){
        for(iprog=0; iprog< 2; iprog++) vProgElem[iprog]->markingMPCSlave();
    }
    uiint ivert;
    if(pElem->isCommMesh2()){
        for(ivert=0; ivert< 2; ivert++){
            if(pElem->isCommEntity(ivert)){
                switch(ivert){
                    case(0):
                        vProgElem[0]->markingCommMesh2(); vProgElem[0]->markingCommEntity(0);
                        break;
                    case(1):
                        vProgElem[1]->markingCommMesh2(); vProgElem[1]->markingCommEntity(1);
                        break;
                }
            }
        };
    }
}
void CMeshFactory::setupBucketMesh(const uiint& mgLevel, const uiint& maxID, const uiint& minID)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTAssyModel->intializeBucket(maxID, minID);
    mpTAssyModel->setMaxMeshID(maxID);
    mpTAssyModel->setMinMeshID(minID);
}
void CMeshFactory::GeneAssyModel(const uiint& nNumOfLevel)
{
    if(nNumOfLevel < 2) return;
    mpGMGModel->resizeAssyModel(nNumOfLevel);
    for(uiint i=1; i< nNumOfLevel; i++){
        mpTAssyModel = new CAssyModel();
        mpTAssyModel->setMGLevel(i);
        mpGMGModel->addModel(mpTAssyModel,i);
    };
}
void CMeshFactory::reserveMesh(const uiint& mgLevel, const uiint& num_of_mesh)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTAssyModel->resizeMesh(num_of_mesh);
}
void CMeshFactory::GeneMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& index, const uiint& nProp)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = new CMesh();
    mpTMesh->setMeshID(mesh_id);
    mpTMesh->setMGLevel(mgLevel);
    mpTMesh->setMaxMGLevel(mMGLevel);
    mpTMesh->setSolutionType(mnSolutionType);
    mpTMesh->setProp(nProp);
    mpTAssyModel->setBucket(mesh_id, index);
    mpTAssyModel->setMesh(mpTMesh, index);
    CBNodeMeshGrp *pBNodeMeshGrp= new CBNodeMeshGrp;
    mpTMesh->setBNodeMeshGrp(pBNodeMeshGrp);
}
void CMeshFactory::GeneNode(const uiint& mgLevel, const uiint& mesh_id, const uiint& id, const vdouble& coord,
                            const uiint& nodeType, const uiint& nNumOfSDOF, const uiint& nNumOfVDOF)
{
    CNode *pNode;
    switch(nodeType){
        case(NodeType::Scalar):
            pNode = new CScalarNode();
            pNode->setScalarDOF(nNumOfSDOF);
            break;
        case(NodeType::Vector):
            pNode = new CVectorNode();
            pNode->setVectorDOF(nNumOfVDOF);
            break;
        case(NodeType::ScalarVector):
            pNode = new CScalarVectorNode();
            pNode->setScalarDOF(nNumOfSDOF);
            pNode->setVectorDOF(nNumOfVDOF);
            break;
        default:
            break;
    }
    pNode->setID(id);
    pNode->setCoord(coord);
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    if(!mpTAssyModel) mpLogger->Info(Utility::LoggerMode::MWDebug, "AssyModel => NULL");
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->setNode(pNode);
}
void CMeshFactory::setupNode(const uiint& mgLevel, const uiint& mesh_id)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh(mesh_id);
    mpTMesh->setupNumOfNode();
}
void CMeshFactory::reserveNode(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_node)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->reserveNode(num_of_node);
}
void CMeshFactory::GeneElement(const uiint& mgLevel, const uiint& mesh_id, const uiint& id, const uiint& type_num, const vint& node_id)
{
    CElement *pElement;
    switch(type_num){
        case(ElementType::Hexa):
            pElement = new CHexa;
            break;
        case(ElementType::Tetra):
            pElement = new CTetra;
            break;
        case(ElementType::Prism):
            pElement = new CPrism;
            break;
        case(ElementType::Quad):
            pElement = new CQuad;
            break;
        case(ElementType::Triangle):
            pElement = new CTriangle;
            break;
        case(ElementType::Beam):
            pElement = new CBeam;
            break;
        case(ElementType::Hexa2):
            pElement = new CHexa2;
            break;
        case(ElementType::Tetra2):
            pElement = new CTetra2;
            break;
        case(ElementType::Prism2):
            pElement = new CPrism2;
            break;
        case(ElementType::Quad2):
            pElement = new CQuad2;
            break;
        case(ElementType::Triangle2):
            pElement = new CTriangle2;
            break;
        case(ElementType::Beam2):
            pElement = new CBeam2;
            break;
        default:
            mpLogger->Info(Utility::LoggerMode::Error,"Error::GeneElement at Factory");
            break;
    }
    pElement->initialize();
    pElement->setID(id);
    CNode *pNode;
    uiint i;
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);    
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    for(i=0; i < node_id.size(); i++){
        pNode = mpTMesh->getNode(node_id[i]);
        pElement->setNode(pNode, i);
    };
    if(pElement->getOrder()==ElementOrder::Second){
        uiint iedge;
        uiint nNumOfVert = pElement->getNumOfVert();
        for(iedge=0; iedge < pElement->getNumOfEdge(); iedge++){
            uiint nNodeID = node_id[nNumOfVert + iedge];
            pNode = mpTMesh->getNode(nNodeID);
            pElement->setEdgeInterNode(pNode, iedge);
        };
    }
    mpTMesh->setElement(pElement);
}
void CMeshFactory::setupElement(const uiint& mgLevel, const uiint& mesh_id)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh(mesh_id);
    mpTMesh->setupNumOfElement();
}
void CMeshFactory::reserveElement(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_element)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->reserveElement(num_of_element);
}
void CMeshFactory::resizeAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_node)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->resizeAggregate(num_of_node);
}
void CMeshFactory::GeneAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_node)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    uiint i;
    if(mnSolutionType==SolutionType::FEM){
        for(i=0; i< num_of_node; i++){
            CAggregateElement *pAggElem = new CAggregateElement;
            mpTMesh->setAggElement(pAggElem, i);
        };
    }
    if(mnSolutionType==SolutionType::FVM){
        for(i=0; i< num_of_node; i++){
            CAggregateNode    *pAggNode = new CAggregateNode;
            mpTMesh->setAggNode(pAggNode, i);
        };
    }
}
void CMeshFactory::reserveBoundaryNodeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_bnd)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    mpTMesh->reserveBndNodeMesh(num_of_bnd);
}
void CMeshFactory::GeneBoundaryNodeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id, const uiint& bnd_type, const string& bnd_name)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    CBoundaryNodeMesh *pBNodeMesh= new CBoundaryNodeMesh;
    pBNodeMesh->setID(bnd_id);
    pBNodeMesh->setBndType(bnd_type);
    pBNodeMesh->setName(bnd_name);
    mpTMesh->setBndNodeMesh(pBNodeMesh);
}
uiint CMeshFactory::getNumOfBounaryNodeMesh(const uiint& mgLevel, const uiint& mesh_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    return mpTMesh->getNumOfBoundaryNodeMesh();
}
void CMeshFactory::reserveBoundaryFaceMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_bnd)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    mpTMesh->reserveBndFaceMesh(num_of_bnd);
}
void CMeshFactory::GeneBoundaryFaceMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id, const uiint& bnd_type, const string& bnd_name,
                                        const uiint& numOfDOF, const vuint& vDOF)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    CBoundaryFaceMesh *pBFaceMesh= new CBoundaryFaceMesh;
    pBFaceMesh->setMGLevel(mgLevel);
    pBFaceMesh->setID(bnd_id);
    pBFaceMesh->setBndType(bnd_type);
    pBFaceMesh->setName(bnd_name);
    pBFaceMesh->resizeDOF(numOfDOF);
    if(numOfDOF != vDOF.size()) mpLogger->Info(Utility::LoggerMode::Error, "CMeshFactory::GeneBoundaryFaceMesh, invalid argument");
    uiint idof, dof;
    for(idof=0; idof < vDOF.size(); idof++){
        dof = vDOF[idof];
        pBFaceMesh->setDOF(idof, dof);
    };
    mpTMesh->setBndFaceMesh(pBFaceMesh);
}
uiint CMeshFactory::getNumOfBounaryFaceMesh(const uiint& mgLevel, const uiint& mesh_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    return mpTMesh->getNumOfBoundaryFaceMesh();
}
CBoundaryFaceMesh* CMeshFactory::getBoundaryFaceMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    return mpTMesh->getBndFaceMeshID(bnd_id);
}
void CMeshFactory::reserveBoundaryVolumeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_bnd)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    mpTMesh->reserveBndVolumeMesh(num_of_bnd);
}
void CMeshFactory::GeneBoundaryVolumeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id, const uiint& bnd_type, const string& bnd_name,
                                          const uiint& numOfDOF, const vuint& vDOF)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    CBoundaryVolumeMesh *pBVolMesh= new CBoundaryVolumeMesh;
    pBVolMesh->setMGLevel(mgLevel);
    pBVolMesh->setID(bnd_id);
    pBVolMesh->setBndType(bnd_type);
    pBVolMesh->setName(bnd_name);
    pBVolMesh->resizeDOF(numOfDOF);
    if(numOfDOF != vDOF.size()) mpLogger->Info(Utility::LoggerMode::Error, "CMeshFactory::GeneBoundaryVolumeMesh, invalid argument");
    uiint idof, dof;
    for(idof=0; idof < vDOF.size(); idof++){
        dof = vDOF[idof];
        pBVolMesh->setDOF(idof, dof);
    }
    mpTMesh->setBndVolumeMesh(pBVolMesh);
}
uiint CMeshFactory::getNumOfBounaryVolumeMesh(const uiint& mgLevel, const uiint& mesh_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    return mpTMesh->getNumOfBoundaryVolumeMesh();
}
CBoundaryVolumeMesh* CMeshFactory::getBoundaryVolumeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    return mpTMesh->getBndVolumeMeshID(bnd_id);
}
void CMeshFactory::reserveBoundaryEdgeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& num_of_bnd)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    mpTMesh->reserveBndEdgeMesh(num_of_bnd);
}
void CMeshFactory::GeneBoundaryEdgeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id, const uiint& bnd_type, const string& bnd_name,
                   const uiint& numOfDOF, const vuint& vDOF)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    CBoundaryEdgeMesh *pBEdgeMesh= new CBoundaryEdgeMesh;
    pBEdgeMesh->setMGLevel(mgLevel);
    pBEdgeMesh->setID(bnd_id);
    pBEdgeMesh->setBndType(bnd_type);
    pBEdgeMesh->setName(bnd_name);
    pBEdgeMesh->resizeDOF(numOfDOF);
    if(numOfDOF != vDOF.size()) mpLogger->Info(Utility::LoggerMode::Error, "CMeshFactory::GeneBoundaryEdgeMesh, invalid argument");
    uiint idof, dof;
    for(idof=0; idof < numOfDOF; idof++){
        dof = vDOF[idof];
        pBEdgeMesh->setDOF(idof, dof);
    }
    mpTMesh->setBndEdgeMesh(pBEdgeMesh);
}
uiint CMeshFactory::getNumOfBounaryEdgeMesh(const uiint& mgLevel, const uiint& mesh_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    return mpTMesh->getNumOfBoundaryEdgeMesh();
}
CBoundaryEdgeMesh* CMeshFactory::getBoundaryEdgeMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    return mpTMesh->getBndEdgeMeshID(bnd_id);
}
void CMeshFactory::GeneBoundaryNode(const uiint& mgLevel, const uiint& bnd_id, const uiint& bndType,
                          const uiint& mesh_id, const uiint& node_id,
                          const uiint& b_node_id, const uiint& dof, const double& val)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    CBoundaryNodeMesh* pBndNodeMesh= mpTMesh->getBndNodeMeshID(bnd_id);
    pBndNodeMesh->setBndType(bndType);
    uiint crIndex= pBndNodeMesh->getNumOfBNode();
    uiint crBNodeID;
    CBoundarySBNode *pCrBNode;
    uiint ib; bool bfind(false);
    if(crIndex > 0){
        if(crIndex==1){
            pCrBNode= pBndNodeMesh->getBNodeIX(0);
            crBNodeID= pCrBNode->getID();
            if(b_node_id==crBNodeID) bfind= true;
        }
        for(ib=crIndex-1; ib > 0; ib--){
            pCrBNode= pBndNodeMesh->getBNodeIX(ib);
            crBNodeID= pCrBNode->getID();
            if(b_node_id==crBNodeID) bfind= true;
        };
    }
    if(bfind){
        pCrBNode= pBndNodeMesh->getBNodeID(b_node_id);
        pCrBNode->addDOF(dof);
        pCrBNode->setValue(dof, val);
    }else{
        CNode *pNode= mpTMesh->getNode(node_id);
        CBoundarySBNode *pBNode = new CBoundarySBNode();
        pBNode->setNode(pNode);
        pBNode->setID(b_node_id);
        pBNode->addDOF(dof);
        pBNode->setValue(dof, val);
        pBndNodeMesh->addBNode(pBNode);
    }
}
void CMeshFactory::GeneBoundaryFaceNode(const uiint& mgLevel, const uiint& bnd_id, const uiint& bndType,
        const uiint& mesh_id, const uiint& node_id, const uiint& b_node_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    CBoundaryFaceMesh *pBFaceMesh= mpTMesh->getBndFaceMeshID(bnd_id);
    CNode *pNode= mpTMesh->getNode(node_id);
    CBoundaryNode *pBNode= new CBoundaryNode;
    pBNode->setNode(pNode);
    pBNode->setID(b_node_id);
    pBNode->setMGLevel(mgLevel);
    pBNode->resizeValue(mMGLevel+1);
    pBFaceMesh->addBNode(pBNode);
}
void CMeshFactory::setValue_BoundaryFaceNode(const uiint& mesh_id, const uiint& bnd_id, const uiint& bnode_id, const uiint& dof, const double& val)
{
    mpTAssyModel = mpGMGModel->getAssyModel(0);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    CBoundaryFaceMesh *pBFaceMesh= mpTMesh->getBndFaceMeshID(bnd_id);
    CBoundaryNode *pBNode= pBFaceMesh->getBNodeID(bnode_id);
    pBNode->addValue(dof, 0, val);
}
void CMeshFactory::initFaceAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    CBoundaryFaceMesh *pBFaceMesh= mpTMesh->getBndFaceMeshID(bnd_id);
    pBFaceMesh->setupAggFace();
}
void CMeshFactory::resizeFaceAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    CBoundaryFaceMesh *pBFaceMesh= mpTMesh->getBndFaceMeshID(bnd_id);
    pBFaceMesh->resizeAggFace();
}
void CMeshFactory::GeneBoundaryFace(const uiint& mgLevel, const uiint& bnd_id, const uiint& bndType, const uiint& elemType,
                            const uiint& mesh_id, const uiint& elem_id, const uiint& face_id, vuint& vBNodeID,
                            const uiint& b_face_id, const uiint& dof, const double& val)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    CBoundaryFaceMesh *pBFaceMesh= mpTMesh->getBndFaceMeshID(bnd_id);
    pBFaceMesh->setBndType(bndType);
    pBFaceMesh->setID(bnd_id);
    uiint crIndex= pBFaceMesh->getNumOfBFace();
    uiint crBFaceID;
    CBoundaryFace *pCrBFace;
    uiint ib; bool bfind(false);
    if(crIndex > 0){
        if(crIndex==1){
            pCrBFace= pBFaceMesh->getBFaceIX(0);
            crBFaceID= pCrBFace->getID();
            if(b_face_id==crBFaceID) bfind= true;
        }
        for(ib=crIndex-1; ib > 0; ib--){
            pCrBFace= pBFaceMesh->getBFaceIX(ib);
            crBFaceID= pCrBFace->getID();
            if(b_face_id==crBFaceID) bfind= true;
        };
    }
    CBoundaryFace *pBFace;
    if(bfind){
        pBFace= pBFaceMesh->getBFaceID(b_face_id);
        if(bndType==BoundaryType::Neumann){
            double dArea= pBFace->getArea();
            double dBndValue= val * dArea;
            pBFace->setBndValue(dof, dBndValue);
        }
    }else{
        pBFace = new CBoundaryFace();
        pBFace->setElementID(elem_id);
        pBFace->setElementFaceID(face_id);
        pBFace->setID(b_face_id);
        pBFace->setBFaceShape(elemType);
        CElement *pElem= mpTMesh->getElement(elem_id);
        pBFace->setElement(pElem);
        switch(elemType){
            case(ElementType::Quad):case(ElementType::Quad2):
                pBFace->resizeBNode(4);
                break;
            case(ElementType::Triangle):case(ElementType::Triangle2):
                pBFace->resizeBNode(3);
                break;
            default:
                break;
        }
        pBFaceMesh->addBFace(pBFace);
    }
    if(!bfind){
        CBoundaryNode *pBNode;
        uiint numOfVert= pBFace->getNumOfVert();
        uiint ivert;
        for(ivert=0; ivert < numOfVert; ivert++){
            pBNode= pBFaceMesh->getBNodeID(vBNodeID[ivert]);
            pBFace->setBNode(ivert, pBNode);
        };
        uiint iedge;
        if(pBFace->getOrder()==ElementOrder::Second){
            uiint nNumOfVert = pBFace->getNumOfVert();
            uiint nNumOfEdge = pBFace->getNumOfEdge();
            uiint id;
            for(iedge=0; iedge < nNumOfEdge; iedge++){
                id = vBNodeID[nNumOfVert + iedge];
                pBNode= pBFaceMesh->getBNodeID(id);
                uiint ibnode = pBFaceMesh->getBNodeIndex(id);
                pBFaceMesh->setAggFace(ibnode, pBFace->getID());
                pBFace->setEdgeBNode(iedge, pBNode);
                pBFace->markingEdge(iedge);         
            };
            pBFace->replaceEdgeBNode();
        }
        pBFace->calcArea();
        if(bndType==BoundaryType::Neumann){
            double dArea= pBFace->getArea();
            double dBndValue= val * dArea;
            pBFace->setBndValue(dof, dBndValue);
        }
    }
}
void CMeshFactory::GeneBoundaryVolumeNode(const uiint& mgLevel, const uiint& bnd_id, const uiint& bndType,
        const uiint& mesh_id, const uiint& node_id, const uiint& b_node_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    CBoundaryVolumeMesh *pBVolumeMesh= mpTMesh->getBndVolumeMeshID(bnd_id);
    CNode *pNode= mpTMesh->getNode(node_id);
    CBoundaryNode *pBNode= new CBoundaryNode;
    pBNode->setNode(pNode);
    pBNode->setID(b_node_id);
    pBNode->setMGLevel(mgLevel);
    pBNode->resizeValue(mMGLevel+1);
    pBVolumeMesh->addBNode(pBNode);
}
void CMeshFactory::setValue_BoundaryVolumeNode(const uiint& mesh_id, const uiint& bnd_id, const uiint& bnode_id, const uiint& dof, const double& val)
{
    mpTAssyModel = mpGMGModel->getAssyModel(0);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    CBoundaryVolumeMesh *pBVolMesh= mpTMesh->getBndVolumeMeshID(bnd_id);
    CBoundaryNode *pBNode= pBVolMesh->getBNodeID(bnode_id);
    pBNode->addValue(dof, 0, val);
}
void CMeshFactory::initVolumeAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    CBoundaryVolumeMesh *pBVolumeMesh= mpTMesh->getBndVolumeMeshID(bnd_id);
    pBVolumeMesh->setupAggVol();
}
void CMeshFactory::resizeVolumeAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    CBoundaryVolumeMesh *pBVolumeMesh= mpTMesh->getBndVolumeMeshID(bnd_id);
    pBVolumeMesh->resizeAggVol();
}
void CMeshFactory::GeneBoundaryVolume(const uiint& mgLevel, const uiint& bnd_id, const uiint& bndType, const uiint& elemType,
                            const uiint& mesh_id, const uiint& elem_id, vuint& vBNodeID,
                            const uiint& b_vol_id, const uiint& dof, const double& val)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    CBoundaryVolumeMesh *pBVolumeMesh= mpTMesh->getBndVolumeMeshID(bnd_id);
    pBVolumeMesh->setBndType(bndType);
    pBVolumeMesh->setID(bnd_id);
    uiint crIndex= pBVolumeMesh->getNumOfVolume();
    uiint crBVolumeID;
    CBoundaryVolume *pCrBVolume;
    uiint ib; bool bfind(false);
    if(crIndex > 0){
        if(crIndex==1){
            pCrBVolume= pBVolumeMesh->getBVolumeIX(0);
            crBVolumeID= pCrBVolume->getID();
            if(b_vol_id==crBVolumeID) bfind= true;
        }
        for(ib=crIndex-1; ib > 0; ib--){
            pCrBVolume= pBVolumeMesh->getBVolumeIX(ib);
            crBVolumeID= pCrBVolume->getID();
            if(b_vol_id==crBVolumeID) bfind= true;
        };
    }
    CBoundaryVolume *pBVolume;
    if(bfind){
        pBVolume= pBVolumeMesh->getBVolumeID(b_vol_id);
        if(bndType==BoundaryType::Neumann){
            double dCubicVol= pBVolume->getCubicVolume();
            double dBndValue= dCubicVol*val;
            pBVolume->setBndValue(dof,dBndValue);
        }
    }else{
        switch(elemType){
            case(ElementType::Hexa):
                pBVolume = new CBoundaryHexa;
                pBVolume->setOrder(ElementOrder::First);
                break;
            case(ElementType::Hexa2):
                pBVolume = new CBoundaryHexa;
                pBVolume->setOrder(ElementOrder::Second);
                break;
            case(ElementType::Tetra):
                pBVolume = new CBoundaryTetra;
                pBVolume->setOrder(ElementOrder::First);
                break;
            case(ElementType::Tetra2):
                pBVolume = new CBoundaryTetra;
                pBVolume->setOrder(ElementOrder::Second);
                break;
            case(ElementType::Prism):
                pBVolume = new CBoundaryPrism;
                pBVolume->setOrder(ElementOrder::First);
                break;
            case(ElementType::Prism2):
                pBVolume = new CBoundaryPrism;
                pBVolume->setOrder(ElementOrder::Second);
                break;
            default:
                mpLogger->Info(Utility::LoggerMode::Error, "invalid ElementType, CMeshFactory::GeneBoundaryVolume");
                break;
        }
        pBVolume->setElementID(elem_id);
        pBVolume->setID(b_vol_id);
        CElement *pElem= mpTMesh->getElement(elem_id);
        pBVolume->setElement(pElem);
        pBVolumeMesh->addBVolume(pBVolume);
    }
    if(!bfind){
        CBoundaryNode *pBNode;
        uiint numOfVert= pBVolume->getNumOfVert();
        uiint ivert;
        for(ivert=0; ivert < numOfVert; ivert++){
            pBNode= pBVolumeMesh->getBNodeID(vBNodeID[ivert]);
            pBVolume->setBNode(ivert, pBNode);
        };
        if(pBVolume->getOrder()==ElementOrder::Second){
            uiint nNumOfVert = pBVolume->getNumOfVert();
            uiint nNumOfEdge = pBVolume->getNumOfEdge();
            uiint id;
            for(uiint iedge=0; iedge < nNumOfEdge; iedge++){
                id = vBNodeID[nNumOfVert + iedge];
                pBNode= pBVolumeMesh->getBNodeID(id);
                uiint ibnode = pBVolumeMesh->getBNodeIndex(id);
                pBVolumeMesh->setAggVol(ibnode, pBVolume->getID());
                pBVolume->setEdgeBNode(iedge, pBNode);
                pBVolume->markingEdge(iedge);         
                pBVolume->replaceEdgeBNode(iedge);
            };
        }
        pBVolume->calcVolume();
        if(bndType==BoundaryType::Neumann){
            double dCubicVol= pBVolume->getCubicVolume();
            double dBndValue= dCubicVol*val;
            pBVolume->setBndValue(dof, dBndValue);
        }
    }
}
void CMeshFactory::GeneBoundaryEdgeNode(const uiint& mgLevel, const uiint& bnd_id, const uiint& bndType,
        const uiint& mesh_id, const uiint& node_id, const uiint& b_node_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    CBoundaryEdgeMesh *pBEdgeMesh= mpTMesh->getBndEdgeMeshID(bnd_id);
    CNode *pNode= mpTMesh->getNode(node_id);
    CBoundaryNode *pBNode= new CBoundaryNode;
    pBNode->setNode(pNode);
    pBNode->setID(b_node_id);
    pBNode->setMGLevel(mgLevel);
    pBNode->resizeValue(mMGLevel+1);
    pBEdgeMesh->addBNode(pBNode);
}
void CMeshFactory::setValue_BoundaryEdgeNode(const uiint& mesh_id, const uiint& bnd_id, const uiint& bnode_id, const uiint& dof, const double& val)
{
    mpTAssyModel = mpGMGModel->getAssyModel(0);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    CBoundaryEdgeMesh *pBEdgeMesh= mpTMesh->getBndEdgeMeshID(bnd_id);
    CBoundaryNode *pBNode= pBEdgeMesh->getBNodeID(bnode_id);
    pBNode->addValue(dof, 0, val);
}
void CMeshFactory::initEdgeAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    CBoundaryEdgeMesh *pBEdgeMesh= mpTMesh->getBndEdgeMeshID(bnd_id);
    pBEdgeMesh->setupAggEdge();
}
void CMeshFactory::resizeEdgeAggregate(const uiint& mgLevel, const uiint& mesh_id, const uiint& bnd_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh_ID(mesh_id);
    CBoundaryEdgeMesh *pBEdgeMesh= mpTMesh->getBndEdgeMeshID(bnd_id);
    pBEdgeMesh->resizeAggEdge();
}
void CMeshFactory::GeneBoundaryEdge(const uiint& mgLevel, const uiint& bnd_id, const uiint& bndType, const uiint& elemType,
        const uiint& mesh_id, const uiint& elem_id, const uiint& edge_id, vuint& vBNodeID,
        const uiint& b_edge_id, const uiint& dof, const double& val)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    CBoundaryEdgeMesh *pBEdgeMesh= mpTMesh->getBndEdgeMeshID(bnd_id);
    pBEdgeMesh->setBndType(bndType);
    pBEdgeMesh->setID(bnd_id);
    uiint crIndex= pBEdgeMesh->getNumOfEdge();
    uiint crBEdgeID;
    CBoundaryEdge *pCrBEdge;
    uiint ib;
    bool bfind(false);
    if(crIndex > 0){
        if(crIndex==1){
            pCrBEdge= pBEdgeMesh->getBEdgeIX(0);
            crBEdgeID= pCrBEdge->getID();
            if(b_edge_id==crBEdgeID) bfind= true;
        }
        for(ib=crIndex-1; ib > 0; ib--){
            pCrBEdge= pBEdgeMesh->getBEdgeIX(ib);
            crBEdgeID= pCrBEdge->getID();
            if(b_edge_id==crBEdgeID) bfind= true;
        };
    }
    CBoundaryEdge *pBEdge;
    if(bfind){
        pBEdge= pBEdgeMesh->getBEdgeID(b_edge_id);
        if(bndType==BoundaryType::Neumann){
            double dLength= pBEdge->getLength();
            double dBndVal= val * dLength;
            pBEdge->setBndValue(dof, dBndVal);
        }
    }else{
        pBEdge = new CBoundaryEdge();
        pBEdge->setElementID(elem_id);
        pBEdge->setElementEdgeID(edge_id);
        pBEdge->setID(b_edge_id);
        pBEdge->setBEdgeShape(elemType);
        CElement *pElem= mpTMesh->getElement(elem_id);
        pBEdge->setElement(pElem);
        switch(elemType){
            case(ElementType::Beam):case(ElementType::Beam2):
                pBEdge->resizeBNode(2);
                break;
            default:
                break;
        }
        pBEdgeMesh->addBEdge(pBEdge);
    }
    if(!bfind){
        CBoundaryNode *pBNode;
        uiint numOfVert= pBEdge->getNumOfVert();
        uiint ivert;
        for(ivert=0; ivert < numOfVert; ivert++){
            pBNode= pBEdgeMesh->getBNodeID(vBNodeID[ivert]);
            pBEdge->setBNode(ivert, pBNode);
        };
        if(pBEdge->getOrder()==ElementOrder::Second){
            pBNode= pBEdgeMesh->getBNodeID(vBNodeID[2]);
            uiint ibnode = pBEdgeMesh->getBNodeIndex(vBNodeID[2]);
            pBEdgeMesh->setAggEdge(ibnode, pBEdge->getID());
            pBEdge->setEdgeBNode(pBNode);
            pBEdge->replaceEdgeBNode();
        }
        pBEdge->calcLength();
        if(bndType==BoundaryType::Neumann){
            double dLength= pBEdge->getLength();
            double dBndVal= val * dLength;
            pBEdge->setBndValue(dof, dBndVal);
        }
    }
}
void CMeshFactory::refineBoundary()
{
    uiint numOfMesh;
    uiint iLevel, iMesh;
    for(iLevel=0; iLevel < mMGLevel; iLevel++){
        mpTAssyModel = mpGMGModel->getAssyModel(iLevel);
        CAssyModel *pProgAssy= mpGMGModel->getAssyModel(iLevel+1);
        numOfMesh= mpTAssyModel->getNumOfMesh();
        for(iMesh=0; iMesh < numOfMesh; iMesh++){
            mpTMesh = mpTAssyModel->getMesh(iMesh);
            CMesh *pProgMesh= pProgAssy->getMesh(iMesh);
            uiint numOfBFaceMesh= mpTMesh->getNumOfBoundaryFaceMesh();
            pProgMesh->reserveBndFaceMesh(numOfBFaceMesh);
            uiint iBFaceMesh;
            for(iBFaceMesh=0; iBFaceMesh < numOfBFaceMesh; iBFaceMesh++){
                CBoundaryFaceMesh *pBFaceMesh= mpTMesh->getBndFaceMeshIX(iBFaceMesh);
                CBoundaryFaceMesh *pProgBFaceMesh= new CBoundaryFaceMesh;
                pProgMesh->setBndFaceMesh(pProgBFaceMesh);
                pProgBFaceMesh->resizeDOF(pBFaceMesh->getNumOfDOF());
                uiint idof, dof;
                for(idof=0; idof < pBFaceMesh->getNumOfDOF(); idof++){
                    dof = pBFaceMesh->getDOF(idof);
                    pProgBFaceMesh->setDOF(idof, dof);
                }
                pBFaceMesh->GeneEdgeBNode();
                pBFaceMesh->GeneFaceBNode();
                pBFaceMesh->refine(pProgBFaceMesh);
                pProgBFaceMesh->resizeAggFace();
                pProgBFaceMesh->setupAggFace();
                uiint nBndType= pBFaceMesh->getBndType();
                pProgBFaceMesh->setBndType(nBndType);
                uiint nID= pBFaceMesh->getID();
                pProgBFaceMesh->setID(nID);
                pProgBFaceMesh->setMGLevel(iLevel+1);
                pProgBFaceMesh->setMaxMGLevel(mMGLevel);
            };
            uiint numOfEdgeMesh= mpTMesh->getNumOfBoundaryEdgeMesh();
            pProgMesh->reserveBndEdgeMesh(numOfEdgeMesh);
            uiint iBEdgeMesh;
            for(iBEdgeMesh=0; iBEdgeMesh < numOfEdgeMesh; iBEdgeMesh++){
                CBoundaryEdgeMesh *pBEdgeMesh= mpTMesh->getBndEdgeMeshIX(iBEdgeMesh);
                CBoundaryEdgeMesh *pProgBEdgeMesh= new CBoundaryEdgeMesh;
                pProgMesh->setBndEdgeMesh(pProgBEdgeMesh);
                pProgBEdgeMesh->resizeDOF(pBEdgeMesh->getNumOfDOF());
                uiint idof, dof;
                for(idof=0; idof < pBEdgeMesh->getNumOfDOF(); idof++){
                    dof = pBEdgeMesh->getDOF(idof);
                    pProgBEdgeMesh->setDOF(idof, dof);
                }
                pBEdgeMesh->GeneEdgeBNode();
                pBEdgeMesh->refine(pProgBEdgeMesh);
                pProgBEdgeMesh->resizeAggEdge();
                pProgBEdgeMesh->setupAggEdge();
                uiint nBndType= pBEdgeMesh->getBndType();
                pProgBEdgeMesh->setBndType(nBndType);
                uiint nID= pBEdgeMesh->getID();
                pProgBEdgeMesh->setID(nID);
                pProgBEdgeMesh->setMGLevel(iLevel+1);
                pProgBEdgeMesh->setMaxMGLevel(mMGLevel);
            };
            uiint numOfVolMesh= mpTMesh->getNumOfBoundaryVolumeMesh();
            pProgMesh->reserveBndVolumeMesh(numOfVolMesh);
            uiint iBVolMesh;
            for(iBVolMesh=0; iBVolMesh < numOfVolMesh; iBVolMesh++){
                CBoundaryVolumeMesh *pBVolMesh= mpTMesh->getBndVolumeMeshIX(iBVolMesh);
                CBoundaryVolumeMesh *pProgBVolMesh= new CBoundaryVolumeMesh;
                pProgMesh->setBndVolumeMesh(pProgBVolMesh);
                pProgBVolMesh->resizeDOF(pBVolMesh->getNumOfDOF());
                uiint idof, dof;
                for(idof=0; idof < pBVolMesh->getNumOfDOF(); idof++){
                    dof = pBVolMesh->getDOF(idof);
                    pProgBVolMesh->setDOF(idof, dof);
                }
                pBVolMesh->GeneEdgeBNode();
                pBVolMesh->GeneFaceBNode();
                pBVolMesh->GeneVolBNode();
                pBVolMesh->refine(pProgBVolMesh);
                pProgBVolMesh->resizeAggVol();
                pProgBVolMesh->setupAggVol();
                uiint nBndType= pBVolMesh->getBndType();
                pProgBVolMesh->setBndType(nBndType);
                uiint nID= pBVolMesh->getID();
                pProgBVolMesh->setID(nID);
                pProgBVolMesh->setMGLevel(iLevel+1);
                pProgBVolMesh->setMaxMGLevel(mMGLevel);
            };
            CBNodeMeshGrp* pBNodeMeshGrp= mpTMesh->getBNodeMeshGrp();
            pProgMesh->setBNodeMeshGrp(pBNodeMeshGrp);
        };
    };
    mpTAssyModel = mpGMGModel->getAssyModel(mMGLevel);
    numOfMesh= mpTAssyModel->getNumOfMesh();
    for(iMesh=0; iMesh < numOfMesh; iMesh++){
        mpTMesh = mpTAssyModel->getMesh(iMesh);
        uiint numOfBFaceMesh= mpTMesh->getNumOfBoundaryFaceMesh();
        uiint iBFaceMesh;
        for(iBFaceMesh=0; iBFaceMesh < numOfBFaceMesh; iBFaceMesh++){
            CBoundaryFaceMesh *pBFaceMesh= mpTMesh->getBndFaceMeshIX(iBFaceMesh);
            pBFaceMesh->GeneEdgeBNode();
        };
        uiint numOfBEdgeMesh= mpTMesh->getNumOfBoundaryEdgeMesh();
        uiint iBEdgeMesh;
        for(iBEdgeMesh=0; iBEdgeMesh < numOfBEdgeMesh; iBEdgeMesh++){
            CBoundaryEdgeMesh *pBEdgeMesh= mpTMesh->getBndEdgeMeshIX(iBEdgeMesh);
            pBEdgeMesh->GeneEdgeBNode();
        };
        uiint numOfBVolMesh= mpTMesh->getNumOfBoundaryVolumeMesh();
        uiint iBVolMesh;
        for(iBVolMesh=0; iBVolMesh < numOfBVolMesh; iBVolMesh++){
            CBoundaryVolumeMesh *pBVolMesh= mpTMesh->getBndVolumeMeshIX(iBVolMesh);
            pBVolMesh->GeneEdgeBNode();
        };
    };
    for(iLevel=0; iLevel < mMGLevel+1; iLevel++){
        mpTAssyModel = mpGMGModel->getAssyModel(iLevel);
        numOfMesh= mpTAssyModel->getNumOfMesh();
        for(iMesh=0; iMesh < numOfMesh; iMesh++){
            mpTMesh = mpTAssyModel->getMesh(iMesh);
            uiint numOfBFaceMesh= mpTMesh->getNumOfBoundaryFaceMesh();
            uiint iBFaceMesh;
            for(iBFaceMesh=0; iBFaceMesh < numOfBFaceMesh; iBFaceMesh++){
                CBoundaryFaceMesh *pBFaceMesh= mpTMesh->getBndFaceMeshIX(iBFaceMesh);
                pBFaceMesh->distValueBNode();
            };
            uiint numOfBEdgeMesh= mpTMesh->getNumOfBoundaryEdgeMesh();
            uiint iBEdgeMesh;
            for(iBEdgeMesh=0; iBEdgeMesh < numOfBEdgeMesh; iBEdgeMesh++){
                CBoundaryEdgeMesh *pBEdgeMesh= mpTMesh->getBndEdgeMeshIX(iBEdgeMesh);
                pBEdgeMesh->distValueBNode();
            };
            uiint numOfBVolMesh= mpTMesh->getNumOfBoundaryVolumeMesh();
            uiint iBVolMesh;
            for(iBVolMesh=0; iBVolMesh < numOfBVolMesh; iBVolMesh++){
                CBoundaryVolumeMesh *pBVolMesh= mpTMesh->getBndVolumeMeshIX(iBVolMesh);
                pBVolMesh->distValueBNode();
            };
        };
    };
}
void CMeshFactory::setupBucketNode(const uiint& mgLevel, const uiint& mesh_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->setupBucketNode();
}
void CMeshFactory::initBucketNode(const uiint& mgLevel, const uiint& mesh_id, const uiint& maxID, const uiint& minID)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->initBucketNode(maxID, minID);
}
void CMeshFactory::setIDBucketNode(const uiint& mgLevel, const uiint& mesh_id, const uiint& id, const uiint& index)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->setupBucketNodeIndex(id, index);
}
void CMeshFactory::setupBucketElement(const uiint& mgLevel, const uiint& mesh_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->setupBucketElement();
}
void CMeshFactory::initBucketElement(const uiint& mgLevel, const uiint& mesh_id, const uiint& maxID, const uiint& minID)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    ;
    mpTMesh->initBucketElement(maxID, minID);
}
void CMeshFactory::setIDBucketElement(const uiint& mgLevel, const uiint& mesh_id, const uiint& id, const uiint& index)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->setupBucketElementIndex(id, index);
}
void CMeshFactory::reserveMaterial(const uiint& res_size)
{
     mpGMGModel->reserveMaterial(res_size);
}
void CMeshFactory::GeneMaterial(const uiint& mesh_id, const uiint& material_id, string& name, vuint& vType, vdouble& vValue)
{
    CMaterial *pMaterial = new CMaterial;
    pMaterial->setID(material_id);
    pMaterial->setName(name);
    pMaterial->setMeshID(mesh_id);
    for(uiint i=0; i< vType.size(); i++) pMaterial->setValue(vType[i],vValue[i]);
    mpGMGModel->setMaterial(pMaterial);
}
void CMeshFactory::reserveCommMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& res_size)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->reserveCommMesh(res_size);
}
void CMeshFactory::GeneCommMesh(const uiint& mgLevel, const uiint& mesh_id, const uiint& comID, const uiint& myRank, const uiint& nTransmitRank)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    CIndexBucket *pBucket= mpTMesh->getBucket();
    mpTCommMesh = new CCommMesh(pBucket);
    mpTCommMesh->setCommID(comID);
    mpTCommMesh->setRankID(myRank);
    mpTCommMesh->setTransmitRankID(nTransmitRank);
    mpTMesh->setCommMesh(mpTCommMesh);
}
void CMeshFactory::reserveCommNode(const uiint& mgLevel, const uiint& mesh_id, const uiint& commesh_id, const uiint& res_size)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh(mesh_id);
    mpTCommMesh= mpTMesh->getCommMesh(commesh_id);
    mpTCommMesh->reserveNode(res_size);   
    mpTCommMesh->resizeNodeRank(res_size);
}
void CMeshFactory::GeneCommNode(const uiint& mgLevel, const uiint& commNodeID,
                                     const uiint& mesh_id, const uiint& commesh_id, const uiint& nodeID, const uiint& rank)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh(mesh_id);
    mpTCommMesh= mpTMesh->getCommMesh(commesh_id);
    CNode* pNode= mpTMesh->getNode(nodeID);
    mpTCommMesh->setNode(pNode);
    mpTCommMesh->setNodeRank(commNodeID, rank);
    if(mpTCommMesh->getRankID()==rank) mpTCommMesh->setSendNode(pNode, commNodeID);
    if(mpTCommMesh->getTransmitRankID()==rank) mpTCommMesh->setRecvNode(pNode, commNodeID);
}
void CMeshFactory::reserveCommElement(const uiint& mgLevel, const uiint& mesh_id, const uiint& commesh_id, const uiint& res_size)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh(mesh_id);
    mpTCommMesh= mpTMesh->getCommMesh(commesh_id);
    mpTCommMesh->reserveCommElementAll(res_size);
}
void CMeshFactory::GeneCommElement(const uiint& mgLevel, const uiint& mesh_id, const uiint& commesh_id,
                                   const uiint& nType, const uiint& elemID, vuint& vCommNodeID)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh(mesh_id);
    mpTCommMesh= mpTMesh->getCommMesh(commesh_id);
    uiint numOfVert= vCommNodeID.size();
    vuint vNodeRank;  vNodeRank.reserve(numOfVert);
    uiint ivert, rank;
    uiint commNodeID;
    for(ivert=0; ivert< numOfVert; ivert++){
        commNodeID= vCommNodeID[ivert];
        rank= mpTCommMesh->getNodeRank(commNodeID);
        vNodeRank.push_back(rank);
    };
    CCommElement *pCommElem;
    switch(nType){
        case(ElementType::Hexa):
            pCommElem = new CCommHexa;
            break;
        case(ElementType::Tetra):
            pCommElem = new CCommTetra;
            break;
        case(ElementType::Prism):
            pCommElem = new CCommPrism;
            break;
        case(ElementType::Quad):
            pCommElem = new CCommQuad;
            break;
        case(ElementType::Triangle):
            pCommElem = new CCommTriangle;
            break;
        case(ElementType::Beam):
            pCommElem = new CCommBeam;
            break;
    }
    CElement* pElem= mpTMesh->getElement(elemID);
    pCommElem->setElement(pElem);
    for(ivert=0; ivert< numOfVert; ivert++){
        rank = vNodeRank[ivert];
        pCommElem->setNodeRank(ivert, rank);
    };
    mpTCommMesh->setCommElementAll(pCommElem);
}
void CMeshFactory::GeneProgCommElem(CCommElement* pCommElem, vector<CCommElement*>& vProgCommElem)
{
    CCommElement *pProgCommElem;
    uiint ivert;
    switch(pCommElem->getShapeType()){
        case(ElementType::Hexa):
            vProgCommElem.reserve(8);
            for(ivert=0; ivert< 8; ivert++){
                pProgCommElem= new CCommHexa;
                vProgCommElem.push_back(pProgCommElem);
            };
            dividCommElem(pCommElem, vProgCommElem);
            mpLogger->Info(Utility::LoggerMode::Debug,"CommElem(Hexa)の分割");
            break;
        case(ElementType::Tetra):
            vProgCommElem.reserve(4);
            for(ivert=0; ivert< 4; ivert++){
                pProgCommElem= new CCommHexa;
                vProgCommElem.push_back(pProgCommElem);
            };
            dividCommElem(pCommElem, vProgCommElem);
            break;
        case(ElementType::Prism):
            vProgCommElem.reserve(6);
            for(ivert=0; ivert< 6; ivert++){
                pProgCommElem= new CCommHexa;
                vProgCommElem.push_back(pProgCommElem);
            };
            dividCommElem(pCommElem, vProgCommElem);
            break;
        case(ElementType::Quad):
            vProgCommElem.reserve(4);
            for(ivert=0; ivert< 4; ivert++){
                pProgCommElem= new CCommQuad;
                vProgCommElem.push_back(pProgCommElem);
            };
            dividCommElem(pCommElem, vProgCommElem);
            break;
        case(ElementType::Triangle):
            vProgCommElem.reserve(3);
            for(ivert=0; ivert< 3; ivert++){
                pProgCommElem= new CCommQuad;
                vProgCommElem.push_back(pProgCommElem);
            };
            dividCommElem(pCommElem, vProgCommElem);
            break;
        case(ElementType::Beam):
            vProgCommElem.reserve(2);
            for(ivert=0; ivert< 2; ivert++){
                pProgCommElem= new CCommBeam;
                vProgCommElem.push_back(pProgCommElem);
            };
            dividCommElem(pCommElem, vProgCommElem);
            break;
        default:
            mpLogger->Info(Utility::LoggerMode::Error, "ShapeType Error @MeshFactory::GeneProgCommElem");
            break;
    }
}
void CMeshFactory::dividCommElem(CCommElement* pCommElem, vector<CCommElement*>& vProgCommElem)
{
    CProgElementTree *pProgTree = CProgElementTree::Instance();
    CElement* pElem= pCommElem->getElement();
    CElement* pProgElem;
    CCommElement* pProgCommElem;
    uiint nRank;
    uiint ivert, progvert;
    uiint iedge, iface;
    uiint invalid= pProgTree->getInvalidNum();
    uiint numOfVert, numOfEdge, numOfFace;
    numOfVert= pElem->getNumOfNode(); numOfEdge= pElem->getNumOfEdge(); numOfFace= pElem->getNumOfFace();
    for(ivert=0; ivert< numOfVert; ivert++){
        pProgElem= pElem->getProgElem(ivert);
        pProgCommElem= vProgCommElem[ivert];
        pProgCommElem->setElement(pProgElem);
        progvert= pProgTree->getVertProgVert(ivert, pCommElem->getShapeType());
        nRank= pCommElem->getNodeRank(ivert);
        pProgCommElem->setNodeRank(progvert,nRank);
        for(iedge=0; iedge< numOfEdge; iedge++){
            progvert= pProgTree->getEdgeProgVert(iedge, ivert, pCommElem->getShapeType());
            if(progvert != invalid){
                nRank= pCommElem->getEdgeRank(iedge);
                pProgCommElem->setNodeRank(progvert, nRank);
            }
        };
        for(iface=0; iface< numOfFace; iface++){
            progvert= pProgTree->getFaceProgVert(iface, ivert, pCommElem->getShapeType());
            if(progvert != invalid){
                nRank= pCommElem->getFaceRank(iface);
                pProgCommElem->setNodeRank(progvert, nRank);
            }
        };
        progvert= pProgTree->getVolProgVert(ivert, pCommElem->getShapeType());
        nRank= pCommElem->getVolRank();
        pProgCommElem->setNodeRank(progvert, nRank);
    };
}
void CMeshFactory::GeneContactMesh(const uiint& contactID, const uiint& myRank, const uiint& transRank, const uiint& nProp)
{
    uiint ilevel=0;
        mpTAssyModel= mpGMGModel->getAssyModel(ilevel);
        CContactMesh *pContactMesh= new CContactMesh;
        pContactMesh->setID(contactID);
        pContactMesh->setLevel(ilevel);
        pContactMesh->setRank(myRank);
        pContactMesh->setTransmitRank(transRank);
        pContactMesh->setProp(nProp);
        mpTAssyModel->addContactMesh(pContactMesh, contactID);
}
void CMeshFactory::GeneContactNode(const uiint& mgLevel, const uiint& contactID, const uiint& conNodeID, const vdouble& vCoord,
        const string& s_param_type, const uiint& numOfVector, const uiint& numOfScalar,
        bool bmesh, const uiint& meshID, const uiint& nodeID,
        const uiint& rank, const uiint& maslave)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    CContactMesh *pConMesh= mpTAssyModel->getContactMesh_ID(contactID);
    CContactNode *pConNode= new CContactNode;
    pConNode->setLevel(mgLevel);
    pConNode->pushLevelMarking();
    pConNode->setID(conNodeID);
    pConNode->setCoord(vCoord);
    if(bmesh){ pConNode->setMeshID(meshID); pConNode->markingSelfMesh();}
    if(bmesh){ pConNode->setNodeID(nodeID); pConNode->markingSelfNode();}
    pConNode->setRank(rank);
    if(s_param_type=="v" || s_param_type=="V" || s_param_type=="sv" || s_param_type=="SV"){
        pConNode->resizeDisp(numOfVector);
        pConNode->initDisp();
    }
    if(s_param_type=="s" || s_param_type=="S" || s_param_type=="sv" || s_param_type=="SV"){
        pConNode->resizeScalar(numOfVector);
        pConNode->initScalar();
    }
    pConMesh->addConNode(pConNode, conNodeID);
    switch(maslave){
    case(0):
        pConMesh->addMasterConNode(pConNode,conNodeID);
        break;
    case(1):
        pConMesh->addSlaveConNode(pConNode,conNodeID);
        break;
    default:
        break;
    }
}
void CMeshFactory::GeneMasterFace(const uiint& contactID, const uiint& shapeType, const uiint& masterFaceID,
        bool bmesh, const uiint& meshID, const uiint& elemID, const uiint& elemFaceID,
        const vuint& vConNodeID, const uiint& face_rank)
{
    uiint mgLevel(0);
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    CSkinFace *pMFace= new CMasterFace;
    pMFace->setShapeType(shapeType);
    if(bmesh){
        CMesh *pMMesh;
        pMMesh= mpTAssyModel->getMesh_ID(meshID);
        CElement *pMasterElem;
        pMasterElem= pMMesh->getElement(elemID);
        pMasterElem->markingMPCMaster();
        pMasterElem->markingMPCFace(elemFaceID);
        pMFace->markingSelf();
        pMFace->setMeshID(meshID);
        pMFace->setElementID(elemID);
        pMFace->setFaceID(elemFaceID);
    }
    CContactMesh *pConMesh= mpTAssyModel->getContactMesh_ID(contactID);
    pMFace->setID(masterFaceID);
    pMFace->setRank(face_rank);
    pMFace->setLevel(mgLevel);
    uiint icnode, numOfConNode = vConNodeID.size();
    for(icnode=0; icnode< numOfConNode; icnode++){
        CContactNode *pConNode= pConMesh->getContactNode_ID(vConNodeID[icnode]);
        pMFace->addNode(pConNode);
    };
    if(pMFace->getOrder()==ElementOrder::Second){
        uiint nNumOfEdge= pMFace->getNumOfEdge();
        uiint nNumOfVert= pMFace->getNumOfVert();
        for(uiint iedge=0; iedge < nNumOfEdge; iedge++){
            CContactNode *pConNode= pMFace->getNode(nNumOfVert + iedge);
            pMFace->setEdgeConNode(pConNode, iedge);
            pMFace->markingEdgeNode(iedge);
        };
    }
    pConMesh->addMasterFace(pMFace);
}
void CMeshFactory::GeneSlaveFace(const uiint& contactID, const uiint& shapeType, const uiint& slaveFaceID,
        bool bmesh, const uiint& meshID, const uiint& elemID, const uiint& elemFaceID,
        const vuint& vConNodeID, const uiint& face_rank)
{
    uiint mgLevel(0);
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    CSkinFace *pSFace= new CSkinFace;
    pSFace->setShapeType(shapeType);
    if(bmesh){
        CMesh *pSMesh;
        pSMesh= mpTAssyModel->getMesh_ID(meshID);
        CElement *pSlaveElem;
        pSlaveElem= pSMesh->getElement(elemID);
        pSlaveElem->markingMPCSlave();
        pSlaveElem->markingMPCFace(elemFaceID);
        pSFace->markingSelf();
        pSFace->setMeshID(meshID);
        pSFace->setElementID(elemID);
        pSFace->setFaceID(elemFaceID);
    }
    CContactMesh *pConMesh= mpTAssyModel->getContactMesh_ID(contactID);
    pSFace->setID(slaveFaceID);
    pSFace->setRank(face_rank);
    pSFace->setLevel(mgLevel);
    uiint icnode, numOfConNode = vConNodeID.size();
    for(icnode=0; icnode< numOfConNode; icnode++){
        CContactNode *pConNode= pConMesh->getContactNode_ID(vConNodeID[icnode]);
        pSFace->addNode(pConNode);
    };
    if(pSFace->getOrder()==ElementOrder::Second){
        uiint nNumOfEdge= pSFace->getNumOfEdge();
        uiint nNumOfVert= pSFace->getNumOfVert();
        for(uiint iedge=0; iedge < nNumOfEdge; iedge++){
            CContactNode *pConNode= pSFace->getNode(nNumOfVert + iedge);
            pSFace->setEdgeConNode(pConNode, iedge);
            pSFace->markingEdgeNode(iedge);
        };
    }
    pConMesh->addSlaveFace(pSFace);
}
void CMeshFactory::refineContactMesh()
{
    CAssyModel  *pAssy,*pProgAssy;
    CContactMesh *pConMesh,*pProgConMesh;
    CSkinFace  *pSkinFace;
    vector<CSkinFace*> vProgFace;
    uiint maslave;                
    uiint meshID,elemID;
    CMesh *pMesh;
    CElement *pElem;
    uiint faceID;
    uiint maxLayer;
    uiint ilevel;
    for(ilevel=0; ilevel< mMGLevel; ilevel++){
        pAssy =  mpGMGModel->getAssyModel(ilevel);    
        pProgAssy= mpGMGModel->getAssyModel(ilevel+1);
        uiint numOfCont= pAssy->getNumOfContactMesh();
        uiint icont;
        for(icont=0; icont< numOfCont; icont++){
            pConMesh= pAssy->getContactMesh(icont);
            pProgConMesh= pProgAssy->getContactMesh(icont);
            pConMesh->setupCoarseConNode(pProgConMesh);      
            pConMesh->setupAggSkinFace();                    
            pConMesh->setupEdgeConNode(pProgConMesh, ilevel);
            pConMesh->setupFaceConNode(pProgConMesh);        
            uiint numOfSkinFace;
            uiint iface;
            for(maslave=0; maslave< 2; maslave++){
                faceID=0;
                if(maslave==0) numOfSkinFace= pConMesh->getNumOfMasterFace();
                if(maslave==1) numOfSkinFace= pConMesh->getNumOfSlaveFace(); 
                for(iface=0; iface< numOfSkinFace; iface++){
                    if(maslave==0)  pSkinFace= pConMesh->getMasterFace(iface);
                    if(maslave==1)  pSkinFace= pConMesh->getSlaveFace(iface); 
                    if(pSkinFace->isSelf() && pSkinFace->getNumOfEdge()!=0){
                        meshID= pSkinFace->getMeshID(); elemID= pSkinFace->getElementID();
                        pMesh= pAssy->getMesh_ID(meshID);
                        pElem= pMesh->getElement(elemID);
                        pSkinFace->refine(pElem, faceID);
                    }else{
                        pSkinFace->refine(NULL, faceID); 
                    }
                    vProgFace= pSkinFace->getProgFace();
                    if(maslave==0) pProgConMesh->addMasterFace(vProgFace);
                    if(maslave==1) pProgConMesh->addSlaveFace(vProgFace); 
                };
            };
        };
    };
    pAssy =  mpGMGModel->getAssyModel(mMGLevel);
    uiint numOfCont= pAssy->getNumOfContactMesh();
    uiint icont;
    for(icont=0; icont < numOfCont; icont++){
        pConMesh= pAssy->getContactMesh(icont);
        pConMesh->setupAggSkinFace();              
        pConMesh->setupEdgeConNode(NULL, mMGLevel);
    };
    for(icont=0; icont < numOfCont; icont++){
        pConMesh= pAssy->getContactMesh(icont);
        for(maslave=0; maslave < 2; maslave++){
            uiint numOfFace;
            if(maslave==0) numOfFace = pConMesh->getNumOfMasterFace();
            if(maslave==1) numOfFace = pConMesh->getNumOfSlaveFace();
            uiint iface;
            for(iface=0; iface< numOfFace; iface++){
                if(maslave==0)  pSkinFace= pConMesh->getMasterFace(iface);
                if(maslave==1)  pSkinFace= pConMesh->getSlaveFace(iface); 
                if(pSkinFace->getOrder()==ElementOrder::Second){
                    if(pSkinFace->isSelf() && pSkinFace->getNumOfEdge()!=0){
                        meshID= pSkinFace->getMeshID(); elemID= pSkinFace->getElementID();
                        pMesh= pAssy->getMesh_ID(meshID);
                        pElem= pMesh->getElement(elemID);
                        pSkinFace->setupNodeID_2nd_LastLevel(pElem);
                    }
                }
            };
        };
    };
    for(ilevel=0; ilevel < mMGLevel+1; ilevel++){
        pAssy= mpGMGModel->getAssyModel(ilevel);
        uiint numOfCont= pAssy->getNumOfContactMesh();
        uiint icont;
        for(icont=0; icont < numOfCont; icont++){
            pConMesh= pAssy->getContactMesh(icont);
            uiint nRange;
            nRange= pConMesh->getNumOfConNode();
            uiint nDigitCount(0);
            while(nRange > 100){
                nRange /= 10;
                nDigitCount++;
            };
            maxLayer= nDigitCount;
            pConMesh->generateOctree(maxLayer);
        };
    };
}
void CMeshFactory::refineCommMesh2()
{
    CAssyModel *pAssy,*pProgAssy;
    CMesh *pMesh,*pProgMesh;
    CCommMesh2 *pCommMesh2,*pProgCommMesh2;
    CCommFace *pCommFace;
    vector<CCommFace*> mvCommFace;
    CCommFace *pProgCommFace;
    uiint countID(0);
    uiint ilevel;
    for(ilevel=0; ilevel< mMGLevel; ilevel++){
        pAssy= mpGMGModel->getAssyModel(ilevel);
        pProgAssy= mpGMGModel->getAssyModel(ilevel+1);
        uiint imesh, numOfMesh;
        numOfMesh= pAssy->getNumOfMesh();
        for(imesh=0; imesh< numOfMesh; imesh++){
            pMesh= pAssy->getMesh(imesh);
            pProgMesh= pProgAssy->getMesh(imesh);
            uiint icomm, numOfComm;
            numOfComm= pMesh->getCommMesh2Size();
            for(icomm=0; icomm< numOfComm; icomm++){
                pCommMesh2= pMesh->getCommMesh2IX(icomm);
                pProgCommMesh2 = new CCommMesh2;
                pProgCommMesh2->setLevel(ilevel+1);
                pProgCommMesh2->setID(pCommMesh2->getID());
                pProgCommMesh2->setRank(pCommMesh2->getRank());
                pProgCommMesh2->setTransmitRank(pCommMesh2->getTrasmitRank());
                pProgMesh->setCommMesh2(pProgCommMesh2);
                pCommMesh2->setupCommNode(pProgCommMesh2);
                pCommMesh2->setupAggFace();
                pCommMesh2->setupEdgeCommNode(pProgCommMesh2, ilevel);
                pCommMesh2->setupFaceCommNode(pProgCommMesh2);
                uiint iface, numOfFace;
                numOfFace= pCommMesh2->getCommFaceSize();
                for(iface=0; iface< numOfFace; iface++){
                    pCommFace= pCommMesh2->getCommFaceIX(iface);
                    uiint elemID = pCommFace->getElementID();
                    CElement *pElem= pMesh->getElement(elemID);
                    mvCommFace= pCommFace->refine(pElem);
                    uiint ipface,numOfProgFace;
                    numOfProgFace= mvCommFace.size();
                    for(ipface=0; ipface< numOfProgFace; ipface++){
                        pProgCommFace= mvCommFace[ipface];
                        pProgCommFace->setID(countID);
                        pProgCommMesh2->addCommFace(pProgCommFace);
                        countID++;
                    };
                };
                uiint elemID,entity_num;
                CElement  *pElem;
                CNode     *pFaceNode;
                CCommNode *pFaceCommNode;
                for(iface=0; iface< numOfFace; iface++){
                    pCommFace= pCommMesh2->getCommFaceIX(iface);
                    if(pCommFace->getNumOfEdge() > 2){
                        pFaceCommNode= pCommFace->getFaceCommNode();
                        elemID= pCommFace->getElementID();
                        entity_num= pCommFace->getElementFaceID();
                        pElem= pMesh->getElement(elemID);
                        pFaceNode= pElem->getFaceNode(entity_num);
                        pFaceCommNode->setNode(pFaceNode);
                    }
                    PairCommNode pairCommNode;
                    CCommNode *pEdgeCommNode;
                    CNode *pNodeFir, *pNodeSec;
                    CNode *pEdgeNode;
                    uiint iedge, numOfEdge;
                    numOfEdge= pCommFace->getNumOfEdge();
                    for(iedge=0; iedge< numOfEdge; iedge++){
                        pairCommNode= pCommFace->getEdgePairCommNode(iedge);
                        pNodeFir= pairCommNode.first->getNode();
                        pNodeSec= pairCommNode.second->getNode();
                        uiint edgeIndex;
                        edgeIndex= pElem->getEdgeIndex(pNodeFir,pNodeSec);
                        pEdgeNode= pElem->getEdgeInterNode(edgeIndex);
                        pEdgeCommNode= pCommFace->getEdgeCommNode(iedge);
                        pEdgeCommNode->setNode(pEdgeNode);
                    };
                };
            };
        };
    };
    pAssy= mpGMGModel->getAssyModel(mMGLevel);
    uiint imesh;
    uiint nNumOfMesh= pAssy->getNumOfMesh();
    for(imesh=0; imesh < nNumOfMesh; imesh++){
        pMesh= pAssy->getMesh(imesh);
        uiint icomm;
        uiint nNumOfComm= pMesh->getCommMesh2Size();
        for(icomm=0; icomm< nNumOfComm; icomm++){
            pCommMesh2= pMesh->getCommMesh2IX(icomm);
            pCommMesh2->setupAggFace();
            pCommMesh2->setupEdgeCommNode(NULL, ilevel);
            uiint iface;
            uiint nNumOfFace = pCommMesh2->getCommFaceSize();
            for(iface=0; iface< nNumOfFace; iface++){
                pCommFace= pCommMesh2->getCommFaceIX(iface);
                uiint elemID = pCommFace->getElementID();
                CElement *pElem= pMesh->getElement(elemID);
                PairCommNode pairCommNode;
                CCommNode *pEdgeCommNode;
                CNode *pNodeFir, *pNodeSec;
                CNode *pEdgeNode;
                uiint iedge, numOfEdge;
                numOfEdge= pCommFace->getNumOfEdge();
                for(iedge=0; iedge< numOfEdge; iedge++){
                    pairCommNode= pCommFace->getEdgePairCommNode(iedge);
                    pNodeFir= pairCommNode.first->getNode();
                    pNodeSec= pairCommNode.second->getNode();
                    uiint edgeIndex;
                    edgeIndex= pElem->getEdgeIndex(pNodeFir,pNodeSec);
                    pEdgeNode= pElem->getEdgeInterNode(edgeIndex);
                    pEdgeCommNode= pCommFace->getEdgeCommNode(iedge);
                    pEdgeCommNode->setNode(pEdgeNode);
                };
            };
        };
    };
}
void CMeshFactory::GeneCommMesh2(const uiint& mgLevel, const uiint& mesh_id, const uiint& comID,
        const uiint& numOfFace, const uiint& numOfCommNode,
        const uiint& myRank, const uiint& nTransmitRank)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh(mesh_id);
    mpTCommMesh2= new CCommMesh2;
    mpTCommMesh2->setLevel(mgLevel);
    mpTCommMesh2->setID(comID);
    mpTCommMesh2->reserveCommFace(numOfFace);
    mpTCommMesh2->reserveCommNode(numOfCommNode);
    mpTCommMesh2->setRank(myRank);
    mpTCommMesh2->setTransmitRank(nTransmitRank);
    mpTMesh->setCommMesh2(mpTCommMesh2);
}
void CMeshFactory::GeneCommFace(const uiint& mgLevel, const uiint& commeshID, const uiint& face_id,
            const uiint& mesh_id,const uiint elem_id, const uiint& elem_ent_num, const uiint& elem_type, const vuint& vCommNodeID)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh_ID(mesh_id);
    CElement *pElement= mpTMesh->getElement(elem_id);
    pElement->markingCommMesh2();
    pElement->markingCommEntity(elem_ent_num);
    mpTCommMesh2= mpTMesh->getCommMesh2(commeshID);
    CCommFace *pCommFace= new CCommFace;
    pCommFace->setID(face_id);
    pCommFace->setElementID(elem_id);
    pCommFace->setElementFaceID(elem_ent_num);
    pCommFace->setMGLevel(mgLevel);
    uiint nNumOfVert, nNumOfEdge, nOrder;
    switch(elem_type)
    {
        case(ElementType::Quad):
            nNumOfVert= 4; 
            nNumOfEdge= 4;
            nOrder = ElementOrder::First;
            break;
        case(ElementType::Quad2):
            nNumOfVert= 4;
            nNumOfEdge= 4;
            nOrder = ElementOrder::Second;
            break;
        case(ElementType::Triangle):
            nNumOfVert= 3;
            nNumOfEdge= 3;
            nOrder = ElementOrder::First;
            break;
        case(ElementType::Triangle2):
            nNumOfVert= 3;
            nNumOfEdge= 3;
            nOrder = ElementOrder::Second;
            break;
        case(ElementType::Beam):
            nNumOfVert= 2;
            nNumOfEdge= 1;
            nOrder = ElementOrder::First;
            break;
        case(ElementType::Beam2):
            nNumOfVert= 2;
            nNumOfEdge= 1;
            nOrder = ElementOrder::Second;
            break;
        default:
            break;
    }
    pCommFace->initialize(nNumOfVert, nNumOfEdge, nOrder);
    uiint nNumOfNode = vCommNodeID.size();
    uiint i, id;
    CCommNode *pCommNode;
    for(i=0; i< nNumOfNode; i++){
        id= vCommNodeID[i];
        pCommNode= mpTCommMesh2->getCommNode(id);
        pCommFace->setCommNode(i, pCommNode);
    };
    bool bSecond(false);
    if(elem_type==ElementType::Quad2)     bSecond=true;
    if(elem_type==ElementType::Triangle2) bSecond=true;
    if(elem_type==ElementType::Beam2)     bSecond=true;
    if(bSecond){
        for(uiint iedge=0; iedge< nNumOfEdge; iedge++){
            id= vCommNodeID[nNumOfVert + iedge];
            pCommNode= mpTCommMesh2->getCommNode(id);
            pCommFace->setEdgeCommNode(pCommNode, iedge);
        };
    }
    mpTCommMesh2->addCommFace(pCommFace);
}
void CMeshFactory::GeneCommNodeCM2(const uiint& mgLevel, const uiint& mesh_id, const uiint& node_id,const uiint& commeshID,
        const uiint& comm_node_id, const vdouble& vCoord)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh_ID(mesh_id);
    mpTCommMesh2= mpTMesh->getCommMesh2(commeshID);
    CNode* pNode= mpTMesh->getNode(node_id);
    CCommNode *pCommNode= new CCommNode;
    pCommNode->setID(comm_node_id);
    pCommNode->setNode(pNode);
    pCommNode->setCoord(vCoord);
    mpTCommMesh2->addCommNode(pCommNode);
}
void CMeshFactory::GeneElemGrpOBJ(const uiint& mgLevel, const uiint& mesh_id, const vuint& vGrpID, vstring& vGrpName)
{
    CAssyModel *pAssyModel = mpGMGModel->getAssyModel(mgLevel);
    CMesh *pMesh = pAssyModel->getMesh_ID(mesh_id);
    uiint i, nNumOfElemGrp = vGrpID.size();
    for(i=0; i < nNumOfElemGrp; i++){
        CElementGroup *pElemGrp = new CElementGroup;
        pElemGrp->setMesh(pMesh);
        uiint nGrpID = vGrpID[i];
        pElemGrp->setID(nGrpID);
        string sGrpName = vGrpName[i];
        pElemGrp->setName(sGrpName);
        pMesh->addElemGrp(pElemGrp);
    }
}
void CMeshFactory::setElemID_with_ElemGrp(const uiint& mgLevel, const uiint& mesh_id, const uiint& nGrpID, const vuint& vElemID)
{
    CAssyModel *pAssyModel = mpGMGModel->getAssyModel(mgLevel);
    CMesh *pMesh = pAssyModel->getMesh_ID(mesh_id);
    CElementGroup *pElemGrp = pMesh->getElemGrpID(nGrpID);
    uiint i, nNumOfElem=vElemID.size();
    for(i=0; i < nNumOfElem; i++){
        pElemGrp->addElementID(vElemID[i]);
    }
}
