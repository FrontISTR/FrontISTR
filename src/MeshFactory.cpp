/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   MeshFactory.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "CommNode.h"
#include "CommMesh2.h"
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
#include "SkinFace.h"
using namespace pmw;
CMeshFactory::CMeshFactory(void)
{
    mpLogger = Utility::CLogger::Instance();
}
CMeshFactory::~CMeshFactory(void)
{
    cout << "~CMeshFactory" << endl;
}
void CMeshFactory::refineMesh()
{
    CAssyModel *pAssy, *pProgAssy;
    CMesh      *pMesh, *pProgMesh;
    CElement     *pElem=NULL;
    vector<CElement*> vProgElem;
    CElement     *pProgElem=NULL;
    CCommMesh  *pCommMesh, *pProgCommMesh;
    CCommElement      *pCommElem;
    vector<CCommElement*> vProgCommElem;
    uint numOfCommMesh,numOfCommElemAll,numOfCommNode;
    uint icommesh,icomelem,iprocom;
    pAssy= mpGMGModel->getAssyModel(0);
    uint numOfMesh= pAssy->getNumOfMesh();
    uint ilevel,imesh,ielem;
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
                pMesh->setupAggregate(); 
                mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->setupAggElement finish at ilevel==0");
            }
            pMesh->presetProgMesh(pProgMesh);
            mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->presetProgMesh finish");
            pMesh->setupEdgeElement(pProgMesh);
            mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->setupEdgeElement finish");
            pMesh->setupFaceElement(pProgMesh);
            mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->setupFaceElement finish");
            pMesh->setupVolumeNode(pProgMesh); 
            mpLogger->Info(Utility::LoggerMode::MWDebug,"pMesh->setupVolumeNode finish");
            uint numOfElem = pMesh->getNumOfElement();
            uint elementID= 0;
            for(ielem=0; ielem< numOfElem; ielem++){
                pElem= pMesh->getElement(ielem);
                vProgElem.clear();
                /*GeneProgElem(ilevel, pElem, pProgElem, vProgElem, elementID, pProgMesh);*/
				GeneProgElem(ilevel, pElem, vProgElem, elementID, pProgMesh);
                uint i, nBaseType;
                for(i=0; i< vProgElem.size(); i++){
                    pProgElem= vProgElem[i];
                    nBaseType= pProgElem->getEntityType();
                    if(nBaseType==BaseElementType::Solid || nBaseType==BaseElementType::Shell)
                        pProgElem->setupFaceCnvNodes();
                };
            };
            pProgMesh->setupNumOfNode();
            pProgMesh->setupNumOfElement();
            uint numOfNode= pProgMesh->getNodeSize();
            CAggregateElement *pAggElem;
            CAggregateNode    *pAggNode;
            pProgMesh->reserveAggregate(numOfNode);
            for(uint iAgg=0; iAgg< numOfNode; iAgg++){
                pAggElem = new CAggregateElement;
                pAggNode = new CAggregateNode;
                pProgMesh->setAggElement(pAggElem);
                pProgMesh->setAggNode(pAggNode);
            };
            pProgAssy->setMesh(pProgMesh,pMesh->getMeshID());
            uint progNodeSize = pProgMesh->getNodeSize();
            pProgMesh->initBucketNode(progNodeSize+1, 0);
            pProgMesh->setupBucketNode();                
            uint progElemSize = pProgMesh->getElementSize();
            pProgMesh->initBucketElement(progElemSize+1,0);
            pProgMesh->setupBucketElement();
            pProgMesh->setupAggregate();
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
        pMesh->setupAggregate();
        pMesh->setupEdgeElement(NULL);
    };
}
void CMeshFactory::GeneProgElem(const uint& ilevel,CElement* pElem, vector<CElement*>& vProgElem, uint& elementID, CMesh* pProgMesh)
{
	CElement* pProgElem;
    uint i;
    switch(pElem->getType()){
        case(ElementType::Hexa):case(ElementType::Hexa2):
            vProgElem.reserve(8);
            for(i=0; i< 8; i++){
                pProgElem= new CHexa; pProgElem->setMGLevel(ilevel+1);
                vProgElem.push_back(pProgElem);
            };
            dividHexa(pElem,vProgElem, elementID, pProgMesh);
            break;
        case(ElementType::Tetra):case(ElementType::Tetra2):
            vProgElem.reserve(4);
            for(i=0; i< 4; i++){
                pProgElem= new CHexa; pProgElem->setMGLevel(ilevel+1);
                vProgElem.push_back(pProgElem);
            };
            dividTetra(pElem,vProgElem, elementID, pProgMesh);
            break;
        case(ElementType::Prism):case(ElementType::Prism2):
            vProgElem.reserve(6);
            for(i=0; i< 6; i++){
                pProgElem= new CHexa; pProgElem->setMGLevel(ilevel+1);
                vProgElem.push_back(pProgElem);
            };
            dividPrism(pElem,vProgElem, elementID,pProgMesh);
            break;
        case(ElementType::Pyramid):case(ElementType::Pyramid2):
            vProgElem.reserve(8);
            for(i=0; i< 4; i++){
                pProgElem= new CHexa; pProgElem->setMGLevel(ilevel+1);
                vProgElem.push_back(pProgElem);
            };
            for(i=0; i< 4; i++){
                pProgElem= new CPyramid; pProgElem->setMGLevel(ilevel+1);
                vProgElem.push_back(pProgElem);
            };
            dividPyramid(pElem,vProgElem, elementID,pProgMesh);
            break;
        case(ElementType::Quad):case(ElementType::Quad2):
            vProgElem.reserve(4);
            for(i=0; i< 4; i++){
                pProgElem= new CQuad; pProgElem->setMGLevel(ilevel+1);
                vProgElem.push_back(pProgElem);
            };
            dividQuad(pElem,vProgElem, elementID,pProgMesh);
            break;
        case(ElementType::Triangle):case(ElementType::Triangle2):
            vProgElem.reserve(3);
            for(i=0; i< 3; i++){
                pProgElem= new CQuad; pProgElem->setMGLevel(ilevel+1);
                vProgElem.push_back(pProgElem);
            };
            dividTriangle(pElem,vProgElem, elementID,pProgMesh);
            break;
        case(ElementType::Beam):case(ElementType::Beam2):
            vProgElem.reserve(2);
            for(i=0; i< 2; i++){
                pProgElem= new CBeam; pProgElem->setMGLevel(ilevel+1);
                vProgElem.push_back(pProgElem);
            };
            dividBeam(pElem,vProgElem, elementID,pProgMesh);
            break;
    }
}
void CMeshFactory::dividHexa(CElement* pElem, vector<CElement*>& vProgElem, uint& elementID, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;
    vector<CNode*> vEdgeNode;
    vector<CNode*> vFaceNode;
    CNode          *pVolNode;
    uint i;
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
        vProgElem[i]->setParentID(pElem->getID());
        vProgElem[i]->setID(elementID);          
        ++elementID;
        pProgMesh->setElement(vProgElem[i]);
    };
    uint iface;
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
void CMeshFactory::setProgHexaMPCMaster(vector<CElement*>& vProgElem, const uint& iface, const uint& i, const uint& j, const uint& k, const uint& l)
{
    vProgElem[i]->markingMPCMaster(); vProgElem[i]->markingMPCFace(iface);
    vProgElem[j]->markingMPCMaster(); vProgElem[j]->markingMPCFace(iface);
    vProgElem[k]->markingMPCMaster(); vProgElem[k]->markingMPCFace(iface);
    vProgElem[l]->markingMPCMaster(); vProgElem[l]->markingMPCFace(iface);
}
void CMeshFactory::setProgHexaMPCSlave(vector<CElement*>& vProgElem, const uint& iface, const uint& i, const uint& j, const uint& k, const uint& l)
{
    vProgElem[i]->markingMPCSlave(); vProgElem[i]->markingMPCFace(iface);
    vProgElem[j]->markingMPCSlave(); vProgElem[j]->markingMPCFace(iface);
    vProgElem[k]->markingMPCSlave(); vProgElem[k]->markingMPCFace(iface);
    vProgElem[l]->markingMPCSlave(); vProgElem[l]->markingMPCFace(iface);
}
void CMeshFactory::setProgHexaCommMesh2Entity(vector<CElement*>& vProgElem, const uint& iface, const uint& i, const uint& j, const uint& k, const uint& l)
{
    vProgElem[i]->markingCommMesh2(); vProgElem[i]->markingCommEntity(iface);
    vProgElem[j]->markingCommMesh2(); vProgElem[j]->markingCommEntity(iface);
    vProgElem[k]->markingCommMesh2(); vProgElem[k]->markingCommEntity(iface);
    vProgElem[l]->markingCommMesh2(); vProgElem[l]->markingCommEntity(iface);
}
void CMeshFactory::dividTetra(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;
    vector<CNode*> vEdgeNode;
    vector<CNode*> vFaceNode;
    CNode          *pVolNode;
    uint numOfVert,numOfEdge,numOfFace;
    numOfVert= NumberOfVertex::Tetra(); numOfFace= NumberOfFace::Tetra(); numOfEdge= NumberOfEdge::Tetra();
    uint i;
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
        vProgElem[i]->setParentID(pElem->getID());
        vProgElem[i]->setID(indexCount);
        ++indexCount;
        pProgMesh->setElement(vProgElem[i]);
    };
    uint iface;
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
void CMeshFactory::dividPrism(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;
    vector<CNode*> vEdgeNode;
    vector<CNode*> vFaceNode;
    CNode          *pVolNode;
    uint numOfVert,numOfEdge,numOfFace;
    numOfVert= NumberOfVertex::Prism(); numOfFace= NumberOfFace::Prism(); numOfEdge= NumberOfEdge::Prism();
    uint i;
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
        vProgElem[i]->setParentID(pElem->getID());
        vProgElem[i]->setID(indexCount);
        ++indexCount;
        pProgMesh->setElement(vProgElem[i]);
    };
    uint iface;
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
void CMeshFactory::dividPyramid(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;
    vector<CNode*> vEdgeNode;
    vector<CNode*> vFaceNode;
    CNode          *pVolNode;
    uint numOfVert,numOfEdge,numOfFace;
    numOfVert= NumberOfVertex::Pyramid(); numOfFace= NumberOfFace::Pyramid(); numOfEdge= NumberOfEdge::Pyramid();
    uint i;
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
        vProgElem[i]->setParentID(pElem->getID());
        vProgElem[i]->setID(indexCount);
        ++indexCount;
        pProgMesh->setElement(vProgElem[i]);
    };
}
void CMeshFactory::dividQuad(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;
    vector<CNode*> vEdgeNode;
    vector<CNode*> vFaceNode;
    uint numOfVert,numOfEdge,numOfFace;
    numOfVert= NumberOfVertex::Quad(); numOfFace= NumberOfFace::Quad(); numOfEdge= NumberOfEdge::Quad();
    uint i;
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
        vProgElem[i]->setParentID(pElem->getID());
        vProgElem[i]->setID(indexCount);
        ++indexCount;
        pProgMesh->setElement(vProgElem[i]);
    };
    uint iprog;
    if(pElem->isMPCMaster()){
        for(iprog=0; iprog< 4; iprog++){ vProgElem[iprog]->markingMPCMaster(); vProgElem[iprog]->markingMPCFace(0);}
    }
    if(pElem->isMPCSlave()){
        for(iprog=0; iprog< 4; iprog++){ vProgElem[iprog]->markingMPCSlave(); vProgElem[iprog]->markingMPCFace(0);}
    }
    uint iedge;
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
void CMeshFactory::dividTriangle(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;
    vector<CNode*> vEdgeNode;
    vector<CNode*> vFaceNode;
    uint numOfVert,numOfEdge,numOfFace;
    numOfVert= NumberOfVertex::Triangle(); numOfFace= NumberOfFace::Triangle(); numOfEdge= NumberOfEdge::Triangle();
    uint i;
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
        vProgElem[i]->setParentID(pElem->getID());
        vProgElem[i]->setID(indexCount);
        ++indexCount;
        pProgMesh->setElement(vProgElem[i]);
    };
    uint iprog;
    if(pElem->isMPCMaster()){
        for(iprog=0; iprog< 3; iprog++){ vProgElem[iprog]->markingMPCMaster(); vProgElem[iprog]->markingMPCFace(0);}
    }
    if(pElem->isMPCSlave()){
        for(iprog=0; iprog< 3; iprog++){ vProgElem[iprog]->markingMPCSlave(); vProgElem[iprog]->markingMPCFace(0);}
    }
    uint iedge;
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
void CMeshFactory::dividBeam(CElement* pElem, vector<CElement*>& vProgElem, uint& indexCount, CMesh* pProgMesh)
{
    vector<CNode*> vVertNode;
    vector<CNode*> vEdgeNode;
    uint numOfVert,numOfEdge;
    numOfVert= NumberOfVertex::Beam(); numOfEdge= NumberOfEdge::Beam();
    uint i;
    vVertNode.resize(numOfVert); for(i=0; i< numOfVert; i++){ vVertNode[i] = pElem->getNode(i);}
    vEdgeNode.resize(numOfEdge); for(i=0; i< numOfEdge; i++){ vEdgeNode[i] = pElem->getEdgeInterNode(i);}
    vProgElem[0]->setNode(vVertNode[0],0); vProgElem[0]->setNode(vEdgeNode[0],1);
    pElem->setProgElem(vProgElem[0], 0);
    vProgElem[1]->setNode(vEdgeNode[0],0); vProgElem[1]->setNode(vVertNode[1],1);
    pElem->setProgElem(vProgElem[1], 1);
    for(i=0; i< 2; i++){
        vProgElem[i]->setParentID(pElem->getID());
        vProgElem[i]->setID(indexCount);
        ++indexCount;
        pProgMesh->setElement(vProgElem[i]);
    };
    uint iprog;
    if(pElem->isMPCMaster()){
        for(iprog=0; iprog< 2; iprog++) vProgElem[iprog]->markingMPCMaster();
    }
    if(pElem->isMPCSlave()){
        for(iprog=0; iprog< 2; iprog++) vProgElem[iprog]->markingMPCSlave();
    }
    uint ivert;
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
void CMeshFactory::setupBucketMesh(const uint& mgLevel, const uint& maxID, const uint& minID)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTAssyModel->intializeBucket(maxID, minID);
    mpTAssyModel->setMaxMeshID(maxID);
    mpTAssyModel->setMinMeshID(minID);
}
void CMeshFactory::GeneAssyModel(const uint& num_of_mgLevel)
{
    mpGMGModel->initAssyModel(num_of_mgLevel);
    mpGMGModel->reserveAssyModel(num_of_mgLevel);
    for(uint i=0; i<num_of_mgLevel; i++){
        mpTAssyModel = new CAssyModel();
        mpTAssyModel->setMGLevel(i);
        mpGMGModel->addModel(mpTAssyModel,i);
    };
}
void CMeshFactory::reserveMesh(const uint& mgLevel, const uint& num_of_mesh)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTAssyModel->resizeMesh(num_of_mesh);
}
void CMeshFactory::GeneMesh(const uint& mgLevel, const uint& mesh_id, const uint& index)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = new CMesh();
    mpTMesh->setMeshID(mesh_id);
    mpTMesh->setMGLevel(mgLevel);
    mpTAssyModel->setBucket(mesh_id, index);
    mpTAssyModel->setMesh(mpTMesh,index);
}
void CMeshFactory::GeneNode(const uint& mgLevel, const uint& mesh_id, const uint& id, const vdouble& coord,
                            const uint& nodeType, const uint& numOfScaParam, const uint& numOfVecParam)
{
    CNode *pNode;
    switch(nodeType){
        case(NodeType::Scalar):
            pNode = new CScalarNode();
            pNode->resizeScalar(numOfScaParam);
            break;
        case(NodeType::Vector):
            pNode = new CVectorNode();
            pNode->resizeVector(numOfVecParam);
            break;
        case(NodeType::ScalarVector):
            pNode = new CScalarVectorNode();
            pNode->resizeScalar(numOfScaParam);
            pNode->resizeVector(numOfVecParam);
            break;
        default:
            break;
    }
    pNode->setID(id);
    pNode->setCoord(coord);
    pNode->setMGLevel(mgLevel);
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    if(!mpTAssyModel) mpLogger->Info(Utility::LoggerMode::MWDebug, "AssyModel => NULL");
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->setNode(pNode);
}
void CMeshFactory::setupNode(const uint& mgLevel, const uint& mesh_id)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh(mesh_id);
    mpTMesh->setupNumOfNode();
}
void CMeshFactory::reserveNode(const uint& mgLevel, const uint& mesh_id, const uint& num_of_node)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->reserveNode(num_of_node);
}
void CMeshFactory::GeneElement(const uint& mgLevel, const uint& mesh_id, const uint& id, const uint& type_num, const vint& node_id)
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
        case(ElementType::Pyramid):
            pElement = new CPyramid;
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
        default:
            mpLogger->Info(Utility::LoggerMode::Error,"Error::GeneElement at Factory");
            break;
    }
    pElement->setID(id);
    CNode *pNode;
    uint i;
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);    
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    for(i=0; i < node_id.size(); i++){
        pNode = mpTMesh->getNode(node_id[i]);
        pElement->setNode(pNode, i);
    };
    pElement->setupFaceCnvNodes();
    mpTMesh->setElement(pElement);
}
void CMeshFactory::setupElement(const uint& mgLevel, const uint& mesh_id)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh(mesh_id);
    mpTMesh->setupNumOfElement();
}
void CMeshFactory::reserveElement(const uint& mgLevel, const uint& mesh_id, const uint& num_of_element)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->reserveElement(num_of_element);
}
void CMeshFactory::reserveAggregate(const uint& mgLevel, const uint& mesh_id, const uint& num_of_node)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->reserveAggregate(num_of_node);
}
void CMeshFactory::GeneAggregate(const uint& mgLevel, const uint& mesh_id, const uint& num_of_node)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    uint i;
    for(i=0; i< num_of_node; i++){
        CAggregateElement *pAggElem = new CAggregateElement;
        CAggregateNode    *pAggNode = new CAggregateNode;
        mpTMesh->setAggElement(pAggElem);
        mpTMesh->setAggNode(pAggNode);
    };
}
void CMeshFactory::reserveBoundaryNode(const uint& mgLevel, const uint& mesh_id, const uint& num_of_bnd)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->reserveBoundaryNode(num_of_bnd);
}
void CMeshFactory::reserveBoundaryFace(const uint& mgLevel, const uint& mesh_id, const uint& num_of_bnd)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->reserveBoundaryFace(num_of_bnd);
}
void CMeshFactory::reserveBoundaryVolume(const uint& mgLevel, const uint& mesh_id, const uint& num_of_bnd)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->reserveBoundaryVolume(num_of_bnd);
}
void CMeshFactory::GeneBoundaryNode(const uint& mgLevel, const uint& mesh_id, const uint& id,
                                    const uint& dof, const uint& bndType, const vdouble& vVal)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    CBoundaryNode *pBNode = new CBoundaryNode();
    uint i;
    switch(bndType){
        case(BoundaryTypeNode::Accel):
            pBNode->initialize(BoundaryTypeNode::Accel, dof);
            for(i=0; i< dof; i++){ pBNode->setValue(vVal[i],i);}
            break;
        case(BoundaryTypeNode::Velo):
            pBNode->initialize(BoundaryTypeNode::Velo, dof);
            for(i=0; i< dof; i++){ pBNode->setValue(vVal[i],i);}
            break;
        case(BoundaryTypeNode::Load):
            pBNode->initialize(BoundaryTypeNode::Load, dof);
            for(i=0; i< dof; i++){ pBNode->setValue(vVal[i],i);}
            break;
        case(BoundaryTypeNode::Disp):
            pBNode->initialize(BoundaryTypeNode::Disp, dof);
            for(i=0; i< dof; i++){ pBNode->setValue(vVal[i],i);}
            break;
        case(BoundaryTypeNode::Temp):
            pBNode->initialize(BoundaryTypeNode::Temp, dof);
            pBNode->setValue(vVal[0],0);
            break;
        case(BoundaryTypeNode::Thermal_Flux):
            pBNode->initialize(BoundaryTypeNode::Thermal_Flux, dof);
            pBNode->setValue(vVal[0], 0);
            break;
        default:
            break;
    }
    pBNode->setID(id);
    mpTMesh->setBoundaryNode(pBNode);
}
void CMeshFactory::GeneBoundaryFace(const uint& mgLevel, const uint& mesh_id, const uint& elem_id, const uint& face_id,
                                    const uint& dof, const uint& bndType, const vdouble& vVal)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    CBoundaryFace *pBFace = new CBoundaryFace();
    uint i;
    switch(bndType){
        case(BoundaryTypeFace::Pressure):
            pBFace->initialize(BoundaryTypeFace::Pressure,dof);
            for(i=0; i< dof; i++){ pBFace->setValue(vVal[i], i);}
            break;
        case(BoundaryTypeFace::Temp):
            pBFace->initialize(BoundaryTypeFace::Temp, dof);
            for(i=0; i< dof; i++){ pBFace->setValue(vVal[i], i);}
            break;
        case(BoundaryTypeFace::Thermal_Flux):
            pBFace->initialize(BoundaryTypeFace::Thermal_Flux, dof);
            for(i=0; i< dof; i++){ pBFace->setValue(vVal[i], i);}
            break;
        case(BoundaryTypeFace::TractionVector):
            pBFace->initialize(BoundaryTypeFace::TractionVector, dof);
            for(i=0; i< dof; i++){ pBFace->setValue(vVal[i], i);}
            break;
        default:
            break;
    }
    pBFace->setElementID(elem_id);
    pBFace->setFaceID(face_id);
    mpTMesh->setBoundaryFace(pBFace);
}
void CMeshFactory::GeneBoundaryVolume(const uint& mgLevel, const uint& mesh_id, const uint& id,
                                     const uint& dof, const uint& bndType, const vdouble& vVal)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    CBoundaryVolume* pBVolume = new CBoundaryVolume();
    uint i;
    switch(bndType){
        case(BoundaryTypeVolume::Accel):
            pBVolume->initialize(BoundaryTypeVolume::Accel, dof);
            for(i=0; i< dof; i++){ pBVolume->setValue(vVal[i], i);}
            break;
        case(BoundaryTypeVolume::Centrifugal_Force):
            pBVolume->initialize(BoundaryTypeVolume::Centrifugal_Force, dof);
            for(i=0; i< dof; i++){ pBVolume->setValue(vVal[i], i);}
            break;
        case(BoundaryTypeVolume::Gravity):
            pBVolume->initialize(BoundaryTypeVolume::Gravity, dof);
            for(i=0; i< dof; i++){ pBVolume->setValue(vVal[i], i);}
            break;
        case(BoundaryTypeVolume::Heat):
            pBVolume->initialize(BoundaryTypeVolume::Heat, dof);
            for(i=0; i< dof; i++){ pBVolume->setValue(vVal[i], i);}
            break;
        default:
            break;
    }
    pBVolume->setID(id);
    mpTMesh->setBoundaryVolume(pBVolume);
}
void CMeshFactory::setupBucketNode(const uint& mgLevel, const uint& mesh_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->setupBucketNode();
}
void CMeshFactory::initBucketNode(const uint& mgLevel, const uint& mesh_id, const uint& maxID, const uint& minID)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->initBucketNode(maxID, minID);
}
void CMeshFactory::setIDBucketNode(const uint& mgLevel, const uint& mesh_id, const uint& id, const uint& index)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->setupBucketNodeIndex(id, index);
}
void CMeshFactory::setupBucketElement(const uint& mgLevel, const uint& mesh_id)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->setupBucketElement();
}
void CMeshFactory::initBucketElement(const uint& mgLevel, const uint& mesh_id, const uint& maxID, const uint& minID)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    ;
    mpTMesh->initBucketElement(maxID, minID);
}
void CMeshFactory::setIDBucketElement(const uint& mgLevel, const uint& mesh_id, const uint& id, const uint& index)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->setupBucketElementIndex(id, index);
}
void CMeshFactory::reserveMaterial(const uint& res_size)
{
     mpGMGModel->reserveMaterial(res_size);
}
void CMeshFactory::GeneMaterial(const uint& mesh_id, const uint& material_id, string& name, vuint& vType, vdouble& vValue)
{
    CMaterial *pMaterial = new CMaterial;
    pMaterial->setID(material_id);
    pMaterial->setName(name);
    pMaterial->setMeshID(mesh_id);
    for(uint i=0; i< vType.size(); i++) pMaterial->setValue(vType[i],vValue[i]);
    mpGMGModel->setMaterial(pMaterial);
}
void CMeshFactory::reserveCommMesh(const uint& mgLevel, const uint& mesh_id, const uint& res_size)
{
    mpTAssyModel = mpGMGModel->getAssyModel(mgLevel);
    mpTMesh = mpTAssyModel->getMesh(mesh_id);
    mpTMesh->reserveCommMesh(res_size);
}
void CMeshFactory::GeneCommMesh(const uint& mgLevel, const uint& mesh_id, const uint& comID, const uint& myRank, const uint& nTransmitRank)
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
void CMeshFactory::reserveCommNode(const uint& mgLevel, const uint& mesh_id, const uint& commesh_id, const uint& res_size)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh(mesh_id);
    mpTCommMesh= mpTMesh->getCommMesh(commesh_id);
    mpTCommMesh->reserveNode(res_size);   
    mpTCommMesh->resizeNodeRank(res_size);
}
void CMeshFactory::GeneCommNode(const uint& mgLevel, const uint& commNodeID,
                                     const uint& mesh_id, const uint& commesh_id, const uint& nodeID, const uint& rank)
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
void CMeshFactory::reserveCommElement(const uint& mgLevel, const uint& mesh_id, const uint& commesh_id, const uint& res_size)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh(mesh_id);
    mpTCommMesh= mpTMesh->getCommMesh(commesh_id);
    mpTCommMesh->reserveCommElementAll(res_size);
}
void CMeshFactory::GeneCommElement(const uint& mgLevel, const uint& mesh_id, const uint& commesh_id, 
                                   const uint& nType, const uint& elemID, vuint& vCommNodeID)
{
    mpTAssyModel= mpGMGModel->getAssyModel(mgLevel);
    mpTMesh= mpTAssyModel->getMesh(mesh_id);
    mpTCommMesh= mpTMesh->getCommMesh(commesh_id);
    uint numOfVert= vCommNodeID.size();
    vuint vNodeRank;  vNodeRank.reserve(numOfVert);
    uint ivert, rank;
    uint commNodeID;
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
        case(ElementType::Pyramid):
            pCommElem = new CCommPyramid;
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
    uint ivert;
    switch(pCommElem->getShapeType()){
        case(ElementType::Hexa):
            vProgCommElem.reserve(8);
            for(ivert=0; ivert< 8; ivert++){
                pProgCommElem= new CCommHexa;
                vProgCommElem.push_back(pProgCommElem);
            };
            dividCommElem(pCommElem, vProgCommElem);
			mpLogger->Info(Utility::LoggerMode::MWDebug,"CommElem(Hexa)の分割");
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
        case(ElementType::Pyramid):
            vProgCommElem.reserve(8);
            for(ivert=0; ivert< 8; ivert++){
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
    uint nRank;
    uint ivert, progvert;
    uint iedge, iface;
    uint invalid= pProgTree->getInvalidNum();
    uint numOfVert, numOfEdge, numOfFace;
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
void CMeshFactory::GeneContactMesh(const uint& contactID, const uint& myRank, const uint& transRank)
{
    uint ilevel;
    for(ilevel=0; ilevel< mMGLevel+1; ilevel++){
        mpTAssyModel= mpGMGModel->getAssyModel(ilevel);
        CContactMesh *pContactMesh= new CContactMesh;
        pContactMesh->setID(contactID);
        pContactMesh->setLevel(ilevel);
        pContactMesh->setRank(myRank);
        pContactMesh->setTransmitRank(transRank);
        mpTAssyModel->addContactMesh(pContactMesh, contactID);
    };
}
void CMeshFactory::GeneContactNode(const uint& mgLevel, const uint& contactID, const uint& conNodeID, const vdouble& vCoord,
        const string& s_param_type, const uint& numOfVector, const uint& numOfScalar,
        bool bmesh, const uint& meshID, const uint& nodeID,
        const uint& rank, const uint& maslave)
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
void CMeshFactory::GeneMasterFace(const uint& contactID, const uint& shapeType, const uint& masterFaceID,
        bool bmesh, const uint& meshID, const uint& elemID, const uint& elemFaceID,
        const vuint& vConNodeID, const uint& face_rank)
{
    uint mgLevel(0);
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
    uint icnode, numOfConNode(vConNodeID.size());
    for(icnode=0; icnode< numOfConNode; icnode++){
        CContactNode *pConNode= pConMesh->getContactNode_ID(vConNodeID[icnode]);
        pMFace->addNode(pConNode);
    };
    pConMesh->addMasterFace(pMFace);
}
void CMeshFactory::GeneSlaveFace(const uint& contactID, const uint& shapeType, const uint& slaveFaceID,
        bool bmesh, const uint& meshID, const uint& elemID, const uint& elemFaceID,
        const vuint& vConNodeID, const uint& face_rank)
{
    uint mgLevel(0);
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
    uint icnode, numOfConNode(vConNodeID.size());
    for(icnode=0; icnode< numOfConNode; icnode++){
        CContactNode *pConNode= pConMesh->getContactNode_ID(vConNodeID[icnode]);
        pSFace->addNode(pConNode);
    };
    pConMesh->addSlaveFace(pSFace);
}
void CMeshFactory::refineContactMesh()
{
    CAssyModel *pAssy,*pProgAssy;
    CContactMesh *pConMesh,*pProgConMesh;
    CSkinFace *pSkinFace;
    vector<CSkinFace*> vProgFace;
    uint maslave;
    uint meshID,elemID;
    CMesh *pMesh;
    CElement *pElem;
    uint faceID;
    uint maxLayer;
    uint ilevel;
    for(ilevel=0; ilevel< mMGLevel; ilevel++){
        pAssy =  mpGMGModel->getAssyModel(ilevel);    
        pProgAssy= mpGMGModel->getAssyModel(ilevel+1);
        uint numOfCont= pAssy->getNumOfContactMesh();
        uint icont;
        for(icont=0; icont< numOfCont; icont++){
            pConMesh= pAssy->getContactMesh(icont);
            pProgConMesh= pProgAssy->getContactMesh(icont);
            pConMesh->setupCoarseConNode(pProgConMesh);
            pConMesh->setupAggSkinFace();            
            pConMesh->setupEdgeConNode(pProgConMesh);
            pConMesh->setupFaceConNode(pProgConMesh);
            uint numOfSkinFace;
            uint iface;
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
    for(ilevel=0; ilevel < mMGLevel+1; ilevel++){
        pAssy= mpGMGModel->getAssyModel(ilevel);
        uint numOfCont= pAssy->getNumOfContactMesh();
        uint icont;
        for(icont=0; icont < numOfCont; icont++){
            pConMesh= pAssy->getContactMesh(icont);
            maxLayer= ilevel+1;
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
    uint countID(0);
    uint ilevel;
    for(ilevel=0; ilevel< mMGLevel; ilevel++){
        pAssy= mpGMGModel->getAssyModel(ilevel);
        pProgAssy= mpGMGModel->getAssyModel(ilevel+1);
        uint imesh, numOfMesh;
        numOfMesh= pAssy->getNumOfMesh();
        for(imesh=0; imesh< numOfMesh; imesh++){
            pMesh= pAssy->getMesh(imesh);
            pProgMesh= pProgAssy->getMesh(imesh);
            uint icomm, numOfComm;
            numOfComm= pMesh->getCommMesh2Size();
            for(icomm=0; icomm< numOfComm; icomm++){
                pCommMesh2= pMesh->getCommMesh2IX(icomm);
                pProgCommMesh2 = new CCommMesh2;
                pProgCommMesh2->setLevel(ilevel+1);
                pProgCommMesh2->setID(pCommMesh2->getID());
                pProgCommMesh2->setRank(pCommMesh2->getRank());
                pProgCommMesh2->setTransmitRank(pCommMesh2->getTrasmitRank());
                pProgMesh->setCommMesh2(pProgCommMesh2);
                pCommMesh2->setupVertCommNode(pProgCommMesh2);
                pCommMesh2->setupAggFace();
                pCommMesh2->setupEdgeCommNode(pProgCommMesh2);
                pCommMesh2->setupFaceCommNode(pProgCommMesh2);
                uint iface, numOfFace;
                numOfFace= pCommMesh2->getCommFaceSize();
                for(iface=0; iface< numOfFace; iface++){
                    pCommFace= pCommMesh2->getCommFaceIX(iface);
                    uint elemID = pCommFace->getElementID();
                    CElement *pElem= pMesh->getElement(elemID);;
                    mvCommFace= pCommFace->refine(pElem);
                    uint ipface,numOfProgFace;
                    numOfProgFace= mvCommFace.size();
                    for(ipface=0; ipface< numOfProgFace; ipface++){
                        pProgCommFace= mvCommFace[ipface];
                        pProgCommFace->setID(countID);
                        pProgCommMesh2->addCommFace(pProgCommFace);
                        countID++;
                    };
                };
                uint elemID,entity_num;
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
                    uint iedge, numOfEdge;
                    numOfEdge= pCommFace->getNumOfEdge();
                    for(iedge=0; iedge< numOfEdge; iedge++){
                        pairCommNode= pCommFace->getEdgePairCommNode(iedge);
                        pNodeFir= pairCommNode.first->getNode();
                        pNodeSec= pairCommNode.second->getNode();
                        uint edgeIndex;
                        edgeIndex= pElem->getEdgeIndex(pNodeFir,pNodeSec);
                        pEdgeNode= pElem->getEdgeInterNode(edgeIndex);
                        pEdgeCommNode= pCommFace->getEdgeCommNode(iedge);
                        pEdgeCommNode->setNode(pEdgeNode);
                    };
                };
            };
        };
    };
}
void CMeshFactory::GeneCommMesh2(const uint& mgLevel, const uint& mesh_id, const uint& comID, 
        const uint& numOfFace, const uint& numOfCommNode,
        const uint& myRank, const uint& nTransmitRank)
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
void CMeshFactory::GeneCommFace(const uint& mgLevel, const uint& commeshID, const uint& face_id,
            const uint& mesh_id,const uint elem_id, const uint& elem_ent_num, const vuint& vCommNodeID)
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
    uint numOfVert= vCommNodeID.size();
    uint numOfEdge;
    switch(numOfVert)
    {
        case(4):
            numOfEdge= 4;
            break;
        case(3):
            numOfEdge= 3;
            break;
        case(2):
            numOfEdge= 1;
            break;
        default:
            break;
    }
    pCommFace->initialize(numOfVert,numOfEdge);
    uint ivert, id;
    numOfVert= vCommNodeID.size();
    CCommNode *pCommNode;
    for(ivert=0; ivert< numOfVert; ivert++){
        id= vCommNodeID[ivert];
        pCommNode= mpTCommMesh2->getCommVertNode(id);
        pCommFace->setVertCommNode(ivert, pCommNode);
    };
    mpTCommMesh2->addCommFace(pCommFace);
}
void CMeshFactory::GeneCommNodeCM2(const uint& mgLevel, const uint& mesh_id, const uint& node_id,const uint& commeshID,
        const uint& comm_node_id, const vdouble& vCoord)
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
