/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/BoundaryNodeMesh.cpp
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
#include "BoundaryNodeMesh.h"
using namespace pmw;
CBoundaryNodeMesh::CBoundaryNodeMesh()
{
    ;
}
CBoundaryNodeMesh::~CBoundaryNodeMesh()
{
    ;
}
void CBoundaryNodeMesh::setBndType(const uiint& boundType)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    switch(boundType){
        case(BoundaryType::Dirichlet):
            mnBndType= boundType;
            break;
        case(BoundaryType::Neumann):
            mnBndType= boundType;
            break;
        default:
            pLogger->Info(Utility::LoggerMode::Error, "BoundaryType Error, CBoundaryNodeMesh::setType");
            break;
    }
}
void CBoundaryNodeMesh::resizeBNode(const uiint& res_size)
{
    mvBNode.resize(res_size);
}
void CBoundaryNodeMesh::setBNode(const uiint& index, CBoundarySBNode* pBNode)
{
    uiint id= pBNode->getID();
    mmID2Index[id]= index;
    CNode *pNode= pBNode->getNode();
    uiint node_id= pNode->getID();
    mmNodeID2BNodeID[node_id]= pBNode->getID();
}
void CBoundaryNodeMesh::addBNode(CBoundarySBNode* pBNode)
{
    mvBNode.push_back(pBNode);
    uiint index= mvBNode.size()-1;
    uiint id= pBNode->getID();
    mmID2Index[id]= index;
    CNode *pNode= pBNode->getNode();
    uiint node_id= pNode->getID();
    mmNodeID2BNodeID[node_id]= pBNode->getID();
}
