/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterBoundaryNode.cpp
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
#include "BoundaryNode.h"
#include "BoundaryNodeMesh.h"
#include "FileWriterBoundaryNode.h"
#include "./BoundaryType.h"
using namespace FileIO;
CFileWriterBoundaryNode::CFileWriterBoundaryNode()
{
    ;
}
CFileWriterBoundaryNode::~CFileWriterBoundaryNode()
{
    ;
}
void CFileWriterBoundaryNode::WriteDebug(ofstream& ofs, const uiint& mgLevel)
{
    pmw::CAssyModel *pAssyModel;
    pmw::CMesh *pMesh;
    pmw::CBoundaryNodeMesh *pBndNodeMesh;
    pmw::CBoundarySBNode *pBndNode;
    pmw::CNode *pNode;
    pAssyModel= mpGMGModel->getAssyModel(mgLevel);
    uiint i,ii,iii;
    uiint numOfPart= pAssyModel->getNumOfMesh();
    for(i=0; i< numOfPart; i++) {
        pMesh= pAssyModel->getMesh(i);
        uiint numOfBndNodeMesh= pMesh->getNumOfBoundaryNodeMesh();
        for(ii=0; ii< numOfBndNodeMesh; ii++) {
            pBndNodeMesh= pMesh->getBndNodeMeshIX(ii);
            ofs << "BoundaryID= " << pBndNodeMesh->getID() << ", BndType= ";
            switch(pBndNodeMesh->getBndType()) {
            case(pmw::BoundaryType::Dirichlet):
                ofs << "Dirichlet";
                break;
            case(pmw::BoundaryType::Neumann):
                ofs << "Neumann";
                break;
            default:
                break;
            }
            ofs << endl;
            uiint numOfBndNode= pBndNodeMesh->getNumOfBNode();
            for(iii=0; iii < numOfBndNode; iii++) {
                pBndNode= pBndNodeMesh->getBNodeIX(iii);
                pNode= pBndNode->getNode();
                ofs << "BNodeID= " << pBndNode->getID() << ", NodeID= " << pNode->getID();
                uiint idof, numOfDOF= pBndNode->getNumOfDOF();
                for(idof=0; idof < numOfDOF; idof++) {
                    uiint dof= pBndNode->getDOF(idof);
                    ofs << ", BndValue[" << dof << "]= " << pBndNode->getValue(dof);
                };
                ofs << endl;
            };
        };
    };
}
