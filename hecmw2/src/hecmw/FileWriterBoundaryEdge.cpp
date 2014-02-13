/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterBoundaryEdge.cpp
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
#include "BoundaryEdge.h"
#include "BoundaryEdgeMesh.h"
#include "FileWriterBoundaryEdge.h"
using namespace FileIO;
CFileWriterBoundaryEdge::CFileWriterBoundaryEdge()
{
    ;
}
CFileWriterBoundaryEdge::~CFileWriterBoundaryEdge()
{
    ;
}
void CFileWriterBoundaryEdge::WriteDebug(ofstream& ofs, const uiint& mgLevel)
{
    pmw::CAssyModel *pAssyModel;
    pmw::CMesh *pMesh;
    pmw::CBoundaryEdgeMesh *pBndEdgeMesh;
    pmw::CBoundaryEdge *pBndEdge;
    pmw::CBoundaryNode *pBndNode;
    pmw::CNode *pNode;
    pAssyModel= mpGMGModel->getAssyModel(mgLevel);
    uiint i,ii,iii;
    uiint numOfPart= pAssyModel->getNumOfMesh();
    for(i=0; i < numOfPart; i++) {
        pMesh= pAssyModel->getMesh(i);
        uiint numOfBndMesh= pMesh->getNumOfBoundaryEdgeMesh();
        for(ii=0; ii < numOfBndMesh; ii++) {
            pBndEdgeMesh= pMesh->getBndEdgeMeshIX(ii);
            ofs << "BoundaryID= " << pBndEdgeMesh->getID() << ", BndType= ";
            switch(pBndEdgeMesh->getBndType()) {
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
            uiint numOfBndNode= pBndEdgeMesh->getNumOfBNode();
            vdouble sumValue;
            sumValue.resize(pBndEdgeMesh->getNumOfDOF());
            for(uiint idof=0; idof < pBndEdgeMesh->getNumOfDOF(); idof++)  sumValue[idof]=0.0;
            for(iii=0; iii < numOfBndNode; iii++) {
                pBndNode= pBndEdgeMesh->getBNodeIX(iii);
                pNode= pBndNode->getNode();
                ofs << "BNodeID= " << pBndNode->getID() << ", NodeID= " << pNode->getID();
                uiint idof, numOfDOF= pBndEdgeMesh->getNumOfDOF();
                for(idof=0; idof < numOfDOF; idof++) {
                    uiint dof = pBndEdgeMesh->getDOF(idof);
                    ofs << ", dof= " << dof << ", Value[" << idof << "]= " << pBndNode->getValue(dof, mgLevel);
                    sumValue[idof] += pBndNode->getValue(dof, mgLevel);
                };
                ofs << endl;
            };
            if(pBndEdgeMesh->getBndType()==pmw::BoundaryType::Neumann) {
                ofs << "BNode Value SumValue :";
                for(uiint idof=0; idof < pBndEdgeMesh->getNumOfDOF(); idof++) {
                    ofs << "  sumValue[" << idof << "]=" << sumValue[idof];
                };
                ofs << endl;
            }
            if(pBndEdgeMesh->getBndType()==pmw::BoundaryType::Neumann) {
                uiint numOfBndEdge= pBndEdgeMesh->getNumOfEdge();
                for(iii=0; iii < numOfBndEdge; iii++) {
                    pBndEdge= pBndEdgeMesh->getBEdgeIX(iii);
                    ofs << "BEdgeID= " << pBndEdge->getID();
                    uiint idof, numOfDOF= pBndEdgeMesh->getNumOfDOF();
                    for(idof=0; idof < numOfDOF; idof++) {
                        uiint dof = pBndEdgeMesh->getDOF(idof);
                        ofs << ", Value[" << idof << "]= " << pBndEdge->getBndValue(dof);
                    };
                    ofs << endl;
                };
            }
        };
    };
}
