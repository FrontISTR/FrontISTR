/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterBoundaryFace.cpp
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
#include "BoundaryMesh.h"
#include "BoundaryNode.h"
#include "BoundaryFace.h"
#include "BoundaryFaceMesh.h"
#include "FileWriterBoundaryFace.h"
using namespace FileIO;
CFileWriterBoundaryFace::CFileWriterBoundaryFace()
{
    ;
}
CFileWriterBoundaryFace::~CFileWriterBoundaryFace()
{
    ;
}
void CFileWriterBoundaryFace::WriteDebug(ofstream& ofs, const uiint& mgLevel)
{
    pmw::CAssyModel *pAssyModel;
    pmw::CMesh *pMesh;
    pmw::CBoundaryFaceMesh *pBndFaceMesh;
    pmw::CBoundaryFace *pBndFace;
    pmw::CBoundaryNode *pBndNode;
    pmw::CNode *pNode;
    pAssyModel= mpGMGModel->getAssyModel(mgLevel);
    uiint i,ii,iii;
    uiint numOfPart= pAssyModel->getNumOfMesh();
    for(i=0; i < numOfPart; i++) {
        pMesh= pAssyModel->getMesh(i);
        uiint numOfBndMesh= pMesh->getNumOfBoundaryFaceMesh();
        for(ii=0; ii < numOfBndMesh; ii++) {
            pBndFaceMesh= pMesh->getBndFaceMeshIX(ii);
            ofs << "BoundaryID= " << pBndFaceMesh->getID() << ", BndType= ";
            switch(pBndFaceMesh->getBndType()) {
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
            uiint numOfBndNode= pBndFaceMesh->getNumOfBNode();
            vdouble sumValue;
            sumValue.resize(pBndFaceMesh->getNumOfDOF());
            for(uiint idof=0; idof < pBndFaceMesh->getNumOfDOF(); idof++) sumValue[idof]=0.0;
            for(iii=0; iii < numOfBndNode; iii++) {
                pBndNode= pBndFaceMesh->getBNodeIX(iii);
                pNode= pBndNode->getNode();
                ofs << "BNodeID= " << pBndNode->getID() << ", NodeID= " << pNode->getID();
                uiint idof, numOfDOF= pBndFaceMesh->getNumOfDOF();
                for(idof=0; idof < numOfDOF; idof++) {
                    uiint dof = pBndFaceMesh->getDOF(idof);
                    ofs << ", dof= " << dof << ", Value[" << idof << "]= " << pBndNode->getValue(dof, mgLevel);
                    sumValue[idof] += pBndNode->getValue(dof, mgLevel);
                };
                ofs << endl;
            };
            if(pBndFaceMesh->getBndType()==pmw::BoundaryType::Neumann) {
                ofs << "BNode Value SumValue :";
                for(uiint idof=0; idof < pBndFaceMesh->getNumOfDOF(); idof++) {
                    ofs << "  sumValue[" << idof << "]=" << sumValue[idof];
                };
                ofs << endl;
            }
            if(pBndFaceMesh->getBndType()==pmw::BoundaryType::Neumann) {
                uiint numOfBndFace= pBndFaceMesh->getNumOfBFace();
                for(iii=0; iii < numOfBndFace; iii++) {
                    pBndFace= pBndFaceMesh->getBFaceIX(iii);
                    ofs << "BFaceID= " << pBndFace->getID();
                    uiint idof, numOfDOF= pBndFaceMesh->getNumOfDOF();
                    for(idof=0; idof < numOfDOF; idof++) {
                        uiint dof = pBndFaceMesh->getDOF(idof);
                        ofs << ", Value[" << idof << "]= " << pBndFace->getBndValue(dof);
                    };
                    ofs << endl;
                };
            }
        };
    };
}
