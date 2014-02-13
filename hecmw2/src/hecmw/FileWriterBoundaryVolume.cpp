/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterBoundaryVolume.cpp
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
#include "BoundaryVolume.h"
#include "BoundaryVolumeMesh.h"
#include "FileWriterBoundaryVolume.h"
using namespace FileIO;
CFileWriterBoundaryVolume::CFileWriterBoundaryVolume()
{
    ;
}
CFileWriterBoundaryVolume::~CFileWriterBoundaryVolume()
{
    ;
}
void CFileWriterBoundaryVolume::WriteDebug(ofstream& ofs, const uiint& mgLevel)
{
    pmw::CAssyModel *pAssyModel;
    pmw::CMesh *pMesh;
    pmw::CBoundaryVolumeMesh *pBndVolumeMesh;
    pmw::CBoundaryVolume *pBndVolume;
    pmw::CBoundaryNode *pBndNode;
    pmw::CNode *pNode;
    pAssyModel= mpGMGModel->getAssyModel(mgLevel);
    uiint i,ii,iii;
    uiint numOfPart= pAssyModel->getNumOfMesh();
    for(i=0; i < numOfPart; i++) {
        pMesh= pAssyModel->getMesh(i);
        uiint numOfBndMesh= pMesh->getNumOfBoundaryVolumeMesh();
        for(ii=0; ii < numOfBndMesh; ii++) {
            pBndVolumeMesh= pMesh->getBndVolumeMeshIX(ii);
            ofs << "BoundaryID= " << pBndVolumeMesh->getID()  << ", BndType= " ;
            switch(pBndVolumeMesh->getBndType()) {
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
            vdouble sumValue;
            sumValue.resize(pBndVolumeMesh->getNumOfDOF());
            for(uiint idof=0; idof < pBndVolumeMesh->getNumOfDOF(); idof++)  sumValue[idof]=0.0;
            uiint numOfBndNode= pBndVolumeMesh->getNumOfBNode();
            for(iii=0; iii < numOfBndNode; iii++) {
                pBndNode= pBndVolumeMesh->getBNodeIX(iii);
                pNode= pBndNode->getNode();
                ofs << "BNodeID= " << pBndNode->getID() << ", NodeID= " << pNode->getID();
                uiint idof, numOfDOF= pBndVolumeMesh->getNumOfDOF();
                for(idof=0; idof < numOfDOF; idof++) {
                    uiint dof = pBndVolumeMesh->getDOF(idof);
                    ofs << ", dof= " << dof << ", Value[" << idof << "]= " << pBndNode->getValue(dof, mgLevel);
                    sumValue[idof] += pBndNode->getValue(dof,mgLevel);
                };
                ofs << endl;
            };
            if(pBndVolumeMesh->getBndType()==pmw::BoundaryType::Neumann) {
                ofs << "BNode Value SumValue :";
                for(uiint idof=0; idof < pBndVolumeMesh->getNumOfDOF(); idof++) {
                    ofs << "  sumValue[" << idof << "]=" << sumValue[idof];
                };
                ofs << endl;
            }
            if(pBndVolumeMesh->getBndType()==pmw::BoundaryType::Neumann) {
                uiint numOfBndVolume= pBndVolumeMesh->getNumOfVolume();
                for(iii=0; iii < numOfBndVolume; iii++) {
                    pBndVolume= pBndVolumeMesh->getBVolumeIX(iii);
                    ofs << "BVolumeID= " << pBndVolume->getID();
                    uiint idof, numOfDOF= pBndVolumeMesh->getNumOfDOF();
                    for(idof=0; idof < numOfDOF; idof++) {
                        uiint dof = pBndVolumeMesh->getDOF(idof);
                        ofs << ", Value[" << idof << "]= " << pBndVolume->getBndValue(dof);
                    };
                    ofs << endl;
                };
            }
        };
    };
}
