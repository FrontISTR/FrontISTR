//
//  FileWriterBoundaryVolume.cpp
//
//              2010.05.07
//              k.Takeda
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


void CFileWriterBoundaryVolume::Write(ofstream& ofs, const uint& mgLevel)
{
   pmw::CAssyModel *pAssyModel;
   pmw::CMesh *pMesh;

   pmw::CBoundaryVolumeMesh *pBndVolumeMesh;
   pmw::CBoundaryVolume *pBndVolume;
   pmw::CBoundaryNode *pBndNode;
   pmw::CNode *pNode;

   pAssyModel= mpGMGModel->getAssyModel(mgLevel);

   uint i,ii,iii;
   uint numOfPart= pAssyModel->getNumOfMesh();

   for(i=0; i < numOfPart; i++){
       pMesh= pAssyModel->getMesh(i);

       uint numOfBndMesh= pMesh->getNumOfBoundaryVolumeMesh();

       for(ii=0; ii < numOfBndMesh; ii++){
           pBndVolumeMesh= pMesh->getBndVolumeMeshIX(ii);

           ofs << "BoundaryID= " << pBndVolumeMesh->getID()  << ", BndType= " ;
           switch(pBndVolumeMesh->getBndType()){
               case(pmw::BoundaryType::Dirichlet):
                   ofs << "Dirichlet";
                   break;
               case(pmw::BoundaryType::Neumann):
                   ofs << "Neumann";
                   break;
               default:
                   //TODO:pLogger
                   break;
           }
           ofs << endl;

           uint numOfBndNode= pBndVolumeMesh->getNumOfBNode();
           for(iii=0; iii < numOfBndNode; iii++){
               pBndNode= pBndVolumeMesh->getBNodeIX(iii);
               pNode= pBndNode->getNode();

               ofs << "BNodeID= " << pBndNode->getID() << ", NodeID= " << pNode->getID();

               uint idof, numOfDOF= pBndVolumeMesh->getNumOfDOF();
               for(idof=0; idof < numOfDOF; idof++){

                   uint dof = pBndVolumeMesh->getDOF(idof);

                   ofs << ", dof= " << dof << ", Value[" << idof << "]= " << pBndNode->getValue(dof, mgLevel);
               };
               ofs << endl;
           };

           uint numOfBndVolume= pBndVolumeMesh->getNumOfVolume();
           for(iii=0; iii < numOfBndVolume; iii++){
               pBndVolume= pBndVolumeMesh->getBVolumeIX(iii);

               ofs << "BVolumeID= " << pBndVolume->getID();

               uint idof, numOfDOF= pBndVolumeMesh->getNumOfDOF();
               for(idof=0; idof < numOfDOF; idof++){

                   uint dof = pBndVolumeMesh->getDOF(idof);

                   ofs << ", Value[" << idof << "]= " << pBndVolume->getBndValue(dof);
               };
               ofs << endl;
           };
       };
   };
}



