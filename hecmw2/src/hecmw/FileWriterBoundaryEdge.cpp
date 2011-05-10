//
//  FileWriterBoundaryEdge.cpp
//
//              2010.05.07
//              k.Takeda
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


void CFileWriterBoundaryEdge::Write(ofstream& ofs, const uint& mgLevel)
{
   pmw::CAssyModel *pAssyModel;
   pmw::CMesh *pMesh;

   pmw::CBoundaryEdgeMesh *pBndEdgeMesh;
   pmw::CBoundaryEdge *pBndEdge;
   pmw::CBoundaryNode *pBndNode;
   pmw::CNode *pNode;

   pAssyModel= mpGMGModel->getAssyModel(mgLevel);

   uint i,ii,iii;
   uint numOfPart= pAssyModel->getNumOfMesh();

   for(i=0; i < numOfPart; i++){
       pMesh= pAssyModel->getMesh(i);

       uint numOfBndMesh= pMesh->getNumOfBoundaryEdgeMesh();

       for(ii=0; ii < numOfBndMesh; ii++){
           pBndEdgeMesh= pMesh->getBndEdgeMeshIX(ii);

           ofs << "BoundaryID= " << pBndEdgeMesh->getID() << ", BndType= ";
           switch(pBndEdgeMesh->getBndType()){
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

           uint numOfBndNode= pBndEdgeMesh->getNumOfBNode();

           //Neumannの値、確認用
           vdouble sumValue;
           sumValue.resize(pBndEdgeMesh->getNumOfDOF());
           for(uint idof=0; idof < pBndEdgeMesh->getNumOfDOF(); idof++)  sumValue[idof]=0.0;

           for(iii=0; iii < numOfBndNode; iii++){
               pBndNode= pBndEdgeMesh->getBNodeIX(iii);
               pNode= pBndNode->getNode();

               ofs << "BNodeID= " << pBndNode->getID() << ", NodeID= " << pNode->getID();

               uint idof, numOfDOF= pBndEdgeMesh->getNumOfDOF();
               for(idof=0; idof < numOfDOF; idof++){

                   uint dof = pBndEdgeMesh->getDOF(idof);

                   ofs << ", dof= " << dof << ", Value[" << idof << "]= " << pBndNode->getValue(dof, mgLevel);

                   //Neumann Σ値
                   sumValue[idof] += pBndNode->getValue(dof, mgLevel);
               };
               ofs << endl;
           };
           //Neumann Σ値 確認
           if(pBndEdgeMesh->getBndType()==pmw::BoundaryType::Neumann){
               ofs << "BNode Value SumValue :";
               for(uint idof=0; idof < pBndEdgeMesh->getNumOfDOF(); idof++){
                   ofs << "  sumValue[" << idof << "]=" << sumValue[idof];
               };
               ofs << endl;
           }

           uint numOfBndEdge= pBndEdgeMesh->getNumOfEdge();
           for(iii=0; iii < numOfBndEdge; iii++){
               pBndEdge= pBndEdgeMesh->getBEdgeIX(iii);

               ofs << "BEdgeID= " << pBndEdge->getID();

               uint idof, numOfDOF= pBndEdgeMesh->getNumOfDOF();
               for(idof=0; idof < numOfDOF; idof++){

                   uint dof = pBndEdgeMesh->getDOF(idof);

                   ofs << ", Value[" << idof << "]= " << pBndEdge->getBndValue(dof);
               };
               ofs << endl;
           };
       };
   };
}






