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

   for(i=0; i < numOfPart; i++){
       pMesh= pAssyModel->getMesh(i);

       uiint numOfBndMesh= pMesh->getNumOfBoundaryEdgeMesh();

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

           uiint numOfBndNode= pBndEdgeMesh->getNumOfBNode();

           //Neumannの値、確認用
           vdouble sumValue;
           sumValue.resize(pBndEdgeMesh->getNumOfDOF());
           for(uiint idof=0; idof < pBndEdgeMesh->getNumOfDOF(); idof++)  sumValue[idof]=0.0;

           for(iii=0; iii < numOfBndNode; iii++){
               pBndNode= pBndEdgeMesh->getBNodeIX(iii);
               pNode= pBndNode->getNode();

               ofs << "BNodeID= " << pBndNode->getID() << ", NodeID= " << pNode->getID();

               uiint idof, numOfDOF= pBndEdgeMesh->getNumOfDOF();
               for(idof=0; idof < numOfDOF; idof++){

                   uiint dof = pBndEdgeMesh->getDOF(idof);

                   ofs << ", dof= " << dof << ", Value[" << idof << "]= " << pBndNode->getValue(dof, mgLevel);

                   //Neumann Σ値
                   sumValue[idof] += pBndNode->getValue(dof, mgLevel);
               };
               ofs << endl;
           };
           //Neumann Σ値 確認
           if(pBndEdgeMesh->getBndType()==pmw::BoundaryType::Neumann){
               ofs << "BNode Value SumValue :";
               for(uiint idof=0; idof < pBndEdgeMesh->getNumOfDOF(); idof++){
                   ofs << "  sumValue[" << idof << "]=" << sumValue[idof];
               };
               ofs << endl;
           }

           //Neumann Edge境界値 確認
           if(pBndEdgeMesh->getBndType()==pmw::BoundaryType::Neumann){
               uiint numOfBndEdge= pBndEdgeMesh->getNumOfEdge();
               for(iii=0; iii < numOfBndEdge; iii++){
                   pBndEdge= pBndEdgeMesh->getBEdgeIX(iii);

                   ofs << "BEdgeID= " << pBndEdge->getID();

                   uiint idof, numOfDOF= pBndEdgeMesh->getNumOfDOF();
                   for(idof=0; idof < numOfDOF; idof++){

                       uiint dof = pBndEdgeMesh->getDOF(idof);

                       ofs << ", Value[" << idof << "]= " << pBndEdge->getBndValue(dof);
                   };
                   ofs << endl;
               };
           }
       };
   };
}






