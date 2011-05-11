
#include "BoundaryMesh.h"

//
//  FileWriterBoundaryFace.cpp
//
//              2010.05.07
//              k.Takeda
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

   for(i=0; i < numOfPart; i++){
       pMesh= pAssyModel->getMesh(i);
       
       uiint numOfBndMesh= pMesh->getNumOfBoundaryFaceMesh();
       
       for(ii=0; ii < numOfBndMesh; ii++){
           pBndFaceMesh= pMesh->getBndFaceMeshIX(ii);

           ofs << "BoundaryID= " << pBndFaceMesh->getID() << ", BndType= ";
           switch(pBndFaceMesh->getBndType()){
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

           uiint numOfBndNode= pBndFaceMesh->getNumOfBNode();
           //cout << "FileWriterBoundaryFace::Write, numOfBndNode = " << numOfBndNode << endl;

           //Neumannの値、確認用
           vdouble sumValue;
           sumValue.resize(pBndFaceMesh->getNumOfDOF());
           for(uiint idof=0; idof < pBndFaceMesh->getNumOfDOF(); idof++) sumValue[idof]=0.0;

           for(iii=0; iii < numOfBndNode; iii++){

               pBndNode= pBndFaceMesh->getBNodeIX(iii);
               pNode= pBndNode->getNode();

               //cout << "FileWriterBoundaryFace::Write, BNodeID " << pBndNode->getID() << endl;
               //cout << "FileWriterBoundaryFace::Write,  NodeID " << pNode->getID() << endl;
               
               ofs << "BNodeID= " << pBndNode->getID() << ", NodeID= " << pNode->getID();

               uiint idof, numOfDOF= pBndFaceMesh->getNumOfDOF();
               for(idof=0; idof < numOfDOF; idof++){

                   uiint dof = pBndFaceMesh->getDOF(idof);

                   ofs << ", dof= " << dof << ", Value[" << idof << "]= " << pBndNode->getValue(dof, mgLevel);

                   //Neumann Σ値
                   sumValue[idof] += pBndNode->getValue(dof, mgLevel);
               };
               ofs << endl;
           };
           //Neumann Σ値 確認
           if(pBndFaceMesh->getBndType()==pmw::BoundaryType::Neumann){
               ofs << "BNode Value SumValue :";
               for(uiint idof=0; idof < pBndFaceMesh->getNumOfDOF(); idof++){
                   ofs << "  sumValue[" << idof << "]=" << sumValue[idof];
               };
               ofs << endl;
           }

           //Neumann Face境界値 確認
           if(pBndFaceMesh->getBndType()==pmw::BoundaryType::Neumann){
               uiint numOfBndFace= pBndFaceMesh->getNumOfBFace();
               for(iii=0; iii < numOfBndFace; iii++){
                   pBndFace= pBndFaceMesh->getBFaceIX(iii);

                   ofs << "BFaceID= " << pBndFace->getID();

                   uiint idof, numOfDOF= pBndFaceMesh->getNumOfDOF();
                   for(idof=0; idof < numOfDOF; idof++){

                       uiint dof = pBndFaceMesh->getDOF(idof);

                       ofs << ", Value[" << idof << "]= " << pBndFace->getBndValue(dof);
                   };
                   ofs << endl;
               };
           }
       };
   };
}





