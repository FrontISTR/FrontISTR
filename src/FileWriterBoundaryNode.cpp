//
//  FileWriterBoundaryNode.cpp
//
//
//
//                  2009.07.23
//                  2009.07.23
//                  k.Takeda
#include "FileWriterBoundaryNode.h"
using namespace FileIO;

CFileWriterBoundaryNode::CFileWriterBoundaryNode()
{
    ;
}

CFileWriterBoundaryNode::~CFileWriterBoundaryNode()
{
    ;
}


// method
// --
void CFileWriterBoundaryNode::Write(ofstream& ofs, const uint& mgLevel)
{
   pmw::CAssyModel *pAssyModel;
   pmw::CMesh *pMesh;
   pmw::CBoundaryNode *pBndNode;
   string white(" ");

   pAssyModel= mpGMGModel->getAssyModel(mgLevel);

   uint numOfPart= pAssyModel->getNumOfMesh();
   uint numOfBndNode,numOfValue;
   uint i,ii,iii;
   for(i=0; i< numOfPart; i++){
       pMesh= pAssyModel->getMesh(i);
       
       ofs << " -- BoundaryNode Block Start -- " << ", mgLevel == " << mgLevel<< ", Mesh ID==" << pMesh->getMeshID() << endl;

       numOfBndNode= pMesh->getNumOfBoundaryNode();
       for(ii=0; ii< numOfBndNode; ii++){
           pBndNode = pMesh->getBoundaryNode_withIndex(ii);

           ofs << pBndNode->getID() << white;

           numOfValue= pBndNode->CGeneralBoundary::numOfValue();

           for(iii=0; iii< numOfValue; iii++){
               ofs << pBndNode->CGeneralBoundary::getValue(iii) << white;
           };

           ofs << endl;
       };
       ofs << " -- BoundaryNode Block end -- " << endl;
   };
}



















