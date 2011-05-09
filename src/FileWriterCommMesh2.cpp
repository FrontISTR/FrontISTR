/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileWriterCommMesh2.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FileWriterCommMesh2.h"
using namespace FileIO;
CFileWriterCommMesh2::CFileWriterCommMesh2()
{
    ;
}
CFileWriterCommMesh2::~CFileWriterCommMesh2()
{
    ;
}
void CFileWriterCommMesh2::Write(ofstream& ofs, const uint& mgLevel)
{
   pmw::CAssyModel *pAssyModel;
   pmw::CMesh *pMesh;
   pmw::CCommMesh2 *pCommMesh2;
   pAssyModel= mpGMGModel->getAssyModel(mgLevel);
   uint numOfPart;
   numOfPart= pAssyModel->getNumOfMesh();
   uint imesh;
   for(imesh=0; imesh< numOfPart; imesh++){
       pMesh= pAssyModel->getMesh(imesh);
       ofs << " CommMesh2 in Mesh_ID:" << pMesh->getMeshID() << endl;
       uint numOfCommMesh2;
       numOfCommMesh2= pMesh->getCommMesh2Size();
       uint icomm;
       for(icomm=0; icomm< numOfCommMesh2; icomm++){
           pCommMesh2= pMesh->getCommMesh2IX(icomm);
           ofs << "CommMesh2 myRank:" << pCommMesh2->getRank() <<
                  ", CommMesh2 transmitRank:" << pCommMesh2->getTrasmitRank() << endl;
           uint numOfCommNode;
           uint icnode;
           numOfCommNode= pCommMesh2->getCommVertNodeSize();
           pmw::CCommNode *pCommNode;
           pmw::CNode     *pNode;
           for(icnode=0; icnode< numOfCommNode; icnode++){
               pCommNode= pCommMesh2->getCommVertNodeIX(icnode);
               pNode= pCommNode->getNode();
               ofs << " Node ID : " << pNode->getID()
                   << ", X= " << pNode->getX()
                   << ", Y= " << pNode->getY()
                   << ", Z= " << pNode->getZ() << endl;
           };
           uint numOfCommFace;
           uint iface;
           numOfCommFace= pCommMesh2->getCommFaceSize();
           pmw::CCommFace *pCommFace;
           for(iface=0; iface < numOfCommFace; iface++){
               pCommFace= pCommMesh2->getCommFaceIX(iface);
               ofs << "CommFace ID= " << pCommFace->getID();
               uint ivert;
               uint numOfVert= pCommFace->getVertCommNodeSize();
               for(ivert=0; ivert < numOfVert; ivert++){
                   pCommNode= pCommFace->getVertCommNode(ivert);
                   ofs << ", " << pCommNode->getID();
               };
               ofs << endl;
               for(ivert=0; ivert < numOfVert; ivert++){
                   pCommNode= pCommFace->getVertCommNode(ivert);
                   ofs << ":: X= " << pCommNode->getX() << ", Y= " << pCommNode->getY() << ", Z= " << pCommNode->getZ();
               };
               ofs << endl;
           };
       };
       ofs << " CommMesh2 END " << endl;
   };
}
