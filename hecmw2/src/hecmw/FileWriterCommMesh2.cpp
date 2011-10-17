/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileWriterCommMesh2.cpp
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "HEC_MPI.h"
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
void CFileWriterCommMesh2::WriteDebug(ofstream& ofs, const uiint& mgLevel)
{
   Utility::CLogger *pLogger= Utility::CLogger::Instance();
   pmw::CAssyModel *pAssyModel;
   pmw::CMesh *pMesh;
   pmw::CCommMesh2 *pCommMesh2;
   pAssyModel= mpGMGModel->getAssyModel(mgLevel);
   uiint nNumOfPart;
   nNumOfPart= pAssyModel->getNumOfMesh();
   uiint imesh;
   for(imesh=0; imesh< nNumOfPart; imesh++){
       pMesh= pAssyModel->getMesh(imesh);
       ofs << " CommMesh2 in Mesh_ID:" << pMesh->getMeshID() << endl;
       uiint nNumOfCommMesh2;
       nNumOfCommMesh2= pMesh->getCommMesh2Size();
       uiint icomm;
       for(icomm=0; icomm< nNumOfCommMesh2; icomm++){
           pCommMesh2= pMesh->getCommMesh2IX(icomm);
           ofs << "CommMesh2 myRank:" << pCommMesh2->getRank() <<
                  ", CommMesh2 transmitRank:" << pCommMesh2->getTrasmitRank() << endl;
           uiint nNumOfCommNode;
           uiint icnode;
           nNumOfCommNode= pCommMesh2->getCommNodeSize();
           pmw::CCommNode *pCommNode;
           pmw::CNode     *pNode;
           pLogger->Info(Utility::LoggerMode::Info, "CommMesh2::nNumOfCommNode ", nNumOfCommNode);
           for(icnode=0; icnode< nNumOfCommNode; icnode++){
               pCommNode = pCommMesh2->getCommNodeIX(icnode);
               pNode = pCommNode->getNode();
               ofs << " Node ID : " << pNode->getID()
                   << ", X= " << pNode->getX()
                   << ", Y= " << pNode->getY()
                   << ", Z= " << pNode->getZ() << endl;
           };
           uiint nNumOfCommFace;
           uiint iface;
           nNumOfCommFace= pCommMesh2->getCommFaceSize();
           pmw::CCommFace *pCommFace;
           for(iface=0; iface < nNumOfCommFace; iface++){
               pCommFace= pCommMesh2->getCommFaceIX(iface);
               ofs << "CommFace ID= " << pCommFace->getID();
               uiint nNumOfCommNodeFace= pCommFace->getCommNodeSize();
               for(icnode=0; icnode < nNumOfCommNodeFace; icnode++){
                   pCommNode= pCommFace->getCommNode(icnode);
                   ofs << ", " << pCommNode->getID(); 
               };
               ofs << endl;
               for(icnode=0; icnode < nNumOfCommNodeFace; icnode++){
                   pCommNode= pCommFace->getCommNode(icnode);
                   ofs << ":: X= " << pCommNode->getX() << ", Y= " << pCommNode->getY() << ", Z= " << pCommNode->getZ();
               };
               ofs << endl;
           };
       };
       ofs << " CommMesh2 END " << endl;
   };
}
