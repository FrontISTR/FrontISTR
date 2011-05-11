//
//  FileWriterCommMesh2.cpp
//
//
//
//
//              2010.03.15
//              k.Takeda
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



// 要素分割型(節点共有型)
//
//  通信テーブル出力サンプル
//
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
   //---
   //Mesh
   //---
   for(imesh=0; imesh< nNumOfPart; imesh++){
       pMesh= pAssyModel->getMesh(imesh);

       // Mesh ID (CommMesh2を所有しているMesh)
       ofs << " CommMesh2 in Mesh_ID:" << pMesh->getMeshID() << endl;

       uiint nNumOfCommMesh2;
       nNumOfCommMesh2= pMesh->getCommMesh2Size();
       uiint icomm;
       //---
       //CommMesh2
       //---
       for(icomm=0; icomm< nNumOfCommMesh2; icomm++){
           pCommMesh2= pMesh->getCommMesh2IX(icomm);

           // myRank(自身のランク) - TransmitRank(相手のランク)
           ofs << "CommMesh2 myRank:" << pCommMesh2->getRank() <<
                  ", CommMesh2 transmitRank:" << pCommMesh2->getTrasmitRank() << endl;

           uiint nNumOfCommNode;
           uiint icnode;
           nNumOfCommNode= pCommMesh2->getCommNodeSize();
           pmw::CCommNode *pCommNode;
           pmw::CNode     *pNode;

           pLogger->Info(Utility::LoggerMode::Info, "CommMesh2::nNumOfCommNode ", nNumOfCommNode);
           
           //---
           //CommNode
           //---
           for(icnode=0; icnode< nNumOfCommNode; icnode++){
               // 通信するNodeID
               pCommNode = pCommMesh2->getCommNodeIX(icnode);
               pNode = pCommNode->getNode();

               ofs << " Node ID : " << pNode->getID()
                   << ", X= " << pNode->getX()
                   << ", Y= " << pNode->getY()
                   << ", Z= " << pNode->getZ() << endl;
           };//comm_node ループエンド
           
           //---
           //CommFace
           //---
           uiint nNumOfCommFace;
           uiint iface;
           nNumOfCommFace= pCommMesh2->getCommFaceSize();
           pmw::CCommFace *pCommFace;
           for(iface=0; iface < nNumOfCommFace; iface++){
               
               pCommFace= pCommMesh2->getCommFaceIX(iface);
               
               ofs << "CommFace ID= " << pCommFace->getID();
               
               //uint ivert;
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
               
           };//comm_face ループエンド
           
       };//comm_mesh2 ループエンド
       
       ofs << " CommMesh2 END " << endl;
       
   };//mesh ループエンド
   
}



