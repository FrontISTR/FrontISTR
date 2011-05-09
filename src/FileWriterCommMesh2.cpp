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
void CFileWriterCommMesh2::Write(ofstream& ofs, const uint& mgLevel)
{
   pmw::CAssyModel *pAssyModel;
   pmw::CMesh *pMesh;
   pmw::CCommMesh2 *pCommMesh2;

   pAssyModel= mpGMGModel->getAssyModel(mgLevel);

   uint numOfPart;
   numOfPart= pAssyModel->getNumOfMesh();

   uint imesh;
   //---
   //Mesh
   //---
   for(imesh=0; imesh< numOfPart; imesh++){
       pMesh= pAssyModel->getMesh(imesh);

       // Mesh ID (CommMesh2を所有しているMesh)
       ofs << " CommMesh2 in Mesh_ID:" << pMesh->getMeshID() << endl;

       uint numOfCommMesh2;
       numOfCommMesh2= pMesh->getCommMesh2Size();
       uint icomm;
       //---
       //CommMesh2
       //---
       for(icomm=0; icomm< numOfCommMesh2; icomm++){
           pCommMesh2= pMesh->getCommMesh2IX(icomm);

           // myRank(自身のランク) - TransmitRank(相手のランク)
           ofs << "CommMesh2 myRank:" << pCommMesh2->getRank() <<
                  ", CommMesh2 transmitRank:" << pCommMesh2->getTrasmitRank() << endl;

           uint numOfCommNode;
           uint icnode;
           numOfCommNode= pCommMesh2->getCommNodeSize();
           pmw::CCommNode *pCommNode;
           pmw::CNode     *pNode;

           cout << "FileWriterCommMesh2::Write,  numOfCommNode==" << numOfCommNode << endl;

           //---
           //CommNode
           //---
           for(icnode=0; icnode< numOfCommNode; icnode++){
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
           uint numOfCommFace;
           uint iface;
           numOfCommFace= pCommMesh2->getCommFaceSize();
           pmw::CCommFace *pCommFace;
           for(iface=0; iface < numOfCommFace; iface++){
               
               pCommFace= pCommMesh2->getCommFaceIX(iface);
               
               ofs << "CommFace ID= " << pCommFace->getID();
               
               //uint ivert;
               uint numOfCommNode= pCommFace->getCommNodeSize();
               
               for(icnode=0; icnode < numOfCommNode; icnode++){
                   pCommNode= pCommFace->getCommNode(icnode);
                   ofs << ", " << pCommNode->getID(); 
               };
               ofs << endl;

               for(icnode=0; icnode < numOfCommNode; icnode++){
                   pCommNode= pCommFace->getCommNode(icnode);
                   ofs << ":: X= " << pCommNode->getX() << ", Y= " << pCommNode->getY() << ", Z= " << pCommNode->getZ();
               };
               ofs << endl;
               
           };//comm_face ループエンド
           
       };//comm_mesh2 ループエンド
       
       ofs << " CommMesh2 END " << endl;
       
   };//mesh ループエンド
   
}



