
#include "Mesh.h"

//
//  FileWriter.cpp
//
//
//
//                      2009.07.23
//                      2009.07.23
//                      k.Takeda
#include "FileWriterNode.h"
using namespace FileIO;


// construct & destruct
//
CFileWriterNode::CFileWriterNode()
{
    ;
}
CFileWriterNode::~CFileWriterNode()
{
    ;
}

// method
// --
void CFileWriterNode::Write(ofstream& ofs, const uint& mgLevel)
{
   pmw::CAssyModel *pAssyModel;
   pmw::CMesh *pMesh;
   pmw::CNode *pNode;
   pmw::CCommMesh *pCommMesh;
   string white(" ");

   // AssyModelの取得(mgLevel <= MultiGrid Level)
   //
   pAssyModel= mpGMGModel->getAssyModel(mgLevel);
   
   
   uint numOfPart= pAssyModel->getNumOfMesh();
   uint numOfNode;
   uint imesh,inode;
   ////debug
   //cout << "CFileWriteNode,  numOfPart == " << numOfPart << endl;

   pmw::CAggregateNode* pAggNode;
   vector<pmw::CAggregateNode*> vAggNode;
   pmw::CNode* pConnNode;
   uint iagg, numOfAggNode;

   pmw::CIndexBucket* pBucket;
   

   for(imesh=0; imesh< numOfPart; imesh++){
       //Meshの取得
       pMesh= pAssyModel->getMesh(imesh);
       pBucket= pMesh->getBucket();

       //AggregateNodeの取得
       vAggNode= pMesh->getAggNodes();

       ofs << " -- Node Block Start -- " << ", mgLevel == " << mgLevel<< ", Mesh ID==" << pMesh->getMeshID() << endl;

       numOfNode= pMesh->getNumOfNode();//CMesh::numOfNode
       ////debug
       //cout << "ノード数 => " << numOfNode << ", CFileWriterNode::Write" << endl;

       for(inode=0; inode< numOfNode; inode++){

           pNode= pMesh->getNodeIX(inode);

           pAggNode= vAggNode[pNode->getID()];//NodeはID基準ー＞AggNode配列もIDに合わせてある

           //ofs << imesh << white
           ofs << white           // <= Visual確認のためimesh出力なし
               << pNode->getID() << white
               << pNode->getX()  << white << pNode->getY() << white << pNode->getZ();

           ofs << white << "Connect_Nodes";
           //AggregateNodeの出力
           //
           numOfAggNode= pAggNode->getNumOfNode();
           for(iagg=0; iagg< numOfAggNode; iagg++){
               pConnNode= pAggNode->getNode(iagg);

               ofs << white << pConnNode->getID();
           };
           ofs << endl;
       };//Nodeループエンド

       uint numOfCommMesh= pMesh->getNumOfCommMesh();
       uint icom, commNodeID;
       for(icom=0; icom< numOfCommMesh; icom++){
           pCommMesh= pMesh->getCommMesh(icom);/////////<-本来はCommID 修正の必要あり.

           uint numOfSend= pCommMesh->getNumOfSendNode();
           uint numOfRecv= pCommMesh->getNumOfRecvNode();
           uint isend,irecv;
           ofs << " -- Send Node -- " << endl;
           for(isend=0; isend< numOfSend; isend++){
               pNode= pCommMesh->getSendNodeIX(isend);
               commNodeID= pCommMesh->getSendCommNodeID(isend);
               ofs << white << "Send Index" << isend << ", id= " << pNode->getID() << ", CommNodeID= " << commNodeID
                   << ", X= "<< pNode->getX() << ", Y= " << pNode->getY() << ", Z= " << pNode->getZ() << endl;
           };
           ofs << " -- Recv Node -- " << endl;
           for(irecv=0; irecv< numOfRecv; irecv++){
               pNode= pCommMesh->getRecvNodeIX(irecv);
               commNodeID= pCommMesh->getRecvCommNodeID(irecv);
               ofs << white << "Recv Index" << irecv << ", id= " << pNode->getID() << ", CommNodeID= " << commNodeID
                   << ", X= "<< pNode->getX() << ", Y= " << pNode->getY() << ", Z= " << pNode->getZ() << endl;
           };
       };
   };//Meshループエンド
   ofs << " -- Node Block End -- " << endl;
}




