
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
void CFileWriterNode::WriteDebug(ofstream& ofs, const uiint& mgLevel)
{
   pmw::CAssyModel *pAssyModel;
   pmw::CMesh *pMesh;
   pmw::CNode *pNode;
   pmw::CCommMesh *pCommMesh;
   string white(" ");

   // AssyModelの取得(mgLevel <= MultiGrid Level)
   //
   pAssyModel= mpGMGModel->getAssyModel(mgLevel);
   
   
   uiint numOfPart= pAssyModel->getNumOfMesh();
   uiint numOfNode;
   uiint imesh,inode;
   ////debug
   //cout << "CFileWriteNode,  numOfPart == " << numOfPart << endl;

   pmw::CAggregateNode* pAggNode;
   vector<pmw::CAggregateNode*> vAggNode;
   pmw::CNode* pConnNode;
   uiint iagg, numOfAggNode;

   pmw::CIndexBucket* pBucket;

   //prolongater関連
   uiint numOfParent;//Refineで生成されたNodeの親Node数
   pmw::CNode* parentNode;//Refineで生成されたNodeの親Node
   uiint ipare;

   for(imesh=0; imesh< numOfPart; imesh++){
       //Meshの取得
       pMesh= pAssyModel->getMesh(imesh);
       pBucket= pMesh->getBucket();

       //AggregateNodeの取得
       if(mnSolutionType==pmw::SolutionType::FVM) vAggNode= pMesh->getAggNodes();

       ofs << " -- Node Block Start -- " << ", mgLevel == " << mgLevel<< ", Mesh ID==" << pMesh->getMeshID() << endl;

       numOfNode= pMesh->getNumOfNode();//CMesh::numOfNode
//       //debug
//       cout << "ノード数 => " << numOfNode << ", CFileWriterNode::Write" << endl;

       for(inode=0; inode< numOfNode; inode++){

           pNode= pMesh->getNodeIX(inode);

           //ofs << imesh << white
           ofs << white           // <= Visual確認のためimesh出力なし
               << pNode->getID() << white
               << pNode->getX()  << white << pNode->getY() << white << pNode->getZ();

           if(mnSolutionType==pmw::SolutionType::FVM){
               pAggNode= vAggNode[pNode->getID()];//NodeはID基準ー＞AggNode配列もIDに合わせてある

               ofs << white << "Connect_Nodes";
               //AggregateNodeの出力
               //
               numOfAggNode= pAggNode->getNumOfNode();
               for(iagg=0; iagg< numOfAggNode; iagg++){
                   pConnNode= pAggNode->getNode(iagg);

                   ofs << white << pConnNode->getID();
               };
           }
           
           ofs << white << "ParentNode";
           // prolongater出力
           //
           numOfParent= pNode->getNumOfParentNode();
           for(ipare=0; ipare< numOfParent; ipare++){
               parentNode= pNode->getParentNode(ipare);

               ofs << white << parentNode->getID();
           };
           ofs << endl;
           
       };//Nodeループエンド

       uiint numOfCommMesh= pMesh->getNumOfCommMesh();
       uiint icom, commNodeID;
       for(icom=0; icom< numOfCommMesh; icom++){
           //pCommMesh= pMesh->getCommMesh(icom);/////////<-CommIDが判明している場合は,CommIDを使用(Mesh内部でHashによりIndex変換される)
           pCommMesh= pMesh->getCommMeshIX(icom);

           uiint numOfSend= pCommMesh->getNumOfSendNode();
           uiint numOfRecv= pCommMesh->getNumOfRecvNode();
           uiint isend,irecv;
           ofs << " -- Send Node -- " << endl;
           for(isend=0; isend< numOfSend; isend++){
               pNode= pCommMesh->getSendNodeIX(isend);
               commNodeID= pCommMesh->getSendCommNodeID(isend);
               ofs << white << "Send Index" << isend << ", id= " << pNode->getID() << ", CommNodeID= " << commNodeID
                   << ", X= "<< pNode->getX() << ", Y= " << pNode->getY() << ", Z= " << pNode->getZ() << endl;

//               //CommNodeIDによるNode出力テスト
//               pNode= pCommMesh->getNode(commNodeID);
//               ofs << white << "Send Index" << isend << ", id= " << pNode->getID() << ", CommNodeIDによる出力テスト" << endl;
           };
           ofs << " -- Recv Node -- " << endl;
           for(irecv=0; irecv< numOfRecv; irecv++){
               pNode= pCommMesh->getRecvNodeIX(irecv);
               commNodeID= pCommMesh->getRecvCommNodeID(irecv);
               ofs << white << "Recv Index" << irecv << ", id= " << pNode->getID() << ", CommNodeID= " << commNodeID
                   << ", X= "<< pNode->getX() << ", Y= " << pNode->getY() << ", Z= " << pNode->getZ() << endl;

//               //CommNodeIDによるNode出力テスト
//               pNode= pCommMesh->getNode(commNodeID);
//               ofs << white << "Recv Index" << irecv << ", id= " << pNode->getID() << ", CommNodeIDによる出力テスト" << endl;
           };
       };
   };//Meshループエンド
   ofs << " -- Node Block End -- " << endl;
}




