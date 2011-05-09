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

   for(imesh=0; imesh< numOfPart; imesh++){
       //Meshの取得
       pMesh= pAssyModel->getMesh(imesh);

       //AggregateNodeの取得
       vAggNode= pMesh->getAggNodes();

       ofs << " -- Node Block Start -- " << ", mgLevel == " << mgLevel<< ", Mesh ID==" << pMesh->getMeshID() << endl;

       numOfNode= pMesh->getNumOfNode();//CMesh::numOfNode
       ////debug
       //cout << "ノード数 => " << numOfNode << ", CFileWriterNode::Write" << endl;

       for(inode=0; inode< numOfNode; inode++){
           pNode= pMesh->getNode(inode);
           pAggNode= vAggNode[inode];

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
       };
   };
   ofs << " -- Node Block End -- " << endl;
}




