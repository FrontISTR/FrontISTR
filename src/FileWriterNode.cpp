/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileWriterNode.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Mesh.h"
#include "FileWriterNode.h"
using namespace FileIO;
CFileWriterNode::CFileWriterNode()
{
    ;
}
CFileWriterNode::~CFileWriterNode()
{
    ;
}
void CFileWriterNode::Write(ofstream& ofs, const uint& mgLevel)
{
   pmw::CAssyModel *pAssyModel;
   pmw::CMesh *pMesh;
   pmw::CNode *pNode;
   pmw::CCommMesh *pCommMesh;
   string white(" ");
   pAssyModel= mpGMGModel->getAssyModel(mgLevel);
   uint numOfPart= pAssyModel->getNumOfMesh();
   uint numOfNode;
   uint imesh,inode;
   pmw::CAggregateNode* pAggNode;
   vector<pmw::CAggregateNode*> vAggNode;
   pmw::CNode* pConnNode;
   uint iagg, numOfAggNode;
   pmw::CIndexBucket* pBucket;
   uint numOfParent;
   pmw::CNode* parentNode;
   uint ipare;
   for(imesh=0; imesh< numOfPart; imesh++){
       pMesh= pAssyModel->getMesh(imesh);
       pBucket= pMesh->getBucket();
       vAggNode= pMesh->getAggNodes();
       ofs << " -- Node Block Start -- " << ", mgLevel == " << mgLevel<< ", Mesh ID==" << pMesh->getMeshID() << endl;
       numOfNode= pMesh->getNumOfNode();
       for(inode=0; inode< numOfNode; inode++){
           pNode= pMesh->getNodeIX(inode);
           pAggNode= vAggNode[pNode->getID()];
           ofs << white           
               << pNode->getID() << white
               << pNode->getX()  << white << pNode->getY() << white << pNode->getZ();
           ofs << white << "Connect_Nodes";
           numOfAggNode= pAggNode->getNumOfNode();
           for(iagg=0; iagg< numOfAggNode; iagg++){
               pConnNode= pAggNode->getNode(iagg);
               ofs << white << pConnNode->getID();
           };
           ofs << white << "ParentNode";
           numOfParent= pNode->getNumOfParentNode();
           for(ipare=0; ipare< numOfParent; ipare++){
               parentNode= pNode->getParentNode(ipare);
               ofs << white << parentNode->getID();
           };
           ofs << endl;
       };
       uint numOfCommMesh= pMesh->getNumOfCommMesh();
       uint icom, commNodeID;
       for(icom=0; icom< numOfCommMesh; icom++){
           pCommMesh= pMesh->getCommMeshIX(icom);
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
   };
   ofs << " -- Node Block End -- " << endl;
}
