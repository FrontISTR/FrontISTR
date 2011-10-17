/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileWriterNode.cpp
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
void CFileWriterNode::WriteDebug(ofstream& ofs, const uiint& mgLevel)
{
   pmw::CAssyModel *pAssyModel;
   pmw::CMesh *pMesh;
   pmw::CNode *pNode;
   pmw::CCommMesh *pCommMesh;
   string white(" ");
   pAssyModel= mpGMGModel->getAssyModel(mgLevel);
   uiint numOfPart= pAssyModel->getNumOfMesh();
   uiint numOfNode;
   uiint imesh,inode;
   pmw::CAggregateNode* pAggNode;
   vector<pmw::CAggregateNode*> vAggNode;
   pmw::CNode* pConnNode;
   uiint iagg, numOfAggNode;
   pmw::CIndexBucket* pBucket;
   uiint numOfParent;
   pmw::CNode* parentNode;
   uiint ipare;
   for(imesh=0; imesh< numOfPart; imesh++){
       pMesh= pAssyModel->getMesh(imesh);
       pBucket= pMesh->getBucket();
       if(mnSolutionType==pmw::SolutionType::FVM) vAggNode= pMesh->getAggNodes();
       ofs << " -- Node Block Start -- " << ", mgLevel == " << mgLevel<< ", Mesh ID==" << pMesh->getMeshID() << endl;
       numOfNode= pMesh->getNumOfNode();
       for(inode=0; inode< numOfNode; inode++){
           pNode= pMesh->getNodeIX(inode);
           ofs << white           
               << pNode->getID() << white
               << pNode->getX()  << white << pNode->getY() << white << pNode->getZ();
           if(mnSolutionType==pmw::SolutionType::FVM){
               pAggNode= vAggNode[pNode->getID()];
               ofs << white << "Connect_Nodes";
               numOfAggNode= pAggNode->getNumOfNode();
               for(iagg=0; iagg< numOfAggNode; iagg++){
                   pConnNode= pAggNode->getNode(iagg);
                   ofs << white << pConnNode->getID();
               };
           }
           ofs << white << "ParentNode";
           numOfParent= pNode->getNumOfParentNode();
           for(ipare=0; ipare< numOfParent; ipare++){
               parentNode= pNode->getParentNode(ipare);
               ofs << white << parentNode->getID();
           };
           ofs << endl;
       };
       uiint numOfCommMesh= pMesh->getNumOfCommMesh();
       uiint icom, commNodeID;
       for(icom=0; icom< numOfCommMesh; icom++){
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
