/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterNode.cpp
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
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
void CFileWriterNode::WriteDebug(ofstream& ofs, const uiint& mgLevel)
{
    pmw::CAssyModel *pAssyModel;
    pmw::CMesh *pMesh;
    pmw::CNode *pNode;
    pmw::CCommMesh *pCommMesh;
    string white(" ");
    pAssyModel= mpGMGModel->getAssyModel(mgLevel);
    uiint nNumOfPart= pAssyModel->getNumOfMesh();
    uiint nNumOfNode;
    uiint imesh,inode;
    pmw::CAggregateNode* pAggNode;
    vector<pmw::CAggregateNode*> vAggNode;
    pmw::CNode* pConnNode;
    uiint iagg, numOfAggNode;
    pmw::CIndexBucket* pBucket;
    uiint nNumOfParent;
    pmw::CNode* parentNode;
    uiint ipare;

    for(imesh=0; imesh< nNumOfPart; imesh++) {
        pMesh= pAssyModel->getMesh(imesh);
        pBucket= pMesh->getBucket();
        if(mnSolutionType==pmw::SolutionType::FVM) vAggNode= pMesh->getAggNodes();
        ofs << " -- Node Block Start -- " << ", mgLevel == " << mgLevel<< ", Mesh ID==" << pMesh->getMeshID() << endl;
        nNumOfNode= pMesh->getNumOfNode();
        for(inode=0; inode< nNumOfNode; inode++) {
            pNode= pMesh->getNodeIX(inode);
            ofs << white
                << pNode->getID() << white
                << pNode->getX()  << white << pNode->getY() << white << pNode->getZ();
            if(mnSolutionType==pmw::SolutionType::FVM) {
                pAggNode= vAggNode[pNode->getID()];
                ofs << white << "Connect_Nodes";
                numOfAggNode= pAggNode->getNumOfNode();
                for(iagg=0; iagg< numOfAggNode; iagg++) {
                    pConnNode= pAggNode->getNode(iagg);
                    ofs << white << pConnNode->getID();
                };
            }
            ofs << white << "ParentNode";
            nNumOfParent= pNode->getNumOfParentNode(mgLevel);

            for(ipare=0; ipare< nNumOfParent; ipare++) {
                parentNode= pNode->getParentNode(mgLevel, ipare);
                ofs << white << parentNode->getID();
            };
            ofs << endl;
        };
        uiint numOfCommMesh= pMesh->getNumOfCommMesh();
        uiint icom, commNodeID;
        for(icom=0; icom< numOfCommMesh; icom++) {
            pCommMesh= pMesh->getCommMeshIX(icom);
            uiint numOfSend= pCommMesh->getNumOfSendNode();
            uiint numOfRecv= pCommMesh->getNumOfRecvNode();
            uiint isend,irecv;
            ofs << " -- Send Node -- " << endl;
            for(isend=0; isend< numOfSend; isend++) {
                pNode= pCommMesh->getSendNodeIX(isend);
                commNodeID= pCommMesh->getSendCommNodeID(isend);
                ofs << white << "Send Index" << isend << ", id= " << pNode->getID() << ", CommNodeID= " << commNodeID
                    << ", X= "<< pNode->getX() << ", Y= " << pNode->getY() << ", Z= " << pNode->getZ() << endl;
            };
            ofs << " -- Recv Node -- " << endl;
            for(irecv=0; irecv< numOfRecv; irecv++) {
                pNode= pCommMesh->getRecvNodeIX(irecv);
                commNodeID= pCommMesh->getRecvCommNodeID(irecv);
                ofs << white << "Recv Index" << irecv << ", id= " << pNode->getID() << ", CommNodeID= " << commNodeID
                    << ", X= "<< pNode->getX() << ", Y= " << pNode->getY() << ", Z= " << pNode->getZ() << endl;
            };
        };
        // BNode Marking
        ofs << " -- isBNode -- " << endl;
        for(inode=0; inode< nNumOfNode; inode++) {
            bool bMarking = pMesh->isBNode(inode);
            pNode= pMesh->getNodeIX(inode);
            uiint id = pNode->getID();
            uiint index= pBucket->getIndexNode(id);

            if(index!=inode) cout << "Error, pBucket_Index mismatch" << endl;

            ofs << "i:" << index << " id:" << id << " " << bMarking << endl;
        };
        // Dirichlet Marking
        ofs << " -- isDirichlet -- " << endl;
        for(inode=0; inode< nNumOfNode; inode++) {
            bool bMarking = pMesh->isDirichletBNode(inode);
            pNode= pMesh->getNodeIX(inode);
            uiint id = pNode->getID();
            uiint index= pBucket->getIndexNode(id);

            if(index!=inode) cout << "Error, pBucket_Index mismatch" << endl;

            ofs << "i:" << index << " id:" << id << " " << bMarking << endl;
        };
        // Neumann Marking
        ofs << " -- isNeumann -- " << endl;
        for(inode=0; inode< nNumOfNode; inode++) {
            bool bMarking = pMesh->isNeumannBNode(inode);
            pNode= pMesh->getNodeIX(inode);
            uiint id = pNode->getID();
            uiint index= pBucket->getIndexNode(id);

            if(index!=inode) cout << "Error, pBucket_Index mismatch" << endl;

            ofs << "i:" << index << " id:" << id << " " << bMarking << endl;
        };
    };

    ofs << " -- Node Block End -- " << endl;
}
