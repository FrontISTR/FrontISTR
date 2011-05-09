/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileWriterContactMesh.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "SkinFace.h"
#include "ContactNode.h"
#include "ContactMesh.h"
#include "AssyModel.h"
#include "FileWriterContactMesh.h"
using namespace FileIO;
CFileWriterContactMesh::CFileWriterContactMesh()
{
    ;
}
CFileWriterContactMesh::~CFileWriterContactMesh()
{
    ;
}
void CFileWriterContactMesh::Write(ofstream& ofs, const uint& mgLevel)
{
    pmw::CAssyModel *pAssy;
    pAssy= mpGMGModel->getAssyModel(mgLevel);
    string white(" ");
    pmw::CContactMesh *pConMesh;
    uint numOfContact;
    uint icont;
    pmw::CSkinFace *pSkinFace;
    uint numOfSkinFace;
    uint iface;
    pmw::CContactNode *pConNode;
    uint numOfConNode;
    uint icnode;
    uint maslave;
    numOfContact= pAssy->getNumOfContactMesh();
    for(icont=0; icont < numOfContact; icont++){
        pConMesh= pAssy->getContactMesh(icont);
        for(maslave=0; maslave< 2; maslave++){
            if(maslave==0) numOfSkinFace= pConMesh->getNumOfMasterFace();
            if(maslave==1) numOfSkinFace= pConMesh->getNumOfSlaveFace();
            for(iface=0; iface< numOfSkinFace; iface++){
                if(maslave==0) pSkinFace= pConMesh->getMasterFace(iface);
                if(maslave==1) pSkinFace= pConMesh->getSlaveFace(iface);
                numOfConNode= pSkinFace->getNumOfNode();
                if(maslave==0) ofs << white << " MasterFaceID=" << pSkinFace->getID();
                if(maslave==1) ofs << white << " SlaveFaceID =" << pSkinFace->getID();
                ofs << ", ConNode:";
                for(icnode=0; icnode< numOfConNode; icnode++){
                    pConNode= pSkinFace->getNode(icnode);
                    ofs << white << pConNode->getID();
                };
                ofs << ", Node:";
                for(icnode=0; icnode< numOfConNode; icnode++){
                    pConNode= pSkinFace->getNode(icnode);
                    if(pSkinFace->isSelf()){
                        ofs << " ID= " << pConNode->getNodeID() << ", MeshID= " << pConNode->getMeshID() << ",";
                    }else{
                        ofs << " ID= " << "-"                   << ", MeshID= " << "-"                   << ",";
                    }
                };
                ofs << ", rank:";
                for(icnode=0; icnode< numOfConNode; icnode++){
                    pConNode= pSkinFace->getNode(icnode);
                    ofs << white << pConNode->getRank();
                };
                ofs << endl;
            };
        };
    };
    numOfContact= pAssy->getNumOfContactMesh();
    for(icont=0; icont < numOfContact; icont++){
        pConMesh= pAssy->getContactMesh(icont);
        numOfSkinFace= pConMesh->getNumOfMasterFace();
        for(iface=0; iface< numOfSkinFace; iface++){
            pSkinFace= pConMesh->getMasterFace(iface);
            ofs << " MasterFaceID= " << pSkinFace->getID();
            ofs << ", SlaveP::";
            uint numOfSlaveP= pSkinFace->getNumOfSlaveNode();
            uint isnode;
            for(isnode=0; isnode < numOfSlaveP; isnode++){
                pConNode= pSkinFace->getSlaveNode(isnode);
                ofs << " ConNodeID= " << pConNode->getID() << ",";
            };
            ofs << endl;
        };
    };
    uint numOfSPoint;
    uint islave;
    int  masterFaceID;
    pmw::CContactNode* pSlaveNode;
    pmw::CSkinFace* pMasterFace;
    pmw::CContactNode* pMasterNode;
    pmw::CHecMPI *pMPI= pmw::CHecMPI::Instance();
    MPI_Status stat;
    MPI_Request request;
    ofs << " -- Coef出力チェック -- " << endl;
    for(icont=0; icont < numOfContact; icont++){
        pConMesh= pAssy->getContactMesh(icont);
        numOfSPoint= pConMesh->getNumOfSlavePoint();
        for(islave=0; islave< numOfSPoint; islave++){
            pSlaveNode = pConMesh->getSlaveConNode(islave);
            ofs  << " スレーブ点 ConID= " << pSlaveNode->getID() ;
            int my_rank;
            pMPI->Comm_rank(MPI_COMM_WORLD, &my_rank);
            if(pSlaveNode->getRank()== my_rank){
                ofs << ", Slave Node_ID= " << pSlaveNode->getNodeID() << endl;
            }else{
                ofs << ", Slave Node_ID= " << " - " << endl;
            }
            if(pSlaveNode->have_MasterFaceID(mgLevel)){
                masterFaceID= pSlaveNode->getMasterFaceID(mgLevel);
                pMasterFace= pConMesh->getMasterFace_ID(masterFaceID);
                ofs << " マスター面 ID= " << masterFaceID << endl;
                int ivert, numOfVert;
                numOfVert= (int)pMasterFace->getNumOfNode();
                for(ivert=0; ivert< numOfVert; ivert++){
                    ofs <<  " coef= " << pMasterFace->getCoef(pSlaveNode->getID(),ivert);
                };
                ofs << endl;
                int qnode_id[4];
                uint rank_count(0);
                for(ivert=0; ivert< 4; ivert++){
                    pMasterNode= pMasterFace->getNode(ivert);
                    qnode_id[ivert]= pMasterNode->getNodeID();
                    if(pMasterNode->getRank()==pConMesh->getRank()) rank_count++;
                };
				if(pConMesh->getRank() != pConMesh->getTransmitRank()){
					if(rank_count==numOfVert){
						int transRank= (int)pConMesh->getTransmitRank();
						pMPI->Send(qnode_id,numOfVert,MPI_INT,transRank,1,MPI_COMM_WORLD);
					}
					if(rank_count==0){
						int rnode_id[4];
						int recvRank= pConMesh->getTransmitRank();
						pMPI->Recv(rnode_id,numOfVert,MPI_INT,recvRank,1,MPI_COMM_WORLD,&stat);
					}
				}
            }
        };
        ofs << " -- Coef出力  end  -- " << endl;
    };
}
