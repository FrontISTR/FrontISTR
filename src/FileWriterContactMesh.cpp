//
//  FileWriterContactMesh.cpp
//
//
//
//                  2009.11.20
//                  2009.11.20
//                  k.Takeda
#include "SkinFace.h"
#include "ContactNode.h"
#include "ContactMesh.h"
#include "AssyModel.h"

#include "FileWriterContactMesh.h"
using namespace FileIO;

// construct & destruct
//
CFileWriterContactMesh::CFileWriterContactMesh()
{
    ;
}
CFileWriterContactMesh::~CFileWriterContactMesh()
{
    ;
}

// MW3標準のファイル出力
//
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

    uint maslave;//マスター,スレーブ切り替え

    // ----
    // 面構成Nodeの出力,Nodeのrank出力
    // ----
    numOfContact= pAssy->getNumOfContactMesh();
    // ContactMeshループ
    // --
    for(icont=0; icont < numOfContact; icont++){
        pConMesh= pAssy->getContactMesh(icont);
        // マスター,スレーブ切り替え
        // --
        for(maslave=0; maslave< 2; maslave++){
            if(maslave==0) numOfSkinFace= pConMesh->getNumOfMasterFace();//マスター面数
            if(maslave==1) numOfSkinFace= pConMesh->getNumOfSlaveFace();//スレーブ面数
            // SkinFaceループ
            // --
            for(iface=0; iface< numOfSkinFace; iface++){
                if(maslave==0) pSkinFace= pConMesh->getMasterFace(iface);//マスター面
                if(maslave==1) pSkinFace= pConMesh->getSlaveFace(iface);//スレーブ面

                numOfConNode= pSkinFace->getNumOfNode();
                // ファイル出力
                // --
                // SkinFaceID  ContactNodeID .........
                if(maslave==0) ofs << white << " MasterFaceID=" << pSkinFace->getID();
                if(maslave==1) ofs << white << " SlaveFaceID =" << pSkinFace->getID();

                ofs << ", ConNode:";
                for(icnode=0; icnode< numOfConNode; icnode++){
                    pConNode= pSkinFace->getNode(icnode);
                    
                    ofs << white << pConNode->getID();
                    
                };//icnodeループ
                ofs << ", Node:";
                for(icnode=0; icnode< numOfConNode; icnode++){
                    pConNode= pSkinFace->getNode(icnode);

                    if(pSkinFace->isSelf()){
                        ofs << " ID= " << pConNode->getNodeID() << ", MeshID= " << pConNode->getMeshID() << ",";
                    }else{
                        ofs << " ID= " << "-"                   << ", MeshID= " << "-"                   << ",";
                    }
                };//icnodeループ
                ofs << ", rank:";
                for(icnode=0; icnode< numOfConNode; icnode++){
                    pConNode= pSkinFace->getNode(icnode);

                    ofs << white << pConNode->getRank();

                };//icnodeループ
                ofs << endl;
            };//ifaceループ
        };//maslave切り替えループ
    };//icontループ


    // ----
    // 形状チェック
    // ----
    ofs << "-- 形状チェック --" << endl;
    for(icont=0; icont < numOfContact; icont++){
        pConMesh= pAssy->getContactMesh(icont);
        // マスター,スレーブ切り替え
        // --
        for(maslave=0; maslave< 2; maslave++){
            if(maslave==0) numOfSkinFace= pConMesh->getNumOfMasterFace();//マスター面数
            if(maslave==1) numOfSkinFace= pConMesh->getNumOfSlaveFace(); //スレーブ面数
            // SkinFaceループ
            // --
            for(iface=0; iface< numOfSkinFace; iface++){
                if(maslave==0) pSkinFace= pConMesh->getMasterFace(iface);//マスター面
                if(maslave==1) pSkinFace= pConMesh->getSlaveFace(iface); //スレーブ面

                numOfConNode= pSkinFace->getNumOfNode();

                // SkinFaceID  ContactNodeID .........
                if(maslave==0) ofs << white << " MasterFaceID=" << pSkinFace->getID();
                if(maslave==1) ofs << white << " SlaveFaceID =" << pSkinFace->getID();

                ofs << ", ConNode:" << endl;
                for(icnode=0; icnode< numOfConNode; icnode++){
                    pConNode= pSkinFace->getNode(icnode);

                    if(pSkinFace->isSelf()){
                        ofs << " conID= " << pConNode->getID()
                            << ", x= " << pConNode->getX() << ", y= " << pConNode->getY() << ", z= " << pConNode->getZ() << endl;
                    }else{
                        ofs << " conID= " << pConNode->getID()
                            << ", x= " << pConNode->getX() << ", y= " << pConNode->getY() << ", z= " << pConNode->getZ() << endl;
                    }
                };//icnodeループ
                ofs << endl;
            };//ifaceループ
        };//maslave切り替えループ
    };//icontループ

    // ----
    // マスター面上のスレーブ点出力
    // ----
    numOfContact= pAssy->getNumOfContactMesh();
    
    for(icont=0; icont < numOfContact; icont++){// ContactMeshループ
        pConMesh= pAssy->getContactMesh(icont);
        numOfSkinFace= pConMesh->getNumOfMasterFace();//マスター面数

        for(iface=0; iface< numOfSkinFace; iface++){
            pSkinFace= pConMesh->getMasterFace(iface);//マスター面

            // MasterFaceID  
            ofs << " MasterFaceID= " << pSkinFace->getID();

            ofs << ", SlaveP::";
            uint numOfSlaveP= pSkinFace->getNumOfSlaveNode();
            uint isnode;
            for(isnode=0; isnode < numOfSlaveP; isnode++){
                pConNode= pSkinFace->getSlaveNode(isnode);

                // SlaveConNode ID
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
    // MPI
    pmw::CHecMPI *pMPI= pmw::CHecMPI::Instance();
    MPI_Status stat;
    MPI_Request request;
    

    // ----
    // スレーブ点を主とした,マスター面のCoef出力
    // ----
    ofs << " -- Coef出力チェック -- " << endl;
    
    for(icont=0; icont < numOfContact; icont++){
        pConMesh= pAssy->getContactMesh(icont);
        numOfSPoint= pConMesh->getNumOfSlavePoint();

        //debug
        cout << "FileWriterContactMesh::Write, numOfSPoint= " << numOfSPoint << endl;
        
        //スレーブ点 :ループ
        for(islave=0; islave< numOfSPoint; islave++){
            
            pSlaveNode = pConMesh->getSlaveConNode(islave);

            //debug
            cout << "FileWriterContactMesh::Write, pSlaveConNode id= " << pSlaveNode->getID()
                 << ", myRank = " << pSlaveNode->getRank() << endl;
            
            ofs << " スレーブ点 ConID= " << pSlaveNode->getID() ;
            
            int my_rank;
            pMPI->Comm_rank(MPI_COMM_WORLD, &my_rank);

            if(pSlaveNode->getRank()== my_rank){//pConMesh->getRank()){
                ofs << ", Slave Node_ID= " << pSlaveNode->getNodeID() << endl;
            }else{
                ofs << ", Slave Node_ID= " << " - " << endl;
            }
            
            //マスター面IDがセットされている.
            //
            if(pSlaveNode->have_MasterFaceID(mgLevel)){
                //スレーブ点の所属するMaster面のIDから,MasterFace*を取得
                masterFaceID= pSlaveNode->getMasterFaceID(mgLevel);
                pMasterFace= pConMesh->getMasterFace_ID(masterFaceID);

                ofs << " マスター面 ID= " << masterFaceID << endl;
                
                int ivert, numOfVert;
                numOfVert= (int)pMasterFace->getNumOfNode();
                ////// Coef出力
                for(ivert=0; ivert< numOfVert; ivert++){
                    ofs <<  " coef= " << pMasterFace->getCoef(pSlaveNode->getID(),ivert);
                };
                ofs << endl;

////                // ここで,通信が良いかも.
////                if(pMasterFace->isSelf()){
////                    for(ivert=0; ivert< numOfVert; ivert++){
////                        pMasterNode= pMasterFace->getNode(ivert);
////                        ofs << ", node_id:" << pMasterNode->getNodeID();
////                    };
////                }
////                ofs << endl;
                


                //MPI
                int qnode_id[4];// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 注意!: MPIのためにスタティックにとってる.
                uint rank_count(0);
                for(ivert=0; ivert< numOfVert; ivert++){
                    pMasterNode= pMasterFace->getNode(ivert);
                    qnode_id[ivert]= pMasterNode->getNodeID();

                    if(pMasterNode->getRank()==pConMesh->getRank()) rank_count++;
                };

                if(pConMesh->getRank() != pConMesh->getTransmitRank()){
                    if(rank_count==numOfVert){
                        int transRank= (int)pConMesh->getTransmitRank();
                        pMPI->Send(qnode_id,numOfVert,MPI_INT,transRank,1,MPI_COMM_WORLD);
                        //mpMPI->Isend(qnode_id,4,MPI_INT,transRank,1,MPI_COMM_WORLD,&request);
                        //mpMPI->Wait(&request, &stat);
                    }
                    if(rank_count==0){
                        int rnode_id[4];
                        int recvRank= pConMesh->getTransmitRank();//相手のrankから受信
                        pMPI->Recv(rnode_id,numOfVert,MPI_INT,recvRank,1,MPI_COMM_WORLD,&stat);
                        //mpMPI->Irecv(rnode_id,4,MPI_INT,recvRank,1,MPI_COMM_WORLD,&request);
                    }
                }
            }
        };//islave ループ
        ofs << " -- Coef出力  end  -- " << endl;
    };//icont ループ (ContactMesh)

}











