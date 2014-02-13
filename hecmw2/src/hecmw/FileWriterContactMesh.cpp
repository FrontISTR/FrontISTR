/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterContactMesh.cpp
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
#include "FileWriterContactMesh.h"
#include "SkinFace.h"
#include "ContactNode.h"
#include "ContactMesh.h"
#include "AssyModel.h"
using namespace FileIO;
CFileWriterContactMesh::CFileWriterContactMesh()
{
    ;
}
CFileWriterContactMesh::~CFileWriterContactMesh()
{
    ;
}
void CFileWriterContactMesh::WriteDebug(ofstream& ofs, const uiint& mgLevel)
{
    pmw::CAssyModel *pAssy;
    pAssy= mpGMGModel->getAssyModel(mgLevel);
    string white(" ");
    pmw::CContactMesh *pConMesh;
    uiint nNumOfContact;
    uiint icont;
    pmw::CSkinFace *pSkinFace;
    uiint nNumOfSkinFace;
    uiint iface;
    pmw::CContactNode *pConNode;
    uiint nNumOfConNode;
    uiint icnode;
    uiint maslave;
    nNumOfContact= pAssy->getNumOfContactMesh();

    for(icont=0; icont < nNumOfContact; icont++) {
        pConMesh= pAssy->getContactMesh(icont);
        //--
        // maslave:0 マスター、maslave:1 スレーブ
        //--
        for(maslave=0; maslave< 2; maslave++) {
            if(maslave==0) nNumOfSkinFace= pConMesh->getNumOfMasterFace();
            if(maslave==1) nNumOfSkinFace= pConMesh->getNumOfSlaveFace();

            ofs << "-- ConNode --" << endl;
            nNumOfConNode= pConMesh->getNumOfConNode();
            for(icnode=0; icnode < nNumOfConNode; icnode++) {
                pConNode= pConMesh->getContactNode(icnode);
                //--
                // オーバーラップ点
                //--
                if(pConNode->isOverlap())
                    ofs << "overlap id:" << pConNode->getID()
                        << " x:" << pConNode->getX() << " y:" << pConNode->getY() << " z:" << pConNode->getZ()
                        << " rank:" << pConNode->getRank() << " overlap_num:" << pConNode->getNumOfOverlapRank() << endl;
            };

            ofs << "-- Face --" << endl;
            for(iface=0; iface< nNumOfSkinFace; iface++) {
                if(maslave==0) pSkinFace= pConMesh->getMasterFace(iface);
                if(maslave==1) pSkinFace= pConMesh->getSlaveFace(iface);

                nNumOfConNode= pSkinFace->getNumOfNode();

                if(maslave==0) ofs << white << " MasterFaceID=" << pSkinFace->getID();
                if(maslave==1) ofs << white << " SlaveFaceID =" << pSkinFace->getID();

                ofs << ", ConNode:";
                for(icnode=0; icnode< nNumOfConNode; icnode++) {
                    pConNode= pSkinFace->getNode(icnode);
                    ofs << white << pConNode->getID();
                };

                ofs << ", Node:";
                for(icnode=0; icnode< nNumOfConNode; icnode++) {
                    pConNode= pSkinFace->getNode(icnode);
                    //--
                    // 面は自身の領域に存在するか.
                    //--
                    if(pSkinFace->isSelf()) {
                        ofs << " ID= " << pConNode->getNodeID() << ", MeshID= " << pConNode->getMeshID() << ",";
                    } else {
                        ofs << " ID= " << "-"                   << ", MeshID= " << "-"                   << ",";
                    }
                };

                ofs << ", rank:";
                for(icnode=0; icnode< nNumOfConNode; icnode++) {
                    pConNode= pSkinFace->getNode(icnode);
                    ofs << white << pConNode->getRank();
                };
                ofs << endl;
            };
        };
    };

    ofs << "-- 形状チェック --" << endl;
    for(icont=0; icont < nNumOfContact; icont++) {
        pConMesh= pAssy->getContactMesh(icont);
        //--
        // maslave:0 マスター、maslave:1 スレーブ
        //--
        for(maslave=0; maslave< 2; maslave++) {
            if(maslave==0) nNumOfSkinFace= pConMesh->getNumOfMasterFace();
            if(maslave==1) nNumOfSkinFace= pConMesh->getNumOfSlaveFace();

            for(iface=0; iface< nNumOfSkinFace; iface++) {
                if(maslave==0) pSkinFace= pConMesh->getMasterFace(iface);
                if(maslave==1) pSkinFace= pConMesh->getSlaveFace(iface);
                nNumOfConNode= pSkinFace->getNumOfNode();
                if(maslave==0) ofs << white << " MasterFaceID=" << pSkinFace->getID();
                if(maslave==1) ofs << white << " SlaveFaceID =" << pSkinFace->getID();
                ofs << ", ConNode:" << endl;

                for(icnode=0; icnode< nNumOfConNode; icnode++) {
                    pConNode= pSkinFace->getNode(icnode);
                    if(pSkinFace->isSelf()) {
                        ofs << " conID= " << pConNode->getID()
                            << ", x= " << pConNode->getX() << ", y= " << pConNode->getY() << ", z= " << pConNode->getZ() << endl;
                    } else {
                        ofs << " conID= " << pConNode->getID()
                            << ", x= " << pConNode->getX() << ", y= " << pConNode->getY() << ", z= " << pConNode->getZ() << endl;
                    }
                };
                ofs << endl;
            };
        };
    };

    for(icont=0; icont < nNumOfContact; icont++) {
        pConMesh= pAssy->getContactMesh(icont);
        nNumOfSkinFace= pConMesh->getNumOfMasterFace();

        for(iface=0; iface< nNumOfSkinFace; iface++) {
            pSkinFace= pConMesh->getMasterFace(iface);
            ofs << " MasterFaceID= " << pSkinFace->getID();
            ofs << ", SlaveP::";
            uiint numOfSlaveP= pSkinFace->getNumOfSlaveNode();
            uiint isnode;

            for(isnode=0; isnode < numOfSlaveP; isnode++) {
                pConNode= pSkinFace->getSlaveNode(isnode);
                ofs << " ConNodeID= " << pConNode->getID() << ",";
            };
            ofs << endl;
        };
    };

    uiint nNumOfSPoint;
    uiint islave;
    uiint masterFaceID;
    pmw::CContactNode* pSlaveNode;
    pmw::CSkinFace* pMasterFace;
    pmw::CContactNode* pMasterNode;
    pmw::CHecMPI *pMPI= pmw::CHecMPI::Instance();
    MPI_Status stat;
    //MPI_Request request;
    MPI_Datatype MPI_UIINT= pMPI->MPI_UIINT();

    ofs << " -- Coef出力チェック -- " << endl;
    for(icont=0; icont < nNumOfContact; icont++) {

        pConMesh= pAssy->getContactMesh(icont);
        nNumOfSPoint= pConMesh->getNumOfSlavePoint();
        //--
        // スレーブ点 ループ
        //--
        for(islave=0; islave< nNumOfSPoint; islave++) {
            pSlaveNode = pConMesh->getSlaveConNode(islave);
            ofs << " スレーブ点 ConID= " << pSlaveNode->getID() ;
            int my_rank;
            pMPI->Comm_rank(MPI_COMM_WORLD, &my_rank);
            //--
            // スレーブ点が自身Rankか.
            //--
            if(pSlaveNode->getRank()== my_rank) {
                ofs << ", Slave Node_ID= " << pSlaveNode->getNodeID() << endl;
            } else {
                ofs << ", Slave Node_ID= " << " - " << endl;
            }
            //--
            // マスター面に載っているスレーブ点か.
            //--
            if(pSlaveNode->have_MasterFaceID(mgLevel)) {

                masterFaceID= pSlaveNode->getMasterFaceID(mgLevel);
                pMasterFace= pConMesh->getMasterFace_ID(masterFaceID);
                ofs << " マスター面 ID= " << masterFaceID << endl;

                uiint ivert, nNumOfVert;
                nNumOfVert= pMasterFace->getNumOfNode();
                for(ivert=0; ivert< nNumOfVert; ivert++) {
                    ofs <<  " coef= " << pMasterFace->getCoef(pSlaveNode->getID(),ivert);
                };
                ofs << endl;
            }
        };
        ofs << " -- Coef出力  end  -- " << endl;
    };
}
