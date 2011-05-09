
#include "ContactMesh.h"


#include "AssyModel.h"

//
//  FileWriterContactMesh.cpp
//
//
//
//                  2009.11.20
//                  2009.11.20
//                  k.Takeda
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

                    if(pSkinFace->isSelf()){
                        ofs << white << pConNode->getID();
                    }else{
                        ofs << white << pConNode->getID();
                    }
                };//icnodeループ
                ofs << ", Node:";
                for(icnode=0; icnode< numOfConNode; icnode++){
                    pConNode= pSkinFace->getNode(icnode);

                    if(pSkinFace->isSelf()){
                        ofs << white << pConNode->getNodeID();
                    }else{
                        ofs << white << "-";
                    }
                };//icnodeループ
                ofs << endl;

            };//ifaceループ
        };//maslave切り替えループ
    };//icontループ
}











