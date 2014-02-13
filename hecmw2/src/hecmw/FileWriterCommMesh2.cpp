/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterCommMesh2.cpp
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
#include "FileWriterCommMesh2.h"
using namespace FileIO;
CFileWriterCommMesh2::CFileWriterCommMesh2()
{
    ;
}
CFileWriterCommMesh2::~CFileWriterCommMesh2()
{
    ;
}
void CFileWriterCommMesh2::WriteDebug(ofstream& ofs, const uiint& mgLevel)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pmw::CHecMPI *pMPI = pmw::CHecMPI::Instance();
    uiint nRank = pMPI->getRank();

    pmw::CAssyModel *pAssyModel;
    pmw::CMesh *pMesh;
    pmw::CCommMesh2 *pCommMesh2;
    pAssyModel= mpGMGModel->getAssyModel(mgLevel);
    uiint nNumOfPart;
    nNumOfPart= pAssyModel->getNumOfMesh();

    for(uiint imesh=0; imesh< nNumOfPart; imesh++) {
        pMesh= pAssyModel->getMesh(imesh);
        ofs << " CommMesh2 in Mesh_ID:" << pMesh->getMeshID() << endl;

        uiint nNumOfCommMesh2;
        nNumOfCommMesh2= pMesh->getCommMesh2Size();

        for(uiint icomm=0; icomm< nNumOfCommMesh2; icomm++) {
            pCommMesh2= pMesh->getCommMesh2IX(icomm);
            ofs << "CommMesh2 myRank:" << pCommMesh2->getRank() <<
                ", CommMesh2 transmitRank:" << pCommMesh2->getTrasmitRank() << endl;

            uiint nNumOfCommNode;
            uiint icnode;
            nNumOfCommNode= pCommMesh2->getCommNodeSize();
            pmw::CCommNode *pCommNode;
            pmw::CNode     *pNode;

            //////// debug
            ////vector<void*> vMsgParam;
            ////string sMssg1 = "FileWriterCommMesh2::WriteDebug  icomm";
            ////string sMssg2 = "nNumOfCommNode ";
            ////vMsgParam.push_back((void*)sMssg1.c_str());
            ////vMsgParam.push_back(&icomm);
            ////vMsgParam.push_back((void*)sMssg2.c_str());
            ////vMsgParam.push_back(&nNumOfCommNode);
            ////vMsgParam.push_back(&nRank);
            ////pLogger->Info(Utility::LoggerMode::Info, "%s%u%s%u%u", vMsgParam);


            for(icnode=0; icnode< nNumOfCommNode; icnode++) {
                pCommNode = pCommMesh2->getCommNodeIX(icnode);
                pNode = pCommNode->getNode();
                ofs << " Node ID : " << pNode->getID()
                    << ", X= " << pNode->getX()
                    << ", Y= " << pNode->getY()
                    << ", Z= " << pNode->getZ() << endl;
            };

            uiint nNumOfCommFace;
            uiint iface;
            nNumOfCommFace= pCommMesh2->getCommFaceSize();
            pmw::CCommFace *pCommFace;
            for(iface=0; iface < nNumOfCommFace; iface++) {
                pCommFace= pCommMesh2->getCommFaceIX(iface);

                ofs << "CommFace ID= " << pCommFace->getID();
                uiint nNumOfCommNodeFace= pCommFace->getCommNodeSize();

                uiint nElemID = pCommFace->getElementID();
                ofs << " ElemID:" << nElemID;

                //面形状の場合：要素の局所面番号を出力
                if(pCommFace->getType()==pmw::ElementType::Quad || pCommFace->getType()==pmw::ElementType::Quad2 ||
                   pCommFace->getType()==pmw::ElementType::Triangle || pCommFace->getType()==pmw::ElementType::Triangle2  ) {

                    ofs << " ElemFace:" << pCommFace->getElementFaceID(pMesh->getElement(nElemID));
                }

                for(icnode=0; icnode < nNumOfCommNodeFace; icnode++) {
                    pCommNode= pCommFace->getCommNode(icnode);
                    ofs << ", " << pCommNode->getID();
                };
                ofs << endl;

                for(icnode=0; icnode < nNumOfCommNodeFace; icnode++) {
                    pCommNode= pCommFace->getCommNode(icnode);
                    ofs << ":: X= " << pCommNode->getX() << ", Y= " << pCommNode->getY() << ", Z= " << pCommNode->getZ();
                };
                ofs << endl;
            };
        };
        ofs << " CommMesh2 END " << endl;
    };
}
