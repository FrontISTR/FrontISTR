/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterAVS.cpp
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
#include "FileWriterAVS.h"
using namespace FileIO;
CFileWriterAVS::CFileWriterAVS()
{
    pmw::CAssyModel *pAssy= mpGMGModel->getAssyModel(0);
    uiint nNumOfMesh = pAssy->getNumOfMesh();
    mvLabel.resize(nNumOfMesh);
    mmUnit.resize(nNumOfMesh);
    mmDOF.resize(nNumOfMesh);
    mmVecVal.resize(nNumOfMesh);

    mvNumOfParameter.resize(nNumOfMesh);
    for(uiint i=0; i < nNumOfMesh; i++) {
        mvNumOfParameter[i]= 0;
    };
}
CFileWriterAVS::~CFileWriterAVS()
{
    ;
}
void CFileWriterAVS::WriteDebug(ofstream& ofs, const uiint& mgLevel)
{
    Utility::CLogger* pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn, "invalid method FileWriterAVS::WriteDebug");
}
void CFileWriterAVS::WriteMesh(ofstream& ofs, const uiint& nNumParam, const uiint& iMesh, const uiint& mgLevel)
{
    pmw::CAssyModel *pAssy= mpGMGModel->getAssyModel(mgLevel);

    pmw::CMesh *pMesh= pAssy->getMesh(iMesh);
    uiint nNumOfNode= pMesh->getNumOfNode();
    stringstream ss;
    ss << nNumOfNode;
    uiint nLength= ss.str().length();
    uiint nIWidth= nLength, nDWidth=15;
    ofs << setw(nIWidth) << right << dec << pAssy->getMesh(iMesh)->getNumOfNode()    << " ";
    ofs << setw(nIWidth) << right << dec << pAssy->getMesh(iMesh)->getNumOfElement() << " ";
    ofs << setw(nIWidth) << right << dec << nNumParam << " " << 0 << " " << 0 << endl;

    bool bNodeIDCheck(false);
    for(uiint i=0; i < pAssy->getMesh(iMesh)->getNumOfNode(); i++) {
        pmw::CNode *pNode = pAssy->getMesh(iMesh)->getNodeIX(i);
        uiint id = pNode->getID();

        if(i==0 && id==0) bNodeIDCheck= true;
        if(bNodeIDCheck) id += 1;//------------ MicroAVS NodeID > 0

        double x = pAssy->getMesh(iMesh)->getNodeIX(i)->getX();
        double y = pAssy->getMesh(iMesh)->getNodeIX(i)->getY();
        double z = pAssy->getMesh(iMesh)->getNodeIX(i)->getZ();
        ofs << setw(nIWidth) << right << dec << id << " ";
        ofs << setw(nDWidth) << right << scientific  << x << " " << y << " " << z << endl;
    }
    vector<pmw::CNode*> vNode;
    bool bElemIDCheck(false);
    for(uiint i=0; i < pAssy->getMesh(iMesh)->getNumOfElement(); i++) {
        vNode.clear();
        pmw::CElement *pElem = pAssy->getMesh(iMesh)->getElementIX(i);
        uiint id = pElem->getID();

        if(i==0 && id==0) bElemIDCheck= true;
        if(bElemIDCheck) id += 1;//------------ MicroAVS ElementID > 0


        vNode = pElem->getNode();
        if(pElem->getType() == pmw::ElementType::Hexa) {
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " hex  ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Hexa(); iENode++) {
                uiint node_id = vNode[iENode]->getID();
                if(bNodeIDCheck) node_id += 1;//---------- MicroAVS NodeID > 0
                ofs << setw(nIWidth) << right << dec << node_id << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Hexa2) {
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " hex2  ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Hexa2(); iENode++) {
                uiint node_id = vNode[iENode]->getID();
                if(bNodeIDCheck) node_id += 1;//---------- MicroAVS NodeID > 0
                ofs << setw(nIWidth) << right << dec << node_id << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Tetra) {
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " tet ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Tetra(); iENode++) {
                uiint node_id = vNode[iENode]->getID();
                if(bNodeIDCheck) node_id += 1;//---------- MicroAVS NodeID > 0
                ofs << setw(nIWidth) << right << dec << node_id << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Tetra2) {
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " tet2 ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Tetra2(); iENode++) {
                uiint node_id = vNode[iENode]->getID();
                if(bNodeIDCheck) node_id += 1;//---------- MicroAVS NodeID > 0
                ofs << setw(nIWidth) << right << dec << node_id << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Prism) {
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " prism ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Prism(); iENode++) {
                uiint node_id = vNode[iENode]->getID();
                if(bNodeIDCheck) node_id += 1;//---------- MicroAVS NodeID > 0
                ofs << setw(nIWidth) << right << dec << node_id << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Prism2) {
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " prism2 ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Prism2(); iENode++) {
                uiint node_id = vNode[iENode]->getID();
                if(bNodeIDCheck) node_id += 1;//---------- MicroAVS NodeID > 0
                ofs << setw(nIWidth) << right << dec << node_id << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Quad) {
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " quad ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Quad(); iENode++) {
                uiint node_id = vNode[iENode]->getID();
                if(bNodeIDCheck) node_id += 1;//---------- MicroAVS NodeID > 0
                ofs << setw(nIWidth) << right << dec << node_id << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Quad2) {
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " quad2 ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Quad2(); iENode++) {
                uiint node_id = vNode[iENode]->getID();
                if(bNodeIDCheck) node_id += 1;//---------- MicroAVS NodeID > 0
                ofs << setw(nIWidth) << right << dec << node_id << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Triangle) {
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " tri ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Triangle(); iENode++) {
                uiint node_id = vNode[iENode]->getID();
                if(bNodeIDCheck) node_id += 1;//---------- MicroAVS NodeID > 0
                ofs << setw(nIWidth) << right << dec << node_id << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Triangle2) {
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " tri2 ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Triangle2(); iENode++) {
                uiint node_id = vNode[iENode]->getID();
                if(bNodeIDCheck) node_id += 1;//---------- MicroAVS NodeID > 0
                ofs << setw(nIWidth) << right << dec << node_id << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Beam) {
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " line ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Beam(); iENode++) {
                uiint node_id = vNode[iENode]->getID();
                if(bNodeIDCheck) node_id += 1;//---------- MicroAVS NodeID > 0
                ofs << setw(nIWidth) << right << dec << node_id << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Beam2) {
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " line2 ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Beam2(); iENode++) {
                uiint node_id = vNode[iENode]->getID();
                if(bNodeIDCheck) node_id += 1;//---------- MicroAVS NodeID > 0
                ofs << setw(nIWidth) << right << dec << node_id << " ";
            };
            ofs << endl;
        }
    };
}
void CFileWriterAVS::WriteBasis(ofstream& ofs, const uiint& ieq, const uiint& iMesh, const uiint& mgLevel)
{
    pmw::CAssyModel *pAssy= mpGMGModel->getAssyModel(mgLevel);
    pmw::CAssyVector *pSolAssyVector= pAssy->getSolutionAssyVector(ieq);
    uiint nDOF = pSolAssyVector->getDOF(iMesh);

    CFileWriterAVS::WriteMesh(ofs, nDOF, iMesh, mgLevel);

    pmw::CMesh *pMesh= pAssy->getMesh(iMesh);
    uiint nNumOfNode= pMesh->getNumOfNode();
    stringstream ss;
    ss << nNumOfNode;
    uiint nLength= ss.str().length();
    uiint nIWidth= nLength, nDWidth=15;

    ofs << 1 << " " << nDOF << endl;
    ofs << "basis, variable" << endl;

    bool bNodeIDCheck(false);
    for(uiint i=0; i < pAssy->getMesh(iMesh)->getNumOfNode(); i++) {
        pmw::CNode *pNode = pAssy->getMesh(iMesh)->getNodeIX(i);
        uiint id = pNode->getID();

        if(i==0 && id==0) bNodeIDCheck=true;
        if(bNodeIDCheck) id += 1;//-------------------- MicroAVS NodeID > 0

        ofs << setw(nIWidth) << right << dec << id << " ";
        for(uiint idof=0; idof < nDOF; idof++) {
            ofs << setw(nDWidth) << right << scientific << pSolAssyVector->getValue(iMesh, i, idof) << " ";
        }
        ofs << endl;
    };
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Info, "Output UCD file");
}
void CFileWriterAVS::recAVS_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF)
{
    string sLabel=cLabel;
    mvLabel[iMesh].push_back(sLabel);
    string sUnit=cUnit;
    mmUnit[iMesh][sLabel]= sUnit;
    mmDOF[iMesh][sLabel]= nNumOfDOF;

    mvNumOfParameter[iMesh] += nNumOfDOF;
}
void CFileWriterAVS::recAVS_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue)
{
    string sLabel=cLabel;
    uiint nDOF= mmDOF[iMesh][sLabel];
    vdouble vValue;
    vValue.resize(nNumOfNode*nDOF);
    for(uiint i=0; i < vValue.size(); i++) vValue[i]=pvValue[i];
    mmVecVal[iMesh][sLabel]= vValue;
}
void CFileWriterAVS::WriteFEM(ofstream& ofs, const uiint& iMesh, const uiint& mgLevel)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();

    CFileWriterAVS::WriteMesh(ofs, mvNumOfParameter[iMesh], iMesh, mgLevel);

    pmw::CAssyModel *pAssy= mpGMGModel->getAssyModel(mgLevel);
    pmw::CMesh *pMesh= pAssy->getMesh(iMesh);
    uiint nNumOfNode= pMesh->getNumOfNode();
    stringstream ss;
    ss << nNumOfNode;
    uiint nLength= ss.str().length();
    uiint nIWidth= nLength, nDWidth=15;
    uiint nNumOfCompo = mvLabel[iMesh].size();
    ofs << nNumOfCompo << " ";
    for(uiint i=0; i < nNumOfCompo; i++) {
        string sLabel= mvLabel[iMesh][i];
        uiint nDOF= mmDOF[iMesh][sLabel];
        ofs << nDOF  << " ";
    }
    ofs << endl;
    for(uiint i=0; i < nNumOfCompo; i++) {
        string sLabel= mvLabel[iMesh][i];
        ofs << mvLabel[iMesh][i] << ", " << mmUnit[iMesh][sLabel] << endl;
    }
    bool bNodeIDCheck(false);
    for(uiint iNode=0; iNode < nNumOfNode; iNode++) {
        pmw::CNode *pNode = pMesh->getNodeIX(iNode);
        uiint id = pNode->getID();

        if(iNode==0 && id==0) bNodeIDCheck=true;
        if(bNodeIDCheck) id += 1;//-------------------- MicroAVS NodeID > 0

        ofs << setw(nIWidth) << right << dec << id << " ";
        for(uiint icompo=0; icompo < nNumOfCompo; icompo++) {
            string sLabel= mvLabel[iMesh][icompo];
            uiint nDOF= mmDOF[iMesh][sLabel];
            for(uiint idof=0; idof < nDOF; idof++)
                ofs << setw(nDWidth) << right << scientific << mmVecVal[iMesh][sLabel][iNode*nDOF + idof] << " ";
        }
        ofs << endl;
    };
    pLogger->Info(Utility::LoggerMode::Info, "Output UCD file");
}
