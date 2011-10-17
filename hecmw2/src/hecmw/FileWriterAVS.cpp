/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileWriterAVS.cpp
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
void CFileWriterAVS::WriteMesh(ofstream& ofs, const uiint& iMesh, const uiint& mgLevel)
{
    pmw::CAssyModel *pAssy= mpGMGModel->getAssyModel(mgLevel);
    pmw::CAssyVector *pSolAssyVector= pAssy->getSolutionAssyVector(iMesh);
    pmw::CMesh *pMesh= pAssy->getMesh(iMesh);
    uiint nNumOfNode= pMesh->getNumOfNode();
    stringstream ss;
    ss << nNumOfNode;
    uiint nLength= ss.str().length();
    uiint nIWidth= nLength, nDWidth=15;
    ofs << setw(nIWidth) << right << dec << pAssy->getMesh(iMesh)->getNumOfNode()    << " ";
    ofs << setw(nIWidth) << right << dec << pAssy->getMesh(iMesh)->getNumOfElement() << " ";
    ofs << setw(nIWidth) << right << dec << pSolAssyVector->getDOF() << " " << 0 << " " << 0 << endl;
    for(uiint i=0; i < pAssy->getMesh(iMesh)->getNumOfNode(); i++)
    {
        pmw::CNode *pNode = pAssy->getMesh(iMesh)->getNodeIX(i);
        uiint id = pNode->getID();
        double x = pAssy->getMesh(iMesh)->getNodeIX(i)->getX();
        double y = pAssy->getMesh(iMesh)->getNodeIX(i)->getY();
        double z = pAssy->getMesh(iMesh)->getNodeIX(i)->getZ();
        ofs << setw(nIWidth) << right << dec << id << " ";
        ofs << setw(nDWidth) << right << scientific  << x << " " << y << " " << z << endl;
    }
    vector<pmw::CNode*> vNode;
    for(uiint i=0; i < pAssy->getMesh(iMesh)->getNumOfElement(); i++)
    {
        vNode.clear();
        pmw::CElement *pElem = pAssy->getMesh(iMesh)->getElementIX(i);
        uiint id = pElem->getID();
        vNode = pElem->getNode();
        if(pElem->getType() == pmw::ElementType::Hexa){
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " hex  ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Hexa(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Hexa2){
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " hex2  ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Hexa2(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Tetra){
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " tet ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Tetra(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Tetra2){
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " tet2 ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Tetra2(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Prism){
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " prism ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Prism(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Prism2){
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " prism2 ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Prism2(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Quad){
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " quad ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Quad(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Quad2){
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " quad2 ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Quad2(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Triangle){
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " tri ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Triangle(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Triangle2){
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " tri2 ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Triangle2(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Beam){
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " line ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Beam(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Beam2){
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " line2 ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Beam2(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
    };
}
void CFileWriterAVS::WriteBasis(ofstream& ofs, const uiint& iMesh, const uiint& mgLevel)
{
    pmw::CAssyModel *pAssy= mpGMGModel->getAssyModel(mgLevel);
    pmw::CAssyVector *pSolAssyVector= pAssy->getSolutionAssyVector(iMesh);
    CFileWriterAVS::WriteMesh(ofs, iMesh, mgLevel);
    pmw::CMesh *pMesh= pAssy->getMesh(iMesh);
    uiint nNumOfNode= pMesh->getNumOfNode();
    stringstream ss;
    ss << nNumOfNode;
    uiint nLength= ss.str().length();
    uiint nIWidth= nLength, nDWidth=15;
    uiint nDOF = pSolAssyVector->getDOF();
    ofs << 1 << " " << nDOF << endl;
    ofs << "basis, variable" << endl;
    for(uiint i=0; i < pAssy->getMesh(iMesh)->getNumOfNode(); i++)
    {
        pmw::CNode *pNode = pAssy->getMesh(iMesh)->getNodeIX(i);
        uiint id = pNode->getID();
        ofs << setw(nIWidth) << right << dec << id << " ";
        for(uiint idof=0; idof < nDOF; idof++){
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
    CFileWriterAVS::WriteMesh(ofs, iMesh, mgLevel);
    pmw::CAssyModel *pAssy= mpGMGModel->getAssyModel(mgLevel);
    pmw::CMesh *pMesh= pAssy->getMesh(iMesh);
    uiint nNumOfNode= pMesh->getNumOfNode();
    stringstream ss;
    ss << nNumOfNode;
    uiint nLength= ss.str().length();
    uiint nIWidth= nLength, nDWidth=15;
    uiint nNumOfCompo = mvLabel[iMesh].size();
    ofs << nNumOfCompo << " ";
    for(uiint i=0; i < nNumOfCompo; i++){
        string sLabel= mvLabel[iMesh][i];
        uiint nDOF= mmDOF[iMesh][sLabel];
        ofs << nDOF  << " ";
    }
    ofs << endl;
    for(uiint i=0; i < nNumOfCompo; i++){
        string sLabel= mvLabel[iMesh][i];
        ofs << mvLabel[iMesh][i] << ", " << mmUnit[iMesh][sLabel] << endl;
    }
    for(uiint iNode=0; iNode < nNumOfNode; iNode++)
    {
        pmw::CNode *pNode = pMesh->getNodeIX(iNode);
        uiint id = pNode->getID();
        ofs << setw(nIWidth) << right << dec << id << " ";
        for(uiint icompo=0; icompo < nNumOfCompo; icompo++){
            string sLabel= mvLabel[iMesh][icompo];
            uiint nDOF= mmDOF[iMesh][sLabel];
            for(uiint idof=0; idof < nDOF; idof++)
            ofs << setw(nDWidth) << right << scientific << mmVecVal[iMesh][sLabel][iNode*nDOF + idof] << " ";
        }
        ofs << endl;
    };
    pLogger->Info(Utility::LoggerMode::Info, "Output UCD file");
}
