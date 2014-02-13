/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterUNS.cpp
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
#include "FileWriterUNS.h"
using namespace FileIO;
CFileWriterUNS::CFileWriterUNS()
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
CFileWriterUNS::~CFileWriterUNS()
{
    ;
}
void CFileWriterUNS::WriteDebug(ofstream& ofs, const uiint& mgLevel)
{
    Utility::CLogger* pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn, "invalid method FileWriterUNS::WriteDebug");
}
void CFileWriterUNS::WriteMesh(ofstream& ofs, const uiint& nNumParam, const uiint& iMesh, const uiint& mgLevel)
{
    pmw::CAssyModel *pAssy= mpGMGModel->getAssyModel(mgLevel);

    pmw::CMesh *pMesh= pAssy->getMesh(iMesh);
    uiint nNumOfNode= pMesh->getNumOfNode();
    pmw::CIndexBucket *pBucket= pMesh->getBucket();
    stringstream ss;
    ss << nNumOfNode;
    uiint nLength= ss.str().length();
    uiint nIWidth= nLength, nDWidth=15;

    ofs << "Nodes " << nNumOfNode << endl;
    for(uiint i=0; i < nNumOfNode; i++) {
        pmw::CNode *pNode = pMesh->getNodeIX(i);
        uiint id = pNode->getID();

        double x = pNode->getX();
        double y = pNode->getY();
        double z = pNode->getZ();
        ofs << setw(nDWidth) << right << scientific  << x << " " << y << " " << z << endl;
    }

    ofs << "Boundary Faces 0" << endl;

    ofs << "Elements" << endl;

    vector<pmw::CNode*> vNode;
    for(uiint i=0; i < pMesh->getNumOfElement(); i++) {
        vNode.clear();
        pmw::CElement *pElem = pMesh->getElementIX(i);
        uiint id = pElem->getID();

        vNode = pElem->getNode();
        if(pElem->getType() == pmw::ElementType::Hexa) {
            ofs << "2 1" << endl;
            uiint iNodeSort[] = {0,1,3,2,4,5,7,6};
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Hexa(); iENode++) {
                uiint iENode_fvuns = iNodeSort[iENode];
                uiint node_id = vNode[iENode_fvuns]->getID();
                uiint index = pBucket->getIndexNode(node_id) + 1;
                ofs << setw(nIWidth) << right << dec << index << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Hexa2) {
            ofs << "2 1" << endl;
            uiint iNodeSort[] = {0,1,3,2,4,5,7,6};
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Hexa2(); iENode++) {
                uiint iENode_fvuns = iNodeSort[iENode];
                uiint node_id = vNode[iENode_fvuns]->getID();
                uiint index = pBucket->getIndexNode(node_id) + 1;
                ofs << setw(nIWidth) << right << dec << index << " ";
                if( iENode == 7 ) break;
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Tetra) {
            ofs << "1 1" << endl;
            uiint iNodeSort[] = {3,2,0,1};
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Tetra(); iENode++) {
                uiint iENode_fvuns = iNodeSort[iENode];
                uiint node_id = vNode[iENode_fvuns]->getID();
                uiint index = pBucket->getIndexNode(node_id) + 1;
                ofs << setw(nIWidth) << right << dec << index << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Tetra2) {
            ofs << "1 1" << endl;
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Tetra2(); iENode++) {
                uiint iNodeSort[] = {3,2,0,1};
                uiint iENode_fvuns = iNodeSort[iENode];
                uiint node_id = vNode[iENode_fvuns]->getID();
                uiint index = pBucket->getIndexNode(node_id) + 1;
                ofs << setw(nIWidth) << right << dec << index << " ";
                if( iENode == 3 ) break;
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Prism) {
            ofs << "3 1" << endl;
            uiint iNodeSort[] = {0,3,4,1,5,2};
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Prism(); iENode++) {
                uiint iENode_fvuns = iNodeSort[iENode];
                uiint node_id = vNode[iENode_fvuns]->getID();
                uiint index = pBucket->getIndexNode(node_id) + 1;
                ofs << setw(nIWidth) << right << dec << index << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Prism2) {
            ofs << "3 1" << endl;
            uiint iNodeSort[] = {0,3,4,1,5,2};
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Prism2(); iENode++) {
                uiint iENode_fvuns = iNodeSort[iENode];
                uiint node_id = vNode[iENode_fvuns]->getID();
                uiint index = pBucket->getIndexNode(node_id) + 1;
                ofs << setw(nIWidth) << right << dec << index << " ";
                if( iENode == 5 ) break;
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Quad) {
            ofs << "2 1" << endl;
            uiint iNodeSort[] = {0,1,3,2};
            for(uiint ii=0; ii < 2; ii++ ) {
                for(uiint iENode=0; iENode < pmw::NumberOfNode::Quad(); iENode++) {
                    uiint iENode_fvuns = iNodeSort[iENode];
                    uiint node_id = vNode[iENode_fvuns]->getID();
                    uiint index = pBucket->getIndexNode(node_id) + 1;
                    ofs << setw(nIWidth) << right << dec << index << " ";
                };
            }
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Quad2) {
            ofs << "2 1" << endl;
            uiint iNodeSort[] = {0,1,3,2};
            for(uiint ii=0; ii < 2; ii++ ) {
                for(uiint iENode=0; iENode < pmw::NumberOfNode::Quad2(); iENode++) {
                    uiint iENode_fvuns = iNodeSort[iENode];
                    uiint node_id = vNode[iENode_fvuns]->getID();
                    uiint index = pBucket->getIndexNode(node_id) + 1;
                    ofs << setw(nIWidth) << right << dec << index << " ";
                    if( iENode == 3 ) break;
                };
            }
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Triangle) {
            ofs << "1 1" << endl;
            uiint iNodeSort[] = {0,1,2};
            for(uiint ii=0; ii < 2; ii++ ) {
                for(uiint iENode=0; iENode < pmw::NumberOfNode::Triangle(); iENode++) {
                    uiint iENode_fvuns = iNodeSort[iENode];
                    uiint node_id = vNode[iENode_fvuns]->getID();
                    uiint index = pBucket->getIndexNode(node_id) + 1;
                    ofs << setw(nIWidth) << right << dec << index << " ";
                };
            }
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Triangle2) {
            ofs << "1 1" << endl;
            uiint iNodeSort[] = {0,1,2};
            for(uiint ii=0; ii < 2; ii++ ) {
                for(uiint iENode=0; iENode < pmw::NumberOfNode::Triangle2(); iENode++) {
                    uiint iENode_fvuns = iNodeSort[iENode];
                    uiint node_id = vNode[iENode_fvuns]->getID();
                    uiint index = pBucket->getIndexNode(node_id) + 1;
                    ofs << setw(nIWidth) << right << dec << index << " ";
                    if( iENode == 2 ) break;
                };
            }
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Beam) {
            ofs << "1 1" << endl;
            uiint iNodeSort[] = {0,1};
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Beam(); iENode++) {
                uiint iENode_fvuns = iNodeSort[iENode];
                uiint node_id = vNode[iENode_fvuns]->getID();
                uiint index = pBucket->getIndexNode(node_id) + 1;
                ofs << setw(nIWidth) << right << dec << index << " ";
                ofs << setw(nIWidth) << right << dec << index << " ";
                ofs << setw(nIWidth) << right << dec << index << " ";
                ofs << setw(nIWidth) << right << dec << index << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Beam2) {
            ofs << "1 1" << endl;
            uiint iNodeSort[] = {0,1};
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Beam2(); iENode++) {
                uiint iENode_fvuns = iNodeSort[iENode];
                uiint node_id = vNode[iENode_fvuns]->getID();
                uiint index = pBucket->getIndexNode(node_id) + 1;
                ofs << setw(nIWidth) << right << dec << index << " ";
                ofs << setw(nIWidth) << right << dec << index << " ";
                ofs << setw(nIWidth) << right << dec << index << " ";
                ofs << setw(nIWidth) << right << dec << index << " ";
                if( iENode == 2 ) break;
            };
            ofs << endl;
        }
    };
}

void CFileWriterUNS::recUNS_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF)
{
    string sLabel=cLabel;
    mvLabel[iMesh].push_back(sLabel);
    string sUnit=cUnit;
    mmUnit[iMesh][sLabel]= sUnit;
    mmDOF[iMesh][sLabel]= nNumOfDOF;

    mvNumOfParameter[iMesh] += nNumOfDOF;
}
void CFileWriterUNS::recUNS_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue)
{
    string sLabel=cLabel;
    uiint nDOF= mmDOF[iMesh][sLabel];
    vdouble vValue;
    vValue.resize(nNumOfNode*nDOF);
    for(uiint i=0; i < vValue.size(); i++) vValue[i]=pvValue[i];
    mmVecVal[iMesh][sLabel]= vValue;
}
void CFileWriterUNS::WriteFEM(ofstream& ofs, const uiint& nNumOfMesh, const uiint& mgLevel)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();

    pmw::CAssyModel *pAssy= mpGMGModel->getAssyModel(mgLevel);

    ofs << "FIELDVIEW 2 5" << endl;
    ofs << "Constants" << endl;
    ofs << "0." << endl;	//Solution time
    ofs << "0." << endl;	//FSMACH
    ofs << "0." << endl;	//ALPHA
    ofs << "0." << endl;	//RE

    ofs << "Grids" << endl;
    ofs << nNumOfMesh << endl;

    ofs << "Boundary Table 0" << endl;

    ofs << "Variable Names " << 3 << endl;
    uiint nNumOfCompo = mvLabel[0].size();
    for(uiint i=0; i < nNumOfCompo; i++) {
        string sLabel= mvLabel[0][i];
        uiint nDOF= mmDOF[0][sLabel];
        if( nDOF == 3 || nDOF == 2 ) {
            ofs << "x_" << sLabel << "; "<< sLabel << endl;
            ofs << "y_" << sLabel << endl;
            ofs << "z_" << sLabel << endl;
        } else {
            ofs << sLabel << endl;
        }
    }

    ofs << "Boundary Variable Names 0" << endl;

    for(uiint iMesh=0; iMesh < nNumOfMesh; iMesh++) {
        CFileWriterUNS::WriteMesh(ofs, mvNumOfParameter[iMesh], iMesh, mgLevel);

        pmw::CMesh *pMesh= pAssy->getMesh(iMesh);
        uiint nNumOfNode= pMesh->getNumOfNode();
        stringstream ss;
        ss << nNumOfNode;
        uiint nLength= ss.str().length();
        uiint nIWidth= nLength, nDWidth=15;
        uiint iNode;

        ofs << "Variables" << endl;
        for(uiint i=0; i < nNumOfCompo; i++) {
            for(uiint icompo=0; icompo < nNumOfCompo; icompo++) {
                string sLabel= mvLabel[iMesh][icompo];
                uiint nDOF= mmDOF[iMesh][sLabel];
                for(uiint idof=0; idof < nDOF; idof++) {
                    for(iNode=0; iNode < nNumOfNode; iNode++) {
                        ofs << setw(nDWidth) << right << scientific << mmVecVal[iMesh][sLabel][iNode*nDOF + idof] << " ";
                        if( iNode % 5 == 4 ) {
                            ofs << endl;
                        }
                    };
                    if( iNode % 5 != 0 ) {
                        ofs << endl;
                    }
                }
            }
        }
    }

    pLogger->Info(Utility::LoggerMode::Info, "Output FVUNS file");
}
