/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterVTK.cpp
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
#include "FileWriterVTK.h"
using namespace FileIO;
CFileWriterVTK::CFileWriterVTK()
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
CFileWriterVTK::~CFileWriterVTK()
{
    ;
}
void CFileWriterVTK::WriteDebug(ofstream& ofs, const uiint& mgLevel)
{
    Utility::CLogger* pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn, "invalid method FileWriterVTK::WriteDebug");
}
void CFileWriterVTK::WriteMesh(ofstream& ofs, const uiint& nNumParam, const uiint& iMesh, const uiint& mgLevel)
{
    pmw::CAssyModel *pAssy= mpGMGModel->getAssyModel(mgLevel);

    pmw::CMesh *pMesh= pAssy->getMesh(iMesh);
    uiint nNumOfNode= pMesh->getNumOfNode();
    pmw::CIndexBucket *pBucket= pMesh->getBucket();
    stringstream ss;
    ss << nNumOfNode;
    uiint nLength= ss.str().length();
    uiint nIWidth= nLength, nDWidth=15;

    ofs << "      <Points>" << endl;
    ofs << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    for(uiint i=0; i < nNumOfNode; i++) {
        pmw::CNode *pNode = pMesh->getNodeIX(i);
        uiint id = pNode->getID();

        double x = pNode->getX();
        double y = pNode->getY();
        double z = pNode->getZ();
        ofs << setw(nDWidth) << right << scientific  << x << " " << y << " " << z << endl;
    }
    ofs << "        </DataArray>" << endl;
    ofs << "      </Points>" << endl;

    ofs << "      <Cells>" << endl;
    ofs << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
    vector<pmw::CNode*> vNode;
    for(uiint i=0; i < pMesh->getNumOfElement(); i++) {
        vNode.clear();
        pmw::CElement *pElem = pMesh->getElementIX(i);
        uiint id = pElem->getID();

        vNode = pElem->getNode();
        if(pElem->getType() == pmw::ElementType::Hexa) {
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Hexa(); iENode++) {
                uiint node_id = vNode[iENode]->getID();
                uiint index = pBucket->getIndexNode(node_id);
                ofs << setw(nIWidth) << right << dec << index << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Hexa2) {
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Hexa2(); iENode++) {
                uiint node_id = vNode[iENode]->getID();
                uiint index = pBucket->getIndexNode(node_id);
                ofs << setw(nIWidth) << right << dec << index << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Tetra) {
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Tetra(); iENode++) {
                uiint node_id = vNode[iENode]->getID();
                uiint index = pBucket->getIndexNode(node_id);
                ofs << setw(nIWidth) << right << dec << index << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Tetra2) {
            uiint conv_table[10] = {0,1,2,3,5,6,4,7,8,9};// '13.02.18
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Tetra2(); iENode++) {
                uiint iENode2 = conv_table[iENode];
                uiint node_id = vNode[iENode2]->getID();
                uiint index = pBucket->getIndexNode(node_id);
                ofs << setw(nIWidth) << right << dec << index << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Prism) {
            uiint conv_table[6] = {0,2,1,3,5,4};
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Prism(); iENode++) {
                uiint iENode2 = conv_table[iENode];
                uiint node_id = vNode[iENode2]->getID();
                uiint index = pBucket->getIndexNode(node_id);
                ofs << setw(nIWidth) << right << dec << index << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Prism2) { // UNSUPPORTED BY VTK
            // Treat as a  Prism element
            uiint conv_table[6] = {0,2,1,3,5,4};
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Prism(); iENode++) {
                uiint iENode2 = conv_table[iENode];
                uiint node_id = vNode[iENode2]->getID();
                uiint index = pBucket->getIndexNode(node_id);
                ofs << setw(nIWidth) << right << dec << index << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Quad) {
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Quad(); iENode++) {
                uiint node_id = vNode[iENode]->getID();
                uiint index = pBucket->getIndexNode(node_id);
                ofs << setw(nIWidth) << right << dec << index << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Quad2) {
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Quad2(); iENode++) {
                uiint node_id = vNode[iENode]->getID();
                uiint index = pBucket->getIndexNode(node_id);
                ofs << setw(nIWidth) << right << dec << index << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Triangle) {
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Triangle(); iENode++) {
                uiint node_id = vNode[iENode]->getID();
                uiint index = pBucket->getIndexNode(node_id);
                ofs << setw(nIWidth) << right << dec << index << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Triangle2) {
            uiint conv_table[6] = {0,1,2,5,3,4};
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Triangle2(); iENode++) {
                uiint iENode2 = conv_table[iENode];
                uiint node_id = vNode[iENode2]->getID();
                uiint index = pBucket->getIndexNode(node_id);
                ofs << setw(nIWidth) << right << dec << index << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Beam) {
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Beam(); iENode++) {
                uiint node_id = vNode[iENode]->getID();
                uiint index = pBucket->getIndexNode(node_id);
                ofs << setw(nIWidth) << right << dec << index << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Beam2) {
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Beam2(); iENode++) {
                uiint node_id = vNode[iENode]->getID();
                uiint index = pBucket->getIndexNode(node_id);
                ofs << setw(nIWidth) << right << dec << index << " ";
            };
            ofs << endl;
        }
    };
    ofs << "        </DataArray>" << endl;

    uiint offset=0;
    ofs << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
    for(uiint i=0; i < pMesh->getNumOfElement(); i++) {
        vNode.clear();
        pmw::CElement *pElem = pMesh->getElementIX(i);

        vNode = pElem->getNode();
        if(pElem->getType() == pmw::ElementType::Hexa) {
            offset += pmw::NumberOfNode::Hexa();
            ofs << offset << endl;
        }
        if(pElem->getType() == pmw::ElementType::Hexa2) {
            offset += pmw::NumberOfNode::Hexa2();
            ofs << offset << endl;
        }
        if(pElem->getType() == pmw::ElementType::Tetra) {
            offset += pmw::NumberOfNode::Tetra();
            ofs << offset << endl;
        }
        if(pElem->getType() == pmw::ElementType::Tetra2) {
            offset += pmw::NumberOfNode::Tetra2();
            ofs << offset << endl;
        }
        if(pElem->getType() == pmw::ElementType::Prism) {
            offset += pmw::NumberOfNode::Prism();
            ofs << offset << endl;
        }
        if(pElem->getType() == pmw::ElementType::Prism2) { // UNSUPPORTED BY VTK
            // Treat as a  Prism element
            offset += pmw::NumberOfNode::Prism();
            ofs << offset << endl;
        }
        if(pElem->getType() == pmw::ElementType::Quad) {
            offset += pmw::NumberOfNode::Quad();
            ofs << offset << endl;
        }
        if(pElem->getType() == pmw::ElementType::Quad2) {
            offset += pmw::NumberOfNode::Quad2();
            ofs << offset << endl;
        }
        if(pElem->getType() == pmw::ElementType::Triangle) {
            offset += pmw::NumberOfNode::Triangle();
            ofs << offset << endl;
        }
        if(pElem->getType() == pmw::ElementType::Triangle2) {
            offset += pmw::NumberOfNode::Triangle2();
            ofs << offset << endl;
        }
        if(pElem->getType() == pmw::ElementType::Beam) {
            offset += pmw::NumberOfNode::Beam();
            ofs << offset << endl;
        }
        if(pElem->getType() == pmw::ElementType::Beam2) {
            offset += pmw::NumberOfNode::Beam2();
            ofs << offset << endl;
        }
    };
    ofs << "        </DataArray>" << endl;

    ofs << "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">" << endl;
    for(uiint i=0; i < pMesh->getNumOfElement(); i++) {
        vNode.clear();
        pmw::CElement *pElem = pMesh->getElementIX(i);

        vNode = pElem->getNode();
        if(pElem->getType() == pmw::ElementType::Hexa) {
            ofs << 12 << endl;
        }
        if(pElem->getType() == pmw::ElementType::Hexa2) {
            ofs << 25 << endl;
        }
        if(pElem->getType() == pmw::ElementType::Tetra) {
            ofs << 10 << endl;
        }
        if(pElem->getType() == pmw::ElementType::Tetra2) {
            ofs << 24 << endl;
        }
        if(pElem->getType() == pmw::ElementType::Prism) {
            ofs << 13 << endl;
        }
        if(pElem->getType() == pmw::ElementType::Prism2) { // UNSUPPORTED BY VTK
            // Treat as a  Prism element
            ofs << 13 << endl;
        }
        if(pElem->getType() == pmw::ElementType::Quad) {
            ofs << 9 << endl;
        }
        if(pElem->getType() == pmw::ElementType::Quad2) {
            ofs << 23 << endl;
        }
        if(pElem->getType() == pmw::ElementType::Triangle) {
            ofs << 5 << endl;
        }
        if(pElem->getType() == pmw::ElementType::Triangle2) {
            ofs << 22 << endl;
        }
        if(pElem->getType() == pmw::ElementType::Beam) {
            ofs << 3 << endl;
        }
        if(pElem->getType() == pmw::ElementType::Beam2) {
            ofs << 21 << endl;
        }
    };
    ofs << "        </DataArray>" << endl;
    ofs << "      </Cells>" << endl;
}
void CFileWriterVTK::recVTK_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF)
{
    string sLabel=cLabel;
    mvLabel[iMesh].push_back(sLabel);
    string sUnit=cUnit;
    mmUnit[iMesh][sLabel]= sUnit;
    mmDOF[iMesh][sLabel]= nNumOfDOF;

    mvNumOfParameter[iMesh] += nNumOfDOF;
}
void CFileWriterVTK::recVTK_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue)
{
    string sLabel=cLabel;
    uiint nDOF= mmDOF[iMesh][sLabel];
    vdouble vValue;
    vValue.resize(nNumOfNode*nDOF);
    for(uiint i=0; i < vValue.size(); i++) vValue[i]=pvValue[i];
    mmVecVal[iMesh][sLabel]= vValue;
}
void CFileWriterVTK::WriteFEM_Rank0(ofstream& ofs, const uiint& iMesh, const uiint& mgLevel, string& basename, const uiint& nNumOfProcs)
{
    ofs << "<?xml version=\"1.0\"?>" << endl;
    ofs << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
    ofs << "  <PUnstructuredGrid GhostLevel=\"1\">" << endl;
    ofs << "    <PPointData>" << endl;
    uiint nNumOfCompo = mvLabel[iMesh].size();
    for(uiint i=0; i < nNumOfCompo; i++) {
        string sLabel= mvLabel[iMesh][i];
        uiint nDOF= mmDOF[iMesh][sLabel];
        ofs << "        <PDataArray type=\"Float32\" Name=\"" << sLabel << "\"";
        if(nDOF > 1) {
            ofs << " NumberOfComponents=\"" << nDOF << "\"";
        }
        ofs << "/>" << endl;
    }
    ofs << "    </PPointData>" << endl;
    ofs << "    <PPoints>" << endl;
    ofs << "      <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>" << endl;
    ofs << "    </PPoints>" << endl;
    for(uiint iRank=0; iRank < nNumOfProcs; iRank++) {
        ofs << "    <Piece Source=\"" << basename << "." << iRank << "." << iMesh << ".vtu\"/>" << endl;
    }
    ofs << "  </PUnstructuredGrid>" << endl;
    ofs << "</VTKFile>" << endl;
}
void CFileWriterVTK::WriteFEM(ofstream& ofs, const uiint& iMesh, const uiint& mgLevel)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();

    pmw::CAssyModel *pAssy= mpGMGModel->getAssyModel(mgLevel);

    ofs << "<?xml version=\"1.0\"?>" << endl;
    ofs << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
    ofs << "  <UnstructuredGrid>" << endl;

    pmw::CMesh *pMesh= pAssy->getMesh(iMesh);
    uiint nNumOfNode= pMesh->getNumOfNode();
    stringstream ss;
    ss << nNumOfNode;
    uiint nLength= ss.str().length();
    uiint nIWidth= nLength, nDWidth=15;
    ofs << "    <Piece NumberOfPoints=\"" << nNumOfNode
        << "\" NumberOfCells=\"" << pMesh->getNumOfElement() << "\">" << endl;

    ofs << "      <PointData>" << endl;
    uiint nNumOfCompo = mvLabel[iMesh].size();
    for(uiint i=0; i < nNumOfCompo; i++) {
        string sLabel= mvLabel[iMesh][i];
        uiint nDOF= mmDOF[iMesh][sLabel];
        ofs << "        <DataArray type=\"Float32\" Name=\"" << sLabel << "\" ";
        if(nDOF > 1) {
            ofs << "NumberOfComponents=\"" << nDOF << "\" ";
        }
        ofs << "format=\"ascii\">" << endl;

        for(uiint iNode=0; iNode < nNumOfNode; iNode++) {
////            pmw::CNode *pNode = pMesh->getNodeIX(iNode);
            for(uiint idof=0; idof < nDOF; idof++)
                ofs << setw(nDWidth) << right << scientific << mmVecVal[iMesh][sLabel][iNode*nDOF + idof] << " ";
            ofs << endl;
        };

        ofs << "        </DataArray>" << endl;
    }
    ofs << "      </PointData>" << endl;

    CFileWriterVTK::WriteMesh(ofs, mvNumOfParameter[iMesh], iMesh, mgLevel);

    ofs << "    </Piece>" << endl;
    ofs << "  </UnstructuredGrid>" << endl;
    ofs << "</VTKFile>" << endl;

    pLogger->Info(Utility::LoggerMode::Info, "Output VTK file");
}
