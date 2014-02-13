/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderBoundaryEdge.cpp
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
#include "ElementType.h"
#include "MeshFactory.h"
#include "FileReaderBoundaryEdge.h"
using namespace FileIO;
CFileReaderBoundaryEdge::CFileReaderBoundaryEdge()
{
    ;
}
CFileReaderBoundaryEdge::~CFileReaderBoundaryEdge()
{
    ;
}
string CFileReaderBoundaryEdge::Name()
{
    return  "FileReaderBoundaryEdge";
}

bool CFileReaderBoundaryEdge::Read(ifstream& ifs, string& sLine)
{
    uiint bnd_id, bnd_type, mesh_id, numOfBNode, numOfBEdge;
    uiint bnode_id, node_id, ndof, dof, mgLevel(0);
    uiint bedge_id, elem_id, ent_id, shape_type;
    vuint  vBNodeID;
    double val(0.0);
    string s_bnd_type, s_shape_type;
    istringstream iss;
    if( TagCheck(sLine, FileBlockName::StartBoundaryEdge()) ) {
        sLine= getLine(ifs);
        iss.clear();
        iss.str(sLine.c_str());
        iss >> bnd_id >> s_bnd_type >> mesh_id >> numOfBNode >> numOfBEdge;

        bnd_type= IntBndType(s_bnd_type);
        uiint ibnode;
        for(ibnode=0; ibnode < numOfBNode; ibnode++) {
            sLine= getLine(ifs);
            iss.clear();
            iss.str(sLine);
            iss >> bnode_id >> node_id;
            mpFactory->GeneBoundaryEdgeNode(mgLevel, bnd_id, bnd_type, mesh_id, node_id, bnode_id);
            if(bnd_type==pmw::BoundaryType::Dirichlet) {
                iss >> ndof;
                vdouble vVal;
                vVal.resize(ndof);
                for(uiint idof=0; idof < ndof; idof++) {
                    iss >> dof >> vVal[idof];
                    mpFactory->setValue_BoundaryEdgeNode(mesh_id, bnd_id, bnode_id, dof, vVal[idof]);
                }
            }
        };
        mpFactory->resizeEdgeAggregate(mgLevel, mesh_id, bnd_id);
        while(!ifs.eof()) {
            sLine= getLine(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryEdge()) ) break;
            iss.clear();
            iss.str(sLine);
            iss >> s_shape_type >> bedge_id >> elem_id >> ent_id;
            if(bnd_type==pmw::BoundaryType::Neumann) {
                iss >> dof;
            }
            shape_type = IntElemType(s_shape_type);
            vBNodeID.clear();
            switch(shape_type) {
            case(pmw::ElementType::Beam):
            case(pmw::ElementType::Line):
                vBNodeID.resize(2);
                iss >> vBNodeID[0] >> vBNodeID[1];
                break;
            case(pmw::ElementType::Beam2):
            case(pmw::ElementType::Line2):
                vBNodeID.resize(3);
                iss >> vBNodeID[0] >> vBNodeID[1] >> vBNodeID[2];
                break;
            default:
                ;
                break;
            }
            if(bnd_type==pmw::BoundaryType::Neumann) {
                iss >> val;
            }
            mpFactory->GeneBoundaryEdge(mgLevel, bnd_id, bnd_type, shape_type,
                                        mesh_id, elem_id, ent_id, vBNodeID, bedge_id, dof, val);
        };
        mpFactory->initEdgeAggregate(mgLevel, mesh_id, bnd_id);
        return true;
    } else {
        return false;
    }
}
bool CFileReaderBoundaryEdge::Read_bin(ifstream& ifs)
{
    CFileReaderBinCheck *pBinCheck= CFileReaderBinCheck::Instance();
    bool bOrder= pBinCheck->isByteOrderSwap();
    bool b32, bCheck;
    string sClassName("FileReaderBoundaryEdge");
    if( !Check_IntSize(b32, bCheck, sClassName) ) return false;
    char cHead='B';
    if( !TagCheck_Bin(ifs, bCheck, cHead, FileBlockName::StartBoundaryEdge(), FileBlockName::BoundaryEdge_Len())) return false;
    uiint bnd_id, bnd_type, mesh_id, nNumOfBNode, nNumOfBEdge;
    uiint bnode_id, node_id, ndof, dof, mgLevel(0);
    uiint bedge_id, elem_id, ent_id, shape_type;
    vuint  vBNodeID;
    double val(0.0);
    string s_bnd_type, s_shape_type;
    char cH;
    ifs.read((char*)&bnd_id, sizeof(uiint));
    if(bOrder) pBinCheck->ByteOrderSwap(bnd_id);
    ifs.read(&cH, 1);
    if(cH=='D') {
        ifs.seekg(8, ios_base::cur);
        s_bnd_type="Dirichlet";
    }
    if(cH=='N') {
        ifs.seekg(6, ios_base::cur);
        s_bnd_type="Neumann";
    }
    ifs.read((char*)&mesh_id, sizeof(uiint));
    if(bOrder) pBinCheck->ByteOrderSwap(mesh_id);
    ifs.read((char*)&nNumOfBNode, sizeof(uiint));
    if(bOrder) pBinCheck->ByteOrderSwap(nNumOfBNode);
    ifs.read((char*)&nNumOfBEdge, sizeof(uiint));
    if(bOrder) pBinCheck->ByteOrderSwap(nNumOfBEdge);
    bnd_type= IntBndType(s_bnd_type);
    uiint ibnode;
    for(ibnode=0; ibnode < nNumOfBNode; ibnode++) {
        ifs.read((char*)&bnode_id, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(bnode_id);
        ifs.read((char*)&node_id, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(node_id);
        mpFactory->GeneBoundaryEdgeNode(mgLevel, bnd_id, bnd_type, mesh_id, node_id, bnode_id);
        if(bnd_type==pmw::BoundaryType::Dirichlet) {
            ifs.read((char*)&ndof, sizeof(uiint));
            if(bOrder) pBinCheck->ByteOrderSwap(ndof);
            vdouble vVal;
            vVal.resize(ndof);
            for(uiint idof=0; idof < ndof; idof++) {
                ifs.read((char*)&dof, sizeof(uiint));
                if(bOrder) pBinCheck->ByteOrderSwap(dof);
                ifs.read((char*)&vVal[idof], sizeof(double));
                if(bOrder) pBinCheck->ByteOrderSwap(vVal[idof]);
                mpFactory->setValue_BoundaryEdgeNode(mesh_id, bnd_id, bnode_id, dof, vVal[idof]);
            }
        }
    };
    mpFactory->resizeEdgeAggregate(mgLevel, mesh_id, bnd_id);
    while(!ifs.eof()) {
        if( Check_End(ifs) ) break;
        short nCase;
        char cHexa[5], cTetra[6], cPrism[6], cQuad[5], cTriangle[9], cBeam[5];
        nCase= Read_ElementType(ifs, cHexa, 5, "Hexa", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cTetra, 6, "Tetra", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cPrism, 6, "Prism", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cQuad, 5, "Quad", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cTriangle, 9, "Triangle", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cBeam, 5, "Beam", s_shape_type);
        ifs.read((char*)&bedge_id, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(bedge_id);
        ifs.read((char*)&elem_id, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(elem_id);
        ifs.read((char*)&ent_id, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(ent_id);
        if(bnd_type==pmw::BoundaryType::Neumann) {
            ifs.read((char*)&dof, sizeof(uiint));
            if(bOrder) pBinCheck->ByteOrderSwap(dof);
        }
        shape_type = IntElemType(s_shape_type);
        vBNodeID.clear();
        switch(shape_type) {
        case(pmw::ElementType::Beam):
        case(pmw::ElementType::Line):
            vBNodeID.resize(2);
            ifs.read((char*)&vBNodeID[0], sizeof(uiint));
            if(bOrder) pBinCheck->ByteOrderSwap(vBNodeID[0]);
            ifs.read((char*)&vBNodeID[1], sizeof(uiint));
            if(bOrder) pBinCheck->ByteOrderSwap(vBNodeID[1]);
            break;
        case(pmw::ElementType::Beam2):
        case(pmw::ElementType::Line2):
            vBNodeID.resize(3);
            ifs.read((char*)&vBNodeID[0], sizeof(uiint));
            if(bOrder) pBinCheck->ByteOrderSwap(vBNodeID[0]);
            ifs.read((char*)&vBNodeID[1], sizeof(uiint));
            if(bOrder) pBinCheck->ByteOrderSwap(vBNodeID[1]);
            ifs.read((char*)&vBNodeID[2], sizeof(uiint));
            if(bOrder) pBinCheck->ByteOrderSwap(vBNodeID[2]);
            break;
        default:
            Utility::CLogger *pLogger= Utility::CLogger::Instance();
            pLogger->Info(Utility::LoggerMode::Error, sClassName+": mismatch Element Type");
            break;
        }
        if(bnd_type==pmw::BoundaryType::Neumann) {
            ifs.read((char*)&val, sizeof(double));
            if(bOrder) pBinCheck->ByteOrderSwap(val);
        }
        mpFactory->GeneBoundaryEdge(mgLevel, bnd_id, bnd_type, shape_type,
                                    mesh_id, elem_id, ent_id, vBNodeID, bedge_id, dof, val);
    };
    mpFactory->initEdgeAggregate(mgLevel, mesh_id, bnd_id);
    return true;
}
