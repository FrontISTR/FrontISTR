/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderBoundaryVolume.cpp
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
#include "FileReaderBoundaryVolume.h"
using namespace FileIO;
CFileReaderBoundaryVolume::CFileReaderBoundaryVolume()
{
    ;
}
CFileReaderBoundaryVolume::~CFileReaderBoundaryVolume()
{
    ;
}
string CFileReaderBoundaryVolume::Name()
{
    return  "FileReaderBoundaryVolume";
}

bool CFileReaderBoundaryVolume::Read(ifstream& ifs, string& sLine)
{
    uiint bnd_id, bnd_type, mesh_id, numOfBNode, numOfBVol;
    uiint bnode_id, node_id, ndof, dof, bvol_id, elem_id, ent_id, shape_type;
    vuint vBNodeID;
    double val;
    string  s_bnd_type, s_shape_type;
    uiint mgLevel(0);
    istringstream iss;
    if(TagCheck(sLine, FileBlockName::StartBoundaryVolume()) ) {
        sLine = getLine(ifs);
        iss.clear();
        iss.str(sLine);
        iss >> bnd_id >> s_bnd_type >> mesh_id >> numOfBNode >> numOfBVol;
        bnd_type = IntBndType(s_bnd_type);
        uiint ibnode;
        for(ibnode=0; ibnode < numOfBNode; ibnode++) {
            sLine= getLine(ifs);
            iss.clear();
            iss.str(sLine);
            iss >> bnode_id >> node_id;

            mpFactory->GeneBoundaryVolumeNode(mgLevel, bnd_id, bnd_type, mesh_id, node_id, bnode_id);

            if(bnd_type==pmw::BoundaryType::Dirichlet) {
                iss >> ndof;
                vdouble vVal;
                vVal.resize(ndof);
                for(uiint idof=0; idof < ndof; idof++) {
                    iss >> dof >> vVal[idof];

                    mpFactory->setValue_BoundaryVolumeNode(mesh_id, bnd_id, bnode_id, dof, vVal[idof]);
                }
            }
        };
        mpFactory->resizeVolumeAggregate(mgLevel, mesh_id, bnd_id);

        while(!ifs.eof()) {
            sLine= getLine(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryVolume())) break;
            iss.clear();
            iss.str(sLine);
            iss >> s_shape_type >> bvol_id >> elem_id >> ent_id;
            if(bnd_type==pmw::BoundaryType::Neumann) {
                iss >> dof;
            }

            ////cout << "FileReaderBoundaryVolume::Read  ShapeType:" << s_shape_type << endl;//debug

            shape_type = IntElemType(s_shape_type);

            vBNodeID.clear();
            switch(shape_type) {
            case(pmw::ElementType::Hexa):
                vBNodeID.resize(pmw::NumberOfNode::Hexa());
                for(ibnode=0; ibnode < pmw::NumberOfNode::Hexa(); ibnode++) {
                    iss >> vBNodeID[ibnode];
                };
                break;
            case(pmw::ElementType::Hexa2):
                vBNodeID.resize(pmw::NumberOfNode::Hexa2());
                for(ibnode=0; ibnode < pmw::NumberOfNode::Hexa2(); ibnode++) {
                    iss >> vBNodeID[ibnode];
                };
                break;
            case(pmw::ElementType::Tetra):
                vBNodeID.resize(pmw::NumberOfNode::Tetra());
                for(ibnode=0; ibnode < pmw::NumberOfNode::Tetra(); ibnode++) {
                    iss >> vBNodeID[ibnode];
                };
                break;
            case(pmw::ElementType::Tetra2):
                vBNodeID.resize(pmw::NumberOfNode::Tetra2());
                for(ibnode=0; ibnode < pmw::NumberOfNode::Tetra2(); ibnode++) {
                    iss >> vBNodeID[ibnode];
                };
                break;
            case(pmw::ElementType::Prism):
                vBNodeID.resize(pmw::NumberOfNode::Prism());
                for(ibnode=0; ibnode < pmw::NumberOfNode::Prism(); ibnode++) {
                    iss >> vBNodeID[ibnode];
                };
                break;
            case(pmw::ElementType::Prism2):
                vBNodeID.resize(pmw::NumberOfNode::Prism2());
                for(ibnode=0; ibnode < pmw::NumberOfNode::Prism2(); ibnode++) {
                    iss >> vBNodeID[ibnode];
                };
                break;
            default:
                break;
            }
            if(bnd_type==pmw::BoundaryType::Neumann) {
                iss >> val;
            }
            mpFactory->GeneBoundaryVolume(mgLevel, bnd_id, bnd_type, shape_type,
                                          mesh_id, elem_id, vBNodeID, bvol_id, dof, val);
        };
        mpFactory->initVolumeAggregate(mgLevel, mesh_id, bnd_id);
        return true;
    } else {
        return false;
    }
}
bool CFileReaderBoundaryVolume::Read_bin(ifstream& ifs)
{
    CFileReaderBinCheck *pBinCheck= CFileReaderBinCheck::Instance();
    bool bOrder= pBinCheck->isByteOrderSwap();
    bool b32, bCheck;
    string sClassName("FileReaderBoundaryVolume");
    if( !Check_IntSize(b32, bCheck, sClassName) ) return false;
    char cHead='B';
    if( !TagCheck_Bin(ifs, bCheck, cHead, FileBlockName::StartBoundaryVolume(), FileBlockName::BoundaryVolume_Len())) return false;
    uiint bnd_id, bnd_type, mesh_id, nNumOfBNode, nNumOfBVol;
    uiint bnode_id, node_id, ndof, dof, bvol_id, elem_id, ent_id, shape_type;
    vuint vBNodeID;
    double val;
    string  s_bnd_type, s_shape_type;
    uiint mgLevel(0);
    ifs.read((char*)&bnd_id, sizeof(uiint));
    if(bOrder) pBinCheck->ByteOrderSwap(bnd_id);
    CFileReader::Read_BndType(ifs, s_bnd_type);
    ifs.read((char*)&mesh_id, sizeof(uiint));
    if(bOrder) pBinCheck->ByteOrderSwap(mesh_id);
    ifs.read((char*)&nNumOfBNode, sizeof(uiint));
    if(bOrder) pBinCheck->ByteOrderSwap(nNumOfBNode);
    ifs.read((char*)&nNumOfBVol, sizeof(uiint));
    if(bOrder) pBinCheck->ByteOrderSwap(nNumOfBVol);
    bnd_type = IntBndType(s_bnd_type);
    uiint ibnode;
    for(ibnode=0; ibnode < nNumOfBNode; ibnode++) {
        ifs.read((char*)&bnode_id, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(bnode_id);
        ifs.read((char*)&node_id, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(node_id);
        mpFactory->GeneBoundaryVolumeNode(mgLevel, bnd_id, bnd_type, mesh_id, node_id, bnode_id);
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
                mpFactory->setValue_BoundaryVolumeNode(mesh_id, bnd_id, bnode_id, dof, vVal[idof]);
            }
        }
    };
    mpFactory->resizeVolumeAggregate(mgLevel, mesh_id, bnd_id);
    while(!ifs.eof()) {
        if( CFileReader::Check_End(ifs) ) break;
        short nCase;
        char cHexa[5], cTetra[6], cPrism[6];
        nCase= Read_ElementType(ifs, cHexa, 5, "Hexa", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cTetra, 6, "Tetra", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cPrism, 6, "Prism", s_shape_type);
        ifs.read((char*)&bvol_id, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(bvol_id);
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
        case(pmw::ElementType::Hexa):
            vBNodeID.resize(pmw::NumberOfNode::Hexa());
            for(ibnode=0; ibnode < pmw::NumberOfNode::Hexa(); ibnode++) {
                ifs.read((char*)&vBNodeID[ibnode], sizeof(uiint));
                if(bOrder) pBinCheck->ByteOrderSwap(vBNodeID[ibnode]);
            };
            break;
        case(pmw::ElementType::Hexa2):
            vBNodeID.resize(pmw::NumberOfNode::Hexa2());
            for(ibnode=0; ibnode < pmw::NumberOfNode::Hexa2(); ibnode++) {
                ifs.read((char*)&vBNodeID[ibnode], sizeof(uiint));
                if(bOrder) pBinCheck->ByteOrderSwap(vBNodeID[ibnode]);
            };
            break;
        case(pmw::ElementType::Tetra):
            vBNodeID.resize(pmw::NumberOfNode::Tetra());
            for(ibnode=0; ibnode < pmw::NumberOfNode::Tetra(); ibnode++) {
                ifs.read((char*)&vBNodeID[ibnode], sizeof(uiint));
                if(bOrder) pBinCheck->ByteOrderSwap(vBNodeID[ibnode]);
            };
            break;
        case(pmw::ElementType::Tetra2):
            vBNodeID.resize(pmw::NumberOfNode::Tetra2());
            for(ibnode=0; ibnode < pmw::NumberOfNode::Tetra2(); ibnode++) {
                ifs.read((char*)&vBNodeID[ibnode], sizeof(uiint));
                if(bOrder) pBinCheck->ByteOrderSwap(vBNodeID[ibnode]);
            };
            break;
        case(pmw::ElementType::Prism):
            vBNodeID.resize(pmw::NumberOfNode::Prism());
            for(ibnode=0; ibnode < pmw::NumberOfNode::Prism(); ibnode++) {
                ifs.read((char*)&vBNodeID[ibnode], sizeof(uiint));
                if(bOrder) pBinCheck->ByteOrderSwap(vBNodeID[ibnode]);
            };
            break;
        case(pmw::ElementType::Prism2):
            vBNodeID.resize(pmw::NumberOfNode::Prism2());
            for(ibnode=0; ibnode < pmw::NumberOfNode::Prism2(); ibnode++) {
                ifs.read((char*)&vBNodeID[ibnode], sizeof(uiint));
                if(bOrder) pBinCheck->ByteOrderSwap(vBNodeID[ibnode]);
            };
            break;
        default:
            mpLogger->Info(Utility::LoggerMode::Error, sClassName+": mismatch Element Type");
            break;
        }
        if(bnd_type==pmw::BoundaryType::Neumann) {
            ifs.read((char*)&val, sizeof(double));
            if(bOrder) pBinCheck->ByteOrderSwap(val);
        }
        mpFactory->GeneBoundaryVolume(mgLevel, bnd_id, bnd_type, shape_type,
                                      mesh_id, elem_id, vBNodeID, bvol_id, dof, val);
    };
    mpFactory->initVolumeAggregate(mgLevel, mesh_id, bnd_id);
    return true;
}
