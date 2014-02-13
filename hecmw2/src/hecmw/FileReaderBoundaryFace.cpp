/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderBoundaryFace.cpp
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
#include "FileReaderBoundaryFace.h"
using namespace FileIO;
CFileReaderBoundaryFace::CFileReaderBoundaryFace()
{
    ;
}
CFileReaderBoundaryFace::~CFileReaderBoundaryFace()
{
    ;
}
string CFileReaderBoundaryFace::Name()
{
    return  "FileReaderBoundaryFace";
}

bool CFileReaderBoundaryFace::Read(ifstream& ifs, string& sLine)
{
    uiint mesh_id, elem_id, ent_id, bface_id, ndof, dof, bnode_id, node_id, numOfBNode, numOfBFace;
    uiint mgLevel(0);
    uiint bface_shape;
    vuint vBNodeID;
    double val;
    uiint bnd_id, bnd_type;
    string s_bnd_type, s_bface_shape;
    istringstream iss;
    if(TagCheck(sLine, FileBlockName::StartBoundaryFace()) ) {
        sLine= getLine(ifs);
        iss.clear();
        iss.str(sLine);
        iss >> bnd_id >> s_bnd_type >> mesh_id >> numOfBNode >> numOfBFace;

        bnd_type= IntBndType(s_bnd_type);
        uiint ibnode;
        for(ibnode=0; ibnode < numOfBNode; ibnode++) {
            sLine= getLine(ifs);
            iss.clear();
            iss.str(sLine);
            iss >> bnode_id >> node_id;
            mpFactory->GeneBoundaryFaceNode(mgLevel, bnd_id, bnd_type, mesh_id, node_id, bnode_id);

            if(bnd_type==pmw::BoundaryType::Dirichlet) {
                iss >> ndof;
                vdouble vVal;
                vVal.resize(ndof);
                for(uiint idof=0; idof < ndof; idof++) {
                    iss >> dof >> vVal[idof];
                    mpFactory->setValue_BoundaryFaceNode(mesh_id, bnd_id, bnode_id, dof, vVal[idof]);
                }
            }
        };
        mpFactory->resizeFaceAggregate(mgLevel, mesh_id, bnd_id);
        while(!ifs.eof()) {
            sLine = getLine(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryFace())) break;
            iss.clear();
            iss.str(sLine);
            iss >> s_bface_shape >> bface_id >> elem_id >> ent_id;
            if(bnd_type==pmw::BoundaryType::Neumann) {
                iss >> dof;
            }
            bface_shape = IntElemType(s_bface_shape);
            vBNodeID.clear();
            uiint nNumOfLocalNode;
            switch(bface_shape) {
            case(pmw::ElementType::Quad):
                nNumOfLocalNode=4;
                break;
            case(pmw::ElementType::Quad2):
                nNumOfLocalNode=8;
                break;
            case(pmw::ElementType::Triangle):
                nNumOfLocalNode=3;
                break;
            case(pmw::ElementType::Triangle2):
                nNumOfLocalNode=6;
                break;
            }
            vBNodeID.resize(nNumOfLocalNode);
            for(ibnode=0; ibnode < nNumOfLocalNode; ibnode++) {
                iss >> bnode_id;
                vBNodeID[ibnode]= bnode_id;
            };
            if(bnd_type==pmw::BoundaryType::Neumann) {
                iss >> val;
            }
            mpFactory->GeneBoundaryFace(mgLevel, bnd_id, bnd_type, bface_shape,
                                        mesh_id, elem_id, ent_id, vBNodeID, bface_id, dof, val);

            //cout << "FileReaderBoundaryFace::Read   val:" << val << "  dof:" << dof << endl;
        };
        mpFactory->initFaceAggregate(mgLevel, mesh_id, bnd_id);
        return true;
    } else {
        return false;
    }
}
bool CFileReaderBoundaryFace::Read_bin(ifstream& ifs)
{
    CFileReaderBinCheck *pBinCheck= CFileReaderBinCheck::Instance();
    bool bOrder= pBinCheck->isByteOrderSwap();
    bool b32, bCheck;
    string sClassName("FileReaderBoundaryFace");
    if( !Check_IntSize(b32, bCheck, sClassName) ) return false;
    char cHead='B';
    if( !TagCheck_Bin(ifs, bCheck, cHead, FileBlockName::StartBoundaryFace(), FileBlockName::BoundaryFace_Len())) return false;
    uiint mesh_id, elem_id, ent_id, bface_id, ndof, dof, bnode_id, node_id, nNumOfBNode, nNumOfBFace;
    uiint mgLevel(0);
    uiint bface_shape;
    vuint vBNodeID;
    double val;
    uiint bnd_id, bnd_type;
    string s_bnd_type, s_shape_type;
    ifs.read((char*)&bnd_id, sizeof(uiint));
    if(bOrder) pBinCheck->ByteOrderSwap(bnd_id);
    Read_BndType(ifs, s_bnd_type);
    ifs.read((char*)&mesh_id, sizeof(uiint));
    if(bOrder) pBinCheck->ByteOrderSwap(mesh_id);
    ifs.read((char*)&nNumOfBNode, sizeof(uiint));
    if(bOrder) pBinCheck->ByteOrderSwap(nNumOfBNode);
    ifs.read((char*)&nNumOfBFace, sizeof(uiint));
    if(bOrder) pBinCheck->ByteOrderSwap(nNumOfBFace);
    bnd_type= IntBndType(s_bnd_type);
    uiint ibnode;
    for(ibnode=0; ibnode < nNumOfBNode; ibnode++) {
        ifs.read((char*)&bnode_id, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(bnode_id);
        ifs.read((char*)&node_id, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(node_id);
        mpFactory->GeneBoundaryFaceNode(mgLevel, bnd_id, bnd_type, mesh_id, node_id, bnode_id);
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
                mpFactory->setValue_BoundaryFaceNode(mesh_id, bnd_id, bnode_id, dof, vVal[idof]);
            }
        }
    };
    mpFactory->resizeFaceAggregate(mgLevel, mesh_id, bnd_id);
    while(!ifs.eof()) {
        if( CFileReader::Check_End(ifs) ) break;
        short nCase;
        char cHexa[5], cTetra[6], cPrism[6], cQuad[5], cTriangle[9], cBeam[5];
        nCase= Read_ElementType(ifs, cHexa, 5, "Hexa", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cTetra, 6, "Tetra", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cPrism, 6, "Prism", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cQuad, 5, "Quad", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cTriangle, 9, "Triangle", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cBeam, 5, "Beam", s_shape_type);
        ifs.read((char*)&bface_id, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(bface_id);
        ifs.read((char*)&elem_id, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(elem_id);
        ifs.read((char*)&ent_id, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(ent_id);
        if(bnd_type==pmw::BoundaryType::Neumann) {
            ifs.read((char*)&dof, sizeof(uiint));
            if(bOrder) pBinCheck->ByteOrderSwap(dof);
        }
        bface_shape = IntElemType(s_shape_type);
        vBNodeID.clear();
        uiint nNumOfLocalNode;
        switch(bface_shape) {
        case(pmw::ElementType::Quad):
            nNumOfLocalNode=4;
            break;
        case(pmw::ElementType::Quad2):
            nNumOfLocalNode=8;
            break;
        case(pmw::ElementType::Triangle):
            nNumOfLocalNode=3;
            break;
        case(pmw::ElementType::Triangle2):
            nNumOfLocalNode=6;
            break;
        default:
            mpLogger->Info(Utility::LoggerMode::Error, sClassName+": mismatch Element Type");
            break;
        }
        vBNodeID.resize(nNumOfLocalNode);
        for(ibnode=0; ibnode < nNumOfLocalNode; ibnode++) {
            ifs.read((char*)&bnode_id, sizeof(uiint));
            if(bOrder) pBinCheck->ByteOrderSwap(bnode_id);
            vBNodeID[ibnode]= bnode_id;
        };
        if(bnd_type==pmw::BoundaryType::Neumann) {
            ifs.read((char*)&val, sizeof(double));
            if(bOrder) pBinCheck->ByteOrderSwap(val);
        }
        mpFactory->GeneBoundaryFace(mgLevel, bnd_id, bnd_type, bface_shape,
                                    mesh_id, elem_id, ent_id, vBNodeID, bface_id, dof, val);
    };
    mpFactory->initFaceAggregate(mgLevel, mesh_id, bnd_id);
    return true;
}
