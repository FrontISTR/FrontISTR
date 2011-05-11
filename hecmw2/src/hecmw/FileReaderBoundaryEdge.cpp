
#include "ElementType.h"
#include "MeshFactory.h"

//
//  FileReaderBoundaryEdge.cpp
//
//              2010.04.28
//              k.Takeda
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


bool CFileReaderBoundaryEdge::Read(ifstream& ifs, string& sLine)
{
    uiint bnd_id, bnd_type, mesh_id, numOfBNode, numOfBEdge;
    uiint bnode_id, node_id, ndof, dof, mgLevel(0);
    uiint bedge_id, elem_id, ent_id, shape_type;
    vuint  vBNodeID;
    double val(0.0);
    string s_bnd_type, s_shape_type;

    istringstream iss;

    if( TagCheck(sLine, FileBlockName::StartBoundaryEdge()) ){
        
        sLine= getLine(ifs);
        iss.clear();
        iss.str(sLine.c_str());

        // 境界ID, 境界種類, MeshID, 境界節点数, 境界辺数
        iss >> bnd_id >> s_bnd_type >> mesh_id >> numOfBNode >> numOfBEdge;

        bnd_type= IntBndType(s_bnd_type);

        uiint ibnode;
        // 境界節点
        for(ibnode=0; ibnode < numOfBNode; ibnode++){
            sLine= getLine(ifs);
            iss.clear();
            iss.str(sLine);

            // BNodeID, NodeID, X, Y, Z  :=> x,y,z は未使用
            //iss >> bnode_id >> node_id >> x >> y >> z;
            iss >> bnode_id >> node_id;

            mpFactory->GeneBoundaryEdgeNode(mgLevel, bnd_id, bnd_type, mesh_id, node_id, bnode_id);

            //Dirichletの場合は、BNodeの後に 自由度 境界値
            if(bnd_type==pmw::BoundaryType::Dirichlet){
                iss >> ndof;
                vdouble vVal; vVal.resize(ndof);
                for(uiint idof=0; idof < ndof; idof++){
                    iss >> dof >> vVal[idof];

                    mpFactory->setValue_BoundaryEdgeNode(mesh_id, bnd_id, bnode_id, dof, vVal[idof]);
                }
            }
        };

        mpFactory->resizeEdgeAggregate(mgLevel, mesh_id, bnd_id);

        // Edge
        //
        while(!ifs.eof()){
            sLine= getLine(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryEdge()) ) break;

            iss.clear();
            iss.str(sLine);

            // 形状タイプ, BEdgeID, ElementID, EdgeID, 
            iss >> s_shape_type >> bedge_id >> elem_id >> ent_id;

            // Neumannの場合は、自由度
            if(bnd_type==pmw::BoundaryType::Neumann){
                iss >> dof;
            }

            shape_type = IntElemType(s_shape_type);
            vBNodeID.clear();
            //   BNodeID, BNodeID
            switch(shape_type){
                case(pmw::ElementType::Beam):case(pmw::ElementType::Line):
                    vBNodeID.resize(2);
                    iss >> vBNodeID[0] >> vBNodeID[1];// BNodeID, BNodeID
                    break;
                case(pmw::ElementType::Beam2):case(pmw::ElementType::Line2):
                    vBNodeID.resize(3);
                    iss >> vBNodeID[0] >> vBNodeID[1] >> vBNodeID[2];// BNodeID, BNodeID
                    break;
                default:
                    ;//TODO: Logger->Error
                    break;
            }
            // Neumannの場合は、境界値
            if(bnd_type==pmw::BoundaryType::Neumann){
                iss >> val;
            }
            mpFactory->GeneBoundaryEdge(mgLevel, bnd_id, bnd_type, shape_type,
                                        mesh_id, elem_id, ent_id, vBNodeID, bedge_id, dof, val);//dof,valは、Neumannと共有のためDirichletでも使用
        };//while end

        mpFactory->initEdgeAggregate(mgLevel, mesh_id, bnd_id);//BNode,BEdgeを全てセットした後に呼び出す

        return true;
    }else{
        return false;
    }
}

bool CFileReaderBoundaryEdge::Read_bin(ifstream& ifs)
{
    CFileReaderBinCheck *pBinCheck= CFileReaderBinCheck::Instance();
    bool bOrder= pBinCheck->isByteOrderSwap();
    
    //BinCheckのサイズ指定との整合性
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

    // 境界ID, 境界種類, MeshID, 境界節点数, 境界辺数
    ifs.read((char*)&bnd_id, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(bnd_id);
    ifs.read(&cH, 1);
    if(cH=='D'){
        ifs.seekg(8, ios_base::cur);//D_irichlet
        s_bnd_type="Dirichlet";
    }
    if(cH=='N'){
        ifs.seekg(6, ios_base::cur);//N_eumann
        s_bnd_type="Neumann";
    }
    ifs.read((char*)&mesh_id, sizeof(uiint));     if(bOrder) pBinCheck->ByteOrderSwap(mesh_id);
    ifs.read((char*)&nNumOfBNode, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(nNumOfBNode);
    ifs.read((char*)&nNumOfBEdge, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(nNumOfBEdge);

    bnd_type= IntBndType(s_bnd_type);

    uiint ibnode;
    // 境界節点
    for(ibnode=0; ibnode < nNumOfBNode; ibnode++){

        // BNodeID, NodeID
        ifs.read((char*)&bnode_id, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(bnode_id);
        ifs.read((char*)&node_id, sizeof(uiint));  if(bOrder) pBinCheck->ByteOrderSwap(node_id);

        mpFactory->GeneBoundaryEdgeNode(mgLevel, bnd_id, bnd_type, mesh_id, node_id, bnode_id);

        //Dirichletの場合は、BNodeの後に 自由度 境界値
        if(bnd_type==pmw::BoundaryType::Dirichlet){
            //NDOF
            ifs.read((char*)&ndof, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(ndof);

            vdouble vVal; vVal.resize(ndof);
            for(uiint idof=0; idof < ndof; idof++){
                //DOF value
                ifs.read((char*)&dof, sizeof(uiint));         if(bOrder) pBinCheck->ByteOrderSwap(dof);
                ifs.read((char*)&vVal[idof], sizeof(double)); if(bOrder) pBinCheck->ByteOrderSwap(vVal[idof]);

                mpFactory->setValue_BoundaryEdgeNode(mesh_id, bnd_id, bnode_id, dof, vVal[idof]);
            }
        }
    };

    mpFactory->resizeEdgeAggregate(mgLevel, mesh_id, bnd_id);

    // Edge
    //
    while(!ifs.eof()){

        if( Check_End(ifs) ) break;

        short nCase;
        char cHexa[5], cTetra[6], cPrism[6], cQuad[5], cTriangle[9], cBeam[5];
        // 形状タイプ, BEdgeID, ElementID, EdgeID,
                     nCase= Read_ElementType(ifs, cHexa, 5, "Hexa", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cTetra, 6, "Tetra", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cPrism, 6, "Prism", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cQuad, 5, "Quad", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cTriangle, 9, "Triangle", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cBeam, 5, "Beam", s_shape_type);

        ifs.read((char*)&bedge_id, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(bedge_id);
        ifs.read((char*)&elem_id, sizeof(uiint));  if(bOrder) pBinCheck->ByteOrderSwap(elem_id);
        ifs.read((char*)&ent_id, sizeof(uiint));   if(bOrder) pBinCheck->ByteOrderSwap(ent_id);

        // Neumannの場合は、自由度
        if(bnd_type==pmw::BoundaryType::Neumann){
            ifs.read((char*)&dof, sizeof(uiint));  if(bOrder) pBinCheck->ByteOrderSwap(dof);
        }

        shape_type = IntElemType(s_shape_type);

        vBNodeID.clear();
        //   BNodeID, BNodeID
        switch(shape_type){
            case(pmw::ElementType::Beam):case(pmw::ElementType::Line):
                vBNodeID.resize(2);// BNodeID, BNodeID
                ifs.read((char*)&vBNodeID[0], sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(vBNodeID[0]);
                ifs.read((char*)&vBNodeID[1], sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(vBNodeID[1]);
                break;
            case(pmw::ElementType::Beam2):case(pmw::ElementType::Line2):
                vBNodeID.resize(3);// BNodeID, BNodeID, BNodeID
                ifs.read((char*)&vBNodeID[0], sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(vBNodeID[0]);
                ifs.read((char*)&vBNodeID[1], sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(vBNodeID[1]);
                ifs.read((char*)&vBNodeID[2], sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(vBNodeID[2]);
                break;
            default:
                Utility::CLogger *pLogger= Utility::CLogger::Instance();
                pLogger->Info(Utility::LoggerMode::Error, sClassName+": mismatch Element Type");
                break;
        }
        // Neumannの場合は、境界値
        if(bnd_type==pmw::BoundaryType::Neumann){
            ifs.read((char*)&val, sizeof(double));  if(bOrder) pBinCheck->ByteOrderSwap(val);
        }
        mpFactory->GeneBoundaryEdge(mgLevel, bnd_id, bnd_type, shape_type,
                                    mesh_id, elem_id, ent_id, vBNodeID, bedge_id, dof, val);//dof,valは、Neumannと共有のためDirichletでも使用
    };//while end

    mpFactory->initEdgeAggregate(mgLevel, mesh_id, bnd_id);//BNode,BEdgeを全てセットした後に呼び出す

    return true;
}

