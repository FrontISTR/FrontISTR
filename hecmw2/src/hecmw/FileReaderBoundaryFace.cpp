
#include "ElementType.h"

//
//  FileReaderBoundaryFace.cpp
//
//
//
//                          2009.05.22
//                          2009.05.22
//                          k.Takeda
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

// method
// --
bool CFileReaderBoundaryFace::Read(ifstream& ifs, string& sLine)
{
    uiint mesh_id, elem_id, ent_id, bface_id, ndof, dof, bnode_id, node_id, numOfBNode, numOfBFace;
    uiint mgLevel(0);// mgLevel=="0" :コースグリッド
    uiint bface_shape;
    vuint vBNodeID;
    double val;
    uiint bnd_id, bnd_type;
    string s_bnd_type, s_bface_shape;
    
    istringstream iss;

    if(TagCheck(sLine, FileBlockName::StartBoundaryFace()) ){

        sLine= getLine(ifs);
        iss.clear();
        iss.str(sLine);

        // 境界ID, 境界種類, MeshID, 境界節点数, 境界面数
        iss >> bnd_id >> s_bnd_type >> mesh_id >> numOfBNode >> numOfBFace;
        
        bnd_type= IntBndType(s_bnd_type);
        
        
        uiint ibnode;
        for(ibnode=0; ibnode < numOfBNode; ibnode++){
            sLine= getLine(ifs);
            iss.clear();
            iss.str(sLine);

            // BNodeID, NodeID, X, Y, Z :=> X,Y,Z は未使用
            //iss >> bnode_id >> node_id >> x >> y >> z;
            iss >> bnode_id >> node_id;

            mpFactory->GeneBoundaryFaceNode(mgLevel, bnd_id, bnd_type, mesh_id, node_id, bnode_id);

            //Dirichletの場合は、自由度 節点境界値
            if(bnd_type==pmw::BoundaryType::Dirichlet){
                iss >> ndof;
                vdouble vVal; vVal.resize(ndof);
                for(uiint idof=0; idof < ndof; idof++){
                    iss >> dof >> vVal[idof];

                    mpFactory->setValue_BoundaryFaceNode(mesh_id, bnd_id, bnode_id, dof, vVal[idof]);
                }
            }
        };


        //初期Aggデータ領域確保
        mpFactory->resizeFaceAggregate(mgLevel, mesh_id, bnd_id);
        
        // Face
        while(!ifs.eof()){
            sLine = getLine(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryFace())) break;

            iss.clear();
            iss.str(sLine);

            // 形状タイプ, BFaceID, ElementID, EntityID
            iss >> s_bface_shape >> bface_id >> elem_id >> ent_id;

            // Neumannの場合、自由度
            if(bnd_type==pmw::BoundaryType::Neumann){
                iss >> dof;
            }
            
            bface_shape = IntElemType(s_bface_shape);
            vBNodeID.clear();
            uiint nNumOfLocalNode;
            // BNodeID....BNodeIDの処理
            //
            switch(bface_shape){
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
            for(ibnode=0; ibnode < nNumOfLocalNode; ibnode++){
                iss >> bnode_id;
                vBNodeID[ibnode]= bnode_id;
            };
            // Neumannの場合、境界値
            if(bnd_type==pmw::BoundaryType::Neumann){
                iss >> val;
            }

            mpFactory->GeneBoundaryFace(mgLevel, bnd_id, bnd_type, bface_shape,
                                            mesh_id, elem_id, ent_id, vBNodeID, bface_id, dof, val);//dof,valはNeumannのみが利用,Dirichletはダミー値
        };

        mpFactory->initFaceAggregate(mgLevel, mesh_id, bnd_id);//BNode,BFaceを全てセットしたあとに呼び出す

        return true;
    }else{
        return false;
    }
}

bool CFileReaderBoundaryFace::Read_bin(ifstream& ifs)
{
    CFileReaderBinCheck *pBinCheck= CFileReaderBinCheck::Instance();
    bool bOrder= pBinCheck->isByteOrderSwap();

    //BinCheckのサイズ指定との整合性
    bool b32, bCheck;
    string sClassName("FileReaderBoundaryFace");

    if( !Check_IntSize(b32, bCheck, sClassName) ) return false;

    char cHead='B';
    if( !TagCheck_Bin(ifs, bCheck, cHead, FileBlockName::StartBoundaryFace(), FileBlockName::BoundaryFace_Len())) return false;

    uiint mesh_id, elem_id, ent_id, bface_id, ndof, dof, bnode_id, node_id, nNumOfBNode, nNumOfBFace;
    uiint mgLevel(0);// mgLevel=="0" :コースグリッド
    uiint bface_shape;
    vuint vBNodeID;
    double val;
    uiint bnd_id, bnd_type;
    string s_bnd_type, s_shape_type;

    // 境界ID, 境界種類, MeshID, 境界節点数, 境界面数
    ifs.read((char*)&bnd_id, sizeof(uiint));      if(bOrder) pBinCheck->ByteOrderSwap(bnd_id);
    Read_BndType(ifs, s_bnd_type);
    ifs.read((char*)&mesh_id, sizeof(uiint));     if(bOrder) pBinCheck->ByteOrderSwap(mesh_id);
    ifs.read((char*)&nNumOfBNode, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(nNumOfBNode);
    ifs.read((char*)&nNumOfBFace, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(nNumOfBFace);

    bnd_type= IntBndType(s_bnd_type);

    uiint ibnode;
    for(ibnode=0; ibnode < nNumOfBNode; ibnode++){
        
        // BNodeID, NodeID
        ifs.read((char*)&bnode_id, sizeof(uiint));  if(bOrder) pBinCheck->ByteOrderSwap(bnode_id);
        ifs.read((char*)&node_id, sizeof(uiint));   if(bOrder) pBinCheck->ByteOrderSwap(node_id);

        mpFactory->GeneBoundaryFaceNode(mgLevel, bnd_id, bnd_type, mesh_id, node_id, bnode_id);

        //Dirichletの場合は、自由度 節点境界値
        if(bnd_type==pmw::BoundaryType::Dirichlet){
            ifs.read((char*)&ndof, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(ndof);
            
            vdouble vVal; vVal.resize(ndof);
            for(uiint idof=0; idof < ndof; idof++){
                ifs.read((char*)&dof, sizeof(uiint));  if(bOrder) pBinCheck->ByteOrderSwap(dof);
                ifs.read((char*)&vVal[idof], sizeof(double));  if(bOrder) pBinCheck->ByteOrderSwap(vVal[idof]);

                mpFactory->setValue_BoundaryFaceNode(mesh_id, bnd_id, bnode_id, dof, vVal[idof]);
            }
        }
    };


    //初期Aggデータ領域確保
    mpFactory->resizeFaceAggregate(mgLevel, mesh_id, bnd_id);

    // Face
    while(!ifs.eof()){
        
        if( CFileReader::Check_End(ifs) ) break;

        // 形状タイプ, BFaceID, ElementID, EntityID
        short nCase;
        char cHexa[5], cTetra[6], cPrism[6], cQuad[5], cTriangle[9], cBeam[5];
        // 形状タイプ, BEdgeID, ElementID, EdgeID,
                     nCase= Read_ElementType(ifs, cHexa, 5, "Hexa", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cTetra, 6, "Tetra", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cPrism, 6, "Prism", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cQuad, 5, "Quad", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cTriangle, 9, "Triangle", s_shape_type);
        if(nCase==0) nCase= Read_ElementType(ifs, cBeam, 5, "Beam", s_shape_type);
        ifs.read((char*)&bface_id, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(bface_id);
        ifs.read((char*)&elem_id, sizeof(uiint));  if(bOrder) pBinCheck->ByteOrderSwap(elem_id);
        ifs.read((char*)&ent_id, sizeof(uiint));   if(bOrder) pBinCheck->ByteOrderSwap(ent_id);

        // Neumannの場合、自由度
        if(bnd_type==pmw::BoundaryType::Neumann){
            ifs.read((char*)&dof, sizeof(uiint));  if(bOrder) pBinCheck->ByteOrderSwap(dof);
        }
        
        bface_shape = IntElemType(s_shape_type);
        
        vBNodeID.clear();
        uiint nNumOfLocalNode;
        // BNodeID....BNodeIDの処理
        switch(bface_shape){
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
        for(ibnode=0; ibnode < nNumOfLocalNode; ibnode++){
            ifs.read((char*)&bnode_id, sizeof(uiint));  if(bOrder) pBinCheck->ByteOrderSwap(bnode_id);
            vBNodeID[ibnode]= bnode_id;
        };
        // Neumannの場合、境界値
        if(bnd_type==pmw::BoundaryType::Neumann){
            ifs.read((char*)&val, sizeof(double));  if(bOrder) pBinCheck->ByteOrderSwap(val);
        }
        //dof,valはNeumannのみが利用,Dirichletはダミー値
        mpFactory->GeneBoundaryFace(mgLevel, bnd_id, bnd_type, bface_shape,
                                        mesh_id, elem_id, ent_id, vBNodeID, bface_id, dof, val);
    };
    mpFactory->initFaceAggregate(mgLevel, mesh_id, bnd_id);//BNode,BFaceを全てセットしたあとに呼び出す

    return true;
}


