
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
    uint bnd_id, bnd_type, mesh_id, numOfBNode, numOfBEdge;
    uint bnode_id, node_id, dof, mgLevel(0);
    uint bedge_id, elem_id, ent_id, shape_type;
    vuint  vBNodeID;
    double val, x, y, z;
    string s_bnd_type, s_shape_type;

    istringstream iss;

    if( TagCheck(sLine, FileBlockName::StartBoundaryEdge()) ){
        mpLogger->Info(Utility::LoggerMode::MWDebug, "FileReaderBoundaryEdge", sLine);
        
        sLine= getLineSt(ifs);
        iss.clear();
        iss.str(sLine.c_str());

        // 境界ID, 境界種類, MeshID, 境界節点数, 境界辺数
        iss >> bnd_id >> s_bnd_type >> mesh_id >> numOfBNode >> numOfBEdge;

        bnd_type= IntBndType(s_bnd_type);

        uint ibnode;
        // 境界節点
        for(ibnode=0; ibnode < numOfBNode; ibnode++){
            sLine= getLineSt(ifs);
            iss.clear();
            iss.str(sLine);

            // BNodeID, NodeID, X, Y, Z  :=> x,y,z は未使用
            iss >> bnode_id >> node_id >> x >> y >> z;

            mpFactory->GeneBoundaryEdgeNode(mgLevel, bnd_id, bnd_type, mesh_id, node_id, bnode_id);
        };

        // Edge
        //
        while(!ifs.eof()){
            sLine= getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryEdge()) ) break;

            iss.clear();
            iss.str(sLine);

            // 形状タイプ, BEdgeID, ElementID, EdgeID, DOF,     BNodeID, BNodeID, Value
            iss >> s_shape_type >> bedge_id >> elem_id >> ent_id >> dof;

            vBNodeID.clear();
            if(s_shape_type=="Line" || s_shape_type=="Beam"){
                shape_type= pmw::ElementType::Beam;
                vBNodeID.resize(2);

                iss >> vBNodeID[0] >> vBNodeID[1] >> val;// BNodeID, BNodeID, Value
            }else{
                ;//TODO: Logger => Error
            }

            mpFactory->GeneBoundaryEdge(mgLevel, bnd_id, bnd_type, shape_type,
                                        mesh_id, elem_id, ent_id, vBNodeID, bedge_id, dof, val);
        };//while end

        mpFactory->initEdgeAggregate(mgLevel, mesh_id, bnd_id);//BNode,BEdgeを全てセットした後に呼び出す

        return true;
    }else{
        return false;
    }
}


