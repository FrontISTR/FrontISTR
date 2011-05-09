//
//  FileReaderBoundaryNode.cpp
//
//
//
//                          2009.05.22
//                          2009.05.22
//                          k.Takeda
#include "FileReaderBoundaryNode.h"
using namespace FileIO;


CFileReaderBoundaryNode::CFileReaderBoundaryNode()
{
    ;
}

CFileReaderBoundaryNode::~CFileReaderBoundaryNode()
{
    ;
}

bool CFileReaderBoundaryNode::Read(ifstream& ifs, string& sLine)
{
    uint bnode_id, node_id, dof;
    uint mgLevel(0);
    uint bnd_id, bnd_type, numOfBNode, mesh_id;
    double val, x, y, z;
    string s_bnd_type;
    
    istringstream iss;
    
    if(TagCheck(sLine, FileBlockName::StartBoundaryNode()) ){
        
        sLine = getLineSt(ifs);
        iss.clear();
        iss.str(sLine);

        // NodeBoundaryID, 境界種類, MeshID, 境界節点数
        iss >> bnd_id >> s_bnd_type >> mesh_id >> numOfBNode;

        // 境界タイプ文字列を uint に変換
        bnd_type= IntBndType(s_bnd_type);

        while(!ifs.eof()){
            sLine= getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryNode())) break;
            iss.clear();
            iss.str(sLine);
            // BNodeID, NodeID, x, y, z, DOF, Value :=> x,y,z は未使用
            iss >> bnode_id >> node_id >> x >> y >> z >> dof >> val;

            mpFactory->GeneBoundaryNode(mgLevel, bnd_id, bnd_type, mesh_id, node_id, bnode_id, dof, val);
        };

        return true;
    }else{
        return false;
    }
}
















