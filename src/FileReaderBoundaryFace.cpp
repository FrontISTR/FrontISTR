
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
    uint mesh_id, elem_id, ent_id, bface_id, dof, bnode_id, node_id, numOfBNode, numOfBFace;
    uint mgLevel(0);// mgLevel=="0" :コースグリッド
    uint bface_shape;
    vuint vBNodeID;
    double val, x, y, z;
    uint bnd_id, bnd_type;
    string s_bnd_type, s_bface_shape;
    
    istringstream iss;

    if(TagCheck(sLine, FileBlockName::StartBoundaryFace()) ){

        sLine= getLineSt(ifs);
        iss.clear();
        iss.str(sLine);

        // 境界ID, 境界種類, MeshID, 境界節点数, 境界面数
        iss >> bnd_id >> s_bnd_type >> mesh_id >> numOfBNode >> numOfBFace;
        
        bnd_type= IntBndType(s_bnd_type);
        
        
        //cout << "FileReaderBoundaryFace::Read, numOfBNode= " << numOfBNode << endl;

        uint ibnode;
        for(ibnode=0; ibnode < numOfBNode; ibnode++){
            sLine= getLineSt(ifs);
            iss.clear();
            iss.str(sLine);

            // BNodeID, NodeID, X, Y, Z :=> X,Y,Z は未使用
            iss >> bnode_id >> node_id >> x >> y >> z;

            //cout << "FileReaderBoundaryFace::Read, bnode_id= " << bnode_id << endl;

            mpFactory->GeneBoundaryFaceNode(mgLevel, bnd_id, bnd_type, mesh_id, node_id, bnode_id);
        };
        
        while(!ifs.eof()){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryFace())) break;

            // Face
            //
            iss.clear();
            iss.str(sLine);

            // 形状タイプ, BFaceID, ElementID, EntityID, DOF, BNodeID.....BNodeID, Value
            iss >> s_bface_shape >> bface_id >> elem_id >> ent_id >> dof;

            // BNodeID....BNodeIDの処理
            vBNodeID.clear();
            if(s_bface_shape=="Quad"){
                bface_shape= pmw::ElementType::Quad;
                vBNodeID.resize(4);

                for(ibnode=0; ibnode < 4; ibnode++){
                    iss >> bnode_id;
                    vBNodeID[ibnode]= bnode_id;
                };

            }else if(s_bface_shape=="Triangle"){
                bface_shape= pmw::ElementType::Triangle;
                vBNodeID.resize(3);

                for(ibnode=0; ibnode < 3; ibnode++){
                    iss >> bnode_id;
                    vBNodeID[ibnode]= bnode_id;
                };

            }else{
                //TODO:Logger => throw
            }

            iss >> val;// 境界値

            //cout << "FileReaderBoundaryFace::read, GeneBoundaryFace A" << endl;

            mpFactory->GeneBoundaryFace(mgLevel, bnd_id, bnd_type, bface_shape,
                                            mesh_id, elem_id, ent_id, vBNodeID, bface_id, dof, val);

        };

        mpFactory->initFaceAggregate(mgLevel, mesh_id, bnd_id);//BNode,BFaceを全てセットしたあとに呼び出す

        return true;
    }else{
        return false;
    }
}


