//
//  FileReaderBoundaryVolume.cpp
//
//
//
//                          2009.05.22
//                          2009.05.22
//                          k.Takeda
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


bool CFileReaderBoundaryVolume::Read(ifstream& ifs, string& sLine)
{
    uint bnd_id, bnd_type, mesh_id, numOfBNode, numOfBVol;
    uint bnode_id, node_id, dof, bvol_id, elem_id, ent_id, shape_type;
    vuint vBNodeID;
    double x, y, z, val;
    string  s_bnd_type, s_shape_type;
    uint mgLevel(0);

    istringstream iss;
    if(TagCheck(sLine, FileBlockName::StartBoundaryVolume()) ){
        sLine = getLineSt(ifs);

        iss.clear();
        iss.str(sLine);

        iss >> bnd_id >> s_bnd_type >> mesh_id >> numOfBNode >> numOfBVol;

        bnd_type= IntBndType(s_bnd_type);

        uint ibnode;
        //BNode
        //----
        for(ibnode=0; ibnode < numOfBNode; ibnode++){
            sLine= getLineSt(ifs);

            iss.clear();
            iss.str(sLine);

            iss >> bnode_id >> node_id >> x >> y >> z;// x,y,z は,未使用

            mpFactory->GeneBoundaryVolumeNode(mgLevel, bnd_id, bnd_type, mesh_id, node_id, bnode_id);
        };

        mpFactory->resizeVolumeAggregate(mgLevel, mesh_id, bnd_id);
        
        //BVolume
        //----
        while(!ifs.eof()){
            sLine= getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryVolume())) break;

            iss.clear();
            iss.str(sLine);

            iss >> s_shape_type >> bvol_id >> elem_id >> ent_id >> dof;

            shape_type = IntElemType(s_shape_type);

            vBNodeID.clear();
            switch(shape_type){
                case(pmw::ElementType::Hexa):
                    vBNodeID.resize(pmw::NumberOfNode::Hexa());
                    for(ibnode=0; ibnode < pmw::NumberOfNode::Hexa(); ibnode++){
                        iss >> vBNodeID[ibnode];
                    };
                    break;
                case(pmw::ElementType::Hexa2):
                    vBNodeID.resize(pmw::NumberOfNode::Hexa2());
                    for(ibnode=0; ibnode < pmw::NumberOfNode::Hexa2(); ibnode++){
                        iss >> vBNodeID[ibnode];
                    };
                    break;
                case(pmw::ElementType::Tetra):
                    vBNodeID.resize(pmw::NumberOfNode::Tetra());
                    for(ibnode=0; ibnode < pmw::NumberOfNode::Tetra(); ibnode++){
                        iss >> vBNodeID[ibnode];
                    };
                    break;
                case(pmw::ElementType::Tetra2):
                    vBNodeID.resize(pmw::NumberOfNode::Tetra2());
                    for(ibnode=0; ibnode < pmw::NumberOfNode::Tetra2(); ibnode++){
                        iss >> vBNodeID[ibnode];
                    };
                    break;
                case(pmw::ElementType::Prism):
                    vBNodeID.resize(pmw::NumberOfNode::Prism());
                    for(ibnode=0; ibnode < pmw::NumberOfNode::Prism(); ibnode++){
                        iss >> vBNodeID[ibnode];
                    };
                    break;
                case(pmw::ElementType::Prism2):
                    vBNodeID.resize(pmw::NumberOfNode::Prism2());
                    for(ibnode=0; ibnode < pmw::NumberOfNode::Prism2(); ibnode++){
                        iss >> vBNodeID[ibnode];
                    };
                    break;
                default:
                    //TODO: Logger->Error
                    break;
            }

            iss >> val;

            //debug
            ///cout << "FileReaderBoundaryVolume::Read, type_name=" << s_shape_type << ", shape_type= " << shape_type << endl;

            mpFactory->GeneBoundaryVolume(mgLevel, bnd_id, bnd_type, shape_type,
                                          mesh_id, elem_id, vBNodeID, bvol_id, dof, val);
        };

        mpFactory->initVolumeAggregate(mgLevel, mesh_id, bnd_id);//BNode,BVolumeを全てセットした後に呼び出す

        return true;
    }else{
        return false;
    }
}


