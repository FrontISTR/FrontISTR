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
        
        //BVolume
        //----
        while(!ifs.eof()){
            sLine= getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryVolume())) break;

            iss.clear();
            iss.str(sLine);

            iss >> s_shape_type >> bvol_id >> elem_id >> ent_id >> dof;

            vBNodeID.clear();
            if(s_shape_type=="Hexa"){
                shape_type= pmw::ElementType::Hexa;
                vBNodeID.resize(pmw::NumberOfVertex::Hexa());

                for(ibnode=0; ibnode < pmw::NumberOfVertex::Hexa(); ibnode++){
                    iss >> vBNodeID[ibnode];
                };
            }else if(s_shape_type=="Tetra"){
                shape_type= pmw::ElementType::Tetra;
                vBNodeID.resize(pmw::NumberOfVertex::Tetra());

                for(ibnode=0; ibnode < pmw::NumberOfVertex::Tetra(); ibnode++){
                    iss >> vBNodeID[ibnode];
                };
            }else if(s_shape_type=="Prism"){
                shape_type= pmw::ElementType::Prism;
                vBNodeID.resize(pmw::NumberOfVertex::Prism());

                for(ibnode=0; ibnode < pmw::NumberOfVertex::Prism(); ibnode++){
                    iss >> vBNodeID[ibnode];
                };
            }else{
                ;//TODO: Logger => Error
            }
            iss >> val;

            ////debug
            //cout << "FileReaderBoundaryVolume::Read, type_name=" << s_shape_type << ", shape_type= " << shape_type << endl;

            mpFactory->GeneBoundaryVolume(mgLevel, bnd_id, bnd_type, shape_type,
                                          mesh_id, elem_id, vBNodeID, bvol_id, dof, val);
        };

        mpFactory->initVolumeAggregate(mgLevel, mesh_id, bnd_id);//BNode,BVolumeを全てセットした後に呼び出す

        return true;
    }else{
        return false;
    }
}


