//
//  FileReaderBoundaryEdgeMesh.cpp
//
//                  2010.04.28
//                  k.Takeda
#include "FileReaderBoundaryEdgeMesh.h"
using namespace FileIO;
using namespace boost;

CFileReaderBoundaryEdgeMesh::CFileReaderBoundaryEdgeMesh()
{
    ;
}
CFileReaderBoundaryEdgeMesh::~CFileReaderBoundaryEdgeMesh()
{
    ;
}


// 境界EdgeMesh => 境界種類数ぶんの領域確保
//
bool CFileReaderBoundaryEdgeMesh::Read(ifstream& ifs, string& sLine)
{
    uint mgLevel(0);
    uint bnd_id, bnd_type, mesh_id, numOfBoundary, numOfDOF;
    vuint vDOF;
    string s_bnd_type, s_bnd_name;
    istringstream iss;

    if( TagCheck(sLine, FileBlockName::StartBoundaryEdgeMesh()) ){
        
        while(!ifs.eof()){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryEdgeMesh()) ) break;

            iss.clear();
            iss.str(sLine);

            iss >> mesh_id >> numOfBoundary;

            mpFactory->reserveBoundaryEdgeMesh(mgLevel, mesh_id, numOfBoundary);

            uint ibound;
            for(ibound=0; ibound < numOfBoundary; ibound++){
                sLine= getLineSt(ifs);

                iss.clear();
                iss.str(sLine);

                iss >> bnd_id >> s_bnd_type >> s_bnd_name >> numOfDOF;

                bnd_type = IntBndType(s_bnd_type);

                vDOF.clear(); vDOF.resize(numOfDOF);
                for(uint i=0; i < numOfDOF; i++){
                    iss >> vDOF[i];
                };

                mpFactory->GeneBoundaryEdgeMesh(mgLevel, mesh_id, bnd_id, bnd_type, s_bnd_name, numOfDOF, vDOF);
            };
        };
        return true;
    }else{
        return false;
    }
}






