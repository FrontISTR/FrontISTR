//
// FileReaderBoundaryVolumeMesh.cpp
//
//              2010.04.28
//              k.Takeda
#include "FileReaderBoundaryVolumeMesh.h"
using namespace FileIO;

CFileReaderBoundaryVolumeMesh::CFileReaderBoundaryVolumeMesh()
{
    ;
}
CFileReaderBoundaryVolumeMesh::~CFileReaderBoundaryVolumeMesh()
{
    ;
}

bool CFileReaderBoundaryVolumeMesh::Read(ifstream& ifs, string& sLine)
{
    uint mgLevel(0);
    uint bnd_id, bnd_type, mesh_id, numOfBoundary, numOfDOF;
    vuint vDOF;
    string s_bnd_type;
    istringstream iss;

    if( TagCheck(sLine, FileBlockName::StartBoundaryVolumeMesh()) ){

        while(!ifs.eof()){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryVolumeMesh()) ) break;

            iss.clear();
            iss.str(sLine);

            iss >> mesh_id >> numOfBoundary;
            
            mpFactory->reserveBoundaryVolumeMesh(mgLevel, mesh_id, numOfBoundary);

            uint ibound;
            for(ibound=0; ibound < numOfBoundary; ibound++){
                sLine= getLineSt(ifs);

                iss.clear();
                iss.str(sLine);
                
                iss >> bnd_id >> s_bnd_type >> numOfDOF;

                bnd_type= IntBndType(s_bnd_type);

                vDOF.clear(); vDOF.resize(numOfDOF);
                for(uint i=0; i < numOfDOF; i++){
                    iss >> vDOF[i];
                }

                mpFactory->GeneBoundaryVolumeMesh(mgLevel, mesh_id, bnd_id, bnd_type, numOfDOF, vDOF);
            };
        };
        return true;
    }else{
        return false;
    }
}








