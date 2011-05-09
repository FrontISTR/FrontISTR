//
//  FileReaderBoundaryFaceMesh.cpp
//
//              2010.04.28
//              k.Takeda
#include "FileReaderBoundaryFaceMesh.h"
using namespace FileIO;


CFileReaderBoundaryFaceMesh::CFileReaderBoundaryFaceMesh()
{
    ;
}
CFileReaderBoundaryFaceMesh::~CFileReaderBoundaryFaceMesh()
{
    ;
}

bool CFileReaderBoundaryFaceMesh::Read(ifstream& ifs, string& sLine)
{
    uint mgLevel(0);
    uint bnd_id, bnd_type, mesh_id, numOfBoundary, numOfDOF;
    vuint vDOF;
    string s_bnd_type;
    istringstream iss;

    if( TagCheck(sLine, FileBlockName::StartBoundaryFaceMesh()) ){

        while(!ifs.eof()){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryFaceMesh()) ) break;

            iss.clear();
            iss.str(sLine);

            iss >> mesh_id >> numOfBoundary;

            mpFactory->reserveBoundaryFaceMesh(mgLevel, mesh_id, numOfBoundary);

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

                mpFactory->GeneBoundaryFaceMesh(mgLevel, mesh_id, bnd_id, bnd_type, numOfDOF, vDOF);
            };
        };
        return true;
    }else{
        return false;
    }
}





