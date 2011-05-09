//
// FileReaderBoundaryNodeMesh.cpp
//
//          2010.04.28
//          k.Takeda
#include "FileReaderBoundaryNodeMesh.h"
using namespace FileIO;

CFileReaderBoundaryNodeMesh::CFileReaderBoundaryNodeMesh()
{
    ;
}
CFileReaderBoundaryNodeMesh::~CFileReaderBoundaryNodeMesh()
{
    ;
}

bool CFileReaderBoundaryNodeMesh::Read(ifstream& ifs, string& sLine)
{
    istringstream iss;
    uint mgLevel(0), bnd_id, bnd_type, mesh_id, numOfBoundary;
    string s_bnd_type;

    if( TagCheck(sLine, FileBlockName::StartBoundaryNodeMesh()) ){

        while(!ifs.eof()){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryNodeMesh()) ) break;

            iss.clear();
            iss.str(sLine.c_str());

            iss >> mesh_id >> numOfBoundary;

            mpFactory->reserveBoundaryNodeMesh(mgLevel, mesh_id, numOfBoundary);
            
            uint ibound;
            for(ibound=0; ibound < numOfBoundary; ibound++){
                sLine= getLineSt(ifs);
                iss.clear();
                iss.str(sLine);
                
                iss >> bnd_id >> s_bnd_type;

                //cout << "BoundaryNodeMesh bnd_type = " << s_bnd_type << endl;

                bnd_type= IntBndType(s_bnd_type);
                
                mpFactory->GeneBoundaryNodeMesh(mgLevel, mesh_id, bnd_id, bnd_type);
            };
            
        };
        return true;
    }else{
        return false;
    }
}








