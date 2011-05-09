//
// FileReaderBoundaryNodeMesh.cpp
//
//          2010.04.28
//          k.Takeda
#include "FileReaderBoundaryNodeMesh.h"
using namespace FileIO;
using namespace boost;
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
    string s_bnd_type, s_bnd_name("nameless");

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

                //iss.clear();
                //iss.str(sLine);
                //
                //iss >> bnd_id >> s_bnd_type;

                // boost トークン分割
                // ----
                char_separator<char> sep(" \t\n");
                tokenizer< char_separator<char> > tokens(sLine, sep);

                uint nCount(0);
                typedef tokenizer< char_separator<char> >::iterator Iter;
                for(Iter it=tokens.begin(); it != tokens.end(); ++it){
                    string str = *it;
                    if(nCount==0){ bnd_id = atoi(str.c_str());}
                    if(nCount==1){ s_bnd_type = str;}
                    if(nCount==2){ s_bnd_name = str;}
                    nCount++;
                };

                bnd_type= IntBndType(s_bnd_type);
                
                mpFactory->GeneBoundaryNodeMesh(mgLevel, mesh_id, bnd_id, bnd_type, s_bnd_name);
            };
            
        };
        return true;
    }else{
        return false;
    }
}








