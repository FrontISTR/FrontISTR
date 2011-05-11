//
// FileReaderBoundaryVolumeMesh.cpp
//
//              2010.04.28
//              k.Takeda
#include "FileReaderBoundaryVolumeMesh.h"
using namespace FileIO;
using namespace boost;

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
    uiint mgLevel(0);
    uiint bnd_id, bnd_type, mesh_id, numOfBoundary, numOfDOF;
    vuint vDOF;
    string s_bnd_type, s_bnd_name;
    istringstream iss;

    if( TagCheck(sLine, FileBlockName::StartBoundaryVolumeMesh()) ){

        while(!ifs.eof()){
            sLine = getLine(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryVolumeMesh()) ) break;

            iss.clear();
            iss.str(sLine);

            iss >> mesh_id >> numOfBoundary;
            
            mpFactory->reserveBoundaryVolumeMesh(mgLevel, mesh_id, numOfBoundary);

            uiint ibound;
            for(ibound=0; ibound < numOfBoundary; ibound++){
                sLine= getLine(ifs);

                iss.clear();
                iss.str(sLine);
                
                iss >> bnd_id >> s_bnd_type >> s_bnd_name >> numOfDOF;

                bnd_type= IntBndType(s_bnd_type);

                vDOF.clear(); vDOF.resize(numOfDOF);
                for(uiint i=0; i < numOfDOF; i++){
                    iss >> vDOF[i];
                }

                mpFactory->GeneBoundaryVolumeMesh(mgLevel, mesh_id, bnd_id, bnd_type, s_bnd_name, numOfDOF, vDOF);
            };
        };
        return true;
    }else{
        return false;
    }
}


bool CFileReaderBoundaryVolumeMesh::Read_bin(ifstream& ifs)
{
    CFileReaderBinCheck *pBinCheck= CFileReaderBinCheck::Instance();
    bool bOrder= pBinCheck->isByteOrderSwap();

    //BinCheckのサイズ指定との整合性
    bool b32, bCheck;
    string sClassName("FileReaderBoundaryVolumeMesh");

    if( !Check_IntSize(b32, bCheck, sClassName) ) return false;

    char cHead='B';
    if( !TagCheck_Bin(ifs, bCheck, cHead, FileBlockName::StartBoundaryVolumeMesh(), FileBlockName::BoundaryVolumeMesh_Len())) return false;

    uiint mgLevel(0);
    uiint bnd_id, bnd_type, mesh_id, nNumOfBoundary, nNumOfDOF;
    vuint vDOF;
    string s_bnd_type, s_bnd_name;

    while(!ifs.eof()){
        if( CFileReader::Check_End(ifs) ) break;

        ifs.read((char*)&mesh_id, sizeof(uiint));  if(bOrder) pBinCheck->ByteOrderSwap(mesh_id);
        ifs.read((char*)&nNumOfBoundary, sizeof(uiint));  if(bOrder) pBinCheck->ByteOrderSwap(nNumOfBoundary);

        mpFactory->reserveBoundaryVolumeMesh(mgLevel, mesh_id, nNumOfBoundary);

        uiint ibound;
        for(ibound=0; ibound < nNumOfBoundary; ibound++){

            ifs.read((char*)&bnd_id, sizeof(uiint));    if(bOrder) pBinCheck->ByteOrderSwap(bnd_id);
            CFileReader::Read_BndType(ifs, s_bnd_type);
            CFileReader::Read_AnyName(ifs, s_bnd_name);
            ifs.read((char*)&nNumOfDOF, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(nNumOfDOF);

            bnd_type= IntBndType(s_bnd_type);

            vDOF.clear(); vDOF.resize(nNumOfDOF);
            for(uiint i=0; i < nNumOfDOF; i++){
                ifs.read((char*)&vDOF[i], sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(vDOF[i]);
            }

            mpFactory->GeneBoundaryVolumeMesh(mgLevel, mesh_id, bnd_id, bnd_type, s_bnd_name, nNumOfDOF, vDOF);
        };
    };

    return true;
}





