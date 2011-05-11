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
    uiint mgLevel(0);
    uiint bnd_id, bnd_type, mesh_id, numOfBoundary, numOfDOF;
    vuint vDOF;
    string s_bnd_type, s_bnd_name;
    istringstream iss;

    if( TagCheck(sLine, FileBlockName::StartBoundaryEdgeMesh()) ){
        
        while(!ifs.eof()){
            sLine = getLine(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryEdgeMesh()) ) break;

            iss.clear();
            iss.str(sLine);

            iss >> mesh_id >> numOfBoundary;

            mpFactory->reserveBoundaryEdgeMesh(mgLevel, mesh_id, numOfBoundary);

            uiint ibound;
            for(ibound=0; ibound < numOfBoundary; ibound++){
                sLine= getLine(ifs);

                iss.clear();
                iss.str(sLine);

                iss >> bnd_id >> s_bnd_type >> s_bnd_name >> numOfDOF;

                bnd_type = IntBndType(s_bnd_type);

                vDOF.clear(); vDOF.resize(numOfDOF);
                for(uiint i=0; i < numOfDOF; i++){
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

bool CFileReaderBoundaryEdgeMesh::Read_bin(ifstream& ifs)
{
    CFileReaderBinCheck *pBinCheck= CFileReaderBinCheck::Instance();
    bool bOrder= pBinCheck->isByteOrderSwap();

    //BinCheckのサイズ指定との整合性
    bool b32, bCheck;
    string sClassName("FileReaderBoundaryEdgeMesh");

    if( !Check_IntSize(b32, bCheck, sClassName) ) return false;

    char cHead='B';
    if( !TagCheck_Bin(ifs, bCheck, cHead, FileBlockName::StartBoundaryEdgeMesh(), FileBlockName::BoundaryEdgeMesh_Len())) return false;


    uiint mgLevel(0);
    uiint bnd_id, bnd_type, mesh_id, nNumOfBoundary, nNumOfDOF;
    vuint vDOF;
    string s_bnd_type, s_bnd_name;

    while(!ifs.eof()){
        // MeshID 境界Mesh数
        ifs.read((char*)&mesh_id, sizeof(uiint));        if(bOrder) pBinCheck->ByteOrderSwap(mesh_id);
        ifs.read((char*)&nNumOfBoundary, sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(nNumOfBoundary);

        mpFactory->reserveBoundaryEdgeMesh(mgLevel, mesh_id, nNumOfBoundary);

        for(uiint ibound=0; ibound < nNumOfBoundary; ibound++){
            // 境界ID   境界タイプ  境界名   自由度
            ifs.read((char*)&bnd_id, sizeof(uiint));  if(bOrder) pBinCheck->ByteOrderSwap(bnd_id);
            Read_BndType(ifs, s_bnd_type);//Dirichlet, Neumann
            Read_AnyName(ifs, s_bnd_name);//英数字の任意名
            ifs.read((char*)&nNumOfDOF, sizeof(uiint));  if(bOrder) pBinCheck->ByteOrderSwap(nNumOfDOF);
            
            bnd_type = IntBndType(s_bnd_type);

            vDOF.clear(); vDOF.resize(nNumOfDOF);
            for(uiint i=0; i < nNumOfDOF; i++){
                ifs.read((char*)&vDOF[i], sizeof(uiint)); if(bOrder) pBinCheck->ByteOrderSwap(vDOF[i]);
            };

            mpFactory->GeneBoundaryEdgeMesh(mgLevel, mesh_id, bnd_id, bnd_type, s_bnd_name, nNumOfDOF, vDOF);
        };
    };

    return true;
}




