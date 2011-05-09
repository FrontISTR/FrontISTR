//
//  FileReaderCnt.cpp
//
//
//
//                  2009.09.22
//                  2009.09.22
//                  k.Takeda
#include "FileReaderCnt.h"
#include "FileBlockName.h"
using namespace FileIO;

// construct & destruct
// --
CFileReaderCnt::CFileReaderCnt()
{
    ;
}
CFileReaderCnt::~CFileReaderCnt()
{
    ;
}

// ＊他のReadメソッドとは,ReaderChunkでの使い方が異なっているので注意.
// --
bool CFileReaderCnt::Read(ifstream& ifs, string& sLine)
{
    // MeshFileName
    if(TagCheck(sLine, FileBlockName::StartMeshFileName()) ){

        while(true){
            sLine = getLineSt(ifs);
            // Endブロックによるブレーク
            if(sLine==FileBlockName::EndMeshFileName()) break;

            istringstream iss(sLine.c_str());
            iss >> msMeshFileBaseName;

            mpLogger->Info(Utility::LoggerMode::Debug,"ベースファイル名=> ",msMeshFileBaseName);
        };
        return true;
    }else{
        return false;
    }
}



















