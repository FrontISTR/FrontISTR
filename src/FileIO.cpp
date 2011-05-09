//
//  FileIO.cpp
//              2008.12.08
//              2008.12.08
//              k.Takeda

#include "FileIO.h"
using namespace FileIO;

CFileIO::CFileIO()
{
    msPathName = "";
}

CFileIO::~CFileIO()
{
}

// Factory を 各Reader に設置
//
void CFileIO::setFactory(pmw::CMeshFactory* pFactory)
{
    moReader.setFactory(pFactory);
}

// Logger を 各Reader に設置
//
void CFileIO::setLogger(Utility::CLogger* pLogger)
{
    moReader.setLogger(pLogger);
}

// HEC_MW3 の標準入力ファイル
//
void CFileIO::ReadFile(string filename)
{
    moReader.Read(filename);
}

// HEC_MW3 の標準出力ファイル
//
void CFileIO::WriteFile(string filename, const uint& numOfLevel)
{
    moWriter.Write(filename, numOfLevel);
}





