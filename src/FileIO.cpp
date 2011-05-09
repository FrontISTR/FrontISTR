//
//  FileIO.cpp
//              2008.12.08
//              2008.12.08
//              k.Takeda

#include "FileIO.h"
using namespace FileIO;

// construct && destruct
// --
CFileIO::CFileIO()
{
    // CntReaderのセット
    moReader.setCntReader(&moCntReader);
}
CFileIO::~CFileIO()
{
    ;
}


// pathの設定
//
void CFileIO::setPathName(const char* cpath)
{
    msPathName= cpath;
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

// メッシュファイルのベースネームが記述されているHEC_MW3.cntを読み込む.
//  => *.mshベースネームを取得
void CFileIO::ReadCntFile()
{
    moReader.setPath(msPathName);
    moReader.ReadCnt();
}


// HEC_MW3 の標準入力ファイル
//
void CFileIO::ReadFile(string filename)
{
    moReader.Read(filename);
}


// Solution Type
void CFileIO::setSolutionType(const uint& nSolutionType)
{
    mnSolutionType = nSolutionType;

    moWriter.setSolutionType(nSolutionType);
}

// HEC_MW3 の標準出力ファイル
//
void CFileIO::WriteFile(string filename, const uint& numOfLevel)
{
    moWriter.Write(filename, numOfLevel);
}





