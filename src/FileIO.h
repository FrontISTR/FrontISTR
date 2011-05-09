//
//	FileIO.h
//
//			2008.12.08
//			2008.12.08
//			k.Takeda
#ifndef FILE_IO_HH_4D6193F1_360F_4f01_A219_725098D1C2FC
#define FILE_IO_HH_4D6193F1_360F_4f01_A219_725098D1C2FC

#include "CommonStd.h"
#include "TypeDef.h"

#include "FileReaderChunk.h"
#include "FileWriterChunk.h"

#include "FileReaderCnt.h"//mw3.cnt

namespace FileIO{
class CFileIO
{
private:
    CFileIO();
public:
    static CFileIO* Instance(){
        static CFileIO  file_io;
        return &file_io;
    }
    virtual ~CFileIO();

private:
    string msPathName;// ファイルのパス

    CFileReaderChunk moReader;
    CFileWriterChunk moWriter;

    CFileReaderCnt   moCntReader;

public:
    void setFactory(pmw::CMeshFactory *pFactory);
    void setLogger(Utility::CLogger *pLogger);

    // cntファイルからメッシュのベースネームを取得
    void ReadCntFile();
    string& getMeshFileBaseName(){return moCntReader.getMeshFileBaseName();}
    string& getPathName(){ return msPathName;}


    // HEC_MW3 書式のメッシュデータ
    void ReadFile(string filename);
    void WriteFile(string filename, const uint& nmgLevel);
};
}
#endif
