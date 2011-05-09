/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileIO.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef FILE_IO_HH_4D6193F1_360F_4f01_A219_725098D1C2FC
#define FILE_IO_HH_4D6193F1_360F_4f01_A219_725098D1C2FC
#include "CommonStd.h"
#include "TypeDef.h"
#include "FileReaderChunk.h"
#include "FileWriterChunk.h"
#include "FileReaderCnt.h"
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
    string msPathName;
    CFileReaderChunk moReader;
    CFileWriterChunk moWriter;
    CFileReaderCnt   moCntReader;
public:
    void setFactory(pmw::CMeshFactory *pFactory);
    void setLogger(Utility::CLogger *pLogger);
    void ReadCntFile();
    string& getMeshFileBaseName(){return moCntReader.getMeshFileBaseName();}
    void setPathName(const char* path);
    string& getPathName(){ return msPathName;}
    void ReadFile(string filename);
    void WriteFile(string filename, const uint& nmgLevel);
};
}
#endif
