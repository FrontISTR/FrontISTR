/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileIO.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FileIO.h"
using namespace FileIO;
CFileIO::CFileIO()
{
    moReader.setCntReader(&moCntReader);
}
CFileIO::~CFileIO()
{
    ;
}
void CFileIO::setPathName(const char* cpath)
{
    msPathName= cpath;
}
void CFileIO::setFactory(pmw::CMeshFactory* pFactory)
{
    moReader.setFactory(pFactory);
}
void CFileIO::setLogger(Utility::CLogger* pLogger)
{
    moReader.setLogger(pLogger);
}
void CFileIO::ReadCntFile()
{
    moReader.setPath(msPathName);
    moReader.ReadCnt();
}
void CFileIO::ReadFile(string filename)
{
    moReader.Read(filename);
}
void CFileIO::WriteFile(string filename, const uint& numOfLevel)
{
    moWriter.Write(filename, numOfLevel);
}
