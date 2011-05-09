/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderChunk.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef FILE_READER_HH_BBAE6774_525B_49ec_8315_9FD9F5052070
#define FILE_READER_HH_BBAE6774_525B_49ec_8315_9FD9F5052070
#include "TypeDef.h"
#include "Logger.h"
#include "FileReaderNode.h"
#include "FileReaderElement.h"
#include "FileReaderAssyModel.h"
#include "FileReaderRefine.h"
#include "FileReaderMaterial.h"
#include "FileReaderCommMesh.h"
#include "FileReaderCommNode.h"
#include "FileReaderCommElement.h"
#include "FileReaderBoundaryNode.h"
#include "FileReaderBoundaryFace.h"
#include "FileReaderBoundaryVolume.h"
#include "FileReaderCnt.h"
#include "FileReaderContactMesh.h"
#include "FileReaderCommMesh2.h"
#include "FileReaderCommFace.h"
#include "FileReaderCommNode_CM2.h"
namespace FileIO{
class CFileReaderChunk
{
public:
    CFileReaderChunk();
    CFileReaderChunk(pmw::CMeshFactory *pFactory);
    virtual ~CFileReaderChunk();
private:
    vector<CFileReader*> mvReader;
    CFileReaderCnt    *mpCntReader;
    string  msCntFileName;
    Utility::CLogger *mpLogger;
public:
    void setCntReader(CFileReaderCnt* pReader){ mpCntReader= pReader;}
    void ReadCnt();
    void Read(string filename);
    void setPath(string& filepath);
    void setFactory(pmw::CMeshFactory* pFactory);
    void setLogger(Utility::CLogger *pLogger){ mpLogger = pLogger;}
};
}
#endif
