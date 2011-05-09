/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderChunk.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include <vector>
#include "FileReaderChunk.h"
#include "FileReaderContactMesh.h"
using namespace FileIO;
CFileReaderChunk::CFileReaderChunk()
{
    msCntFileName= "mw3.cnt";
    mpLogger = Utility::CLogger::Instance();
    mvReader.reserve(15);
    mvReader.push_back(new CFileReaderNode());
    mvReader.push_back(new CFileReaderElement());
    mvReader.push_back(new CFileReaderAssyModel);
    mvReader.push_back(new CFileReaderRefine);
    mvReader.push_back(new CFileReaderBoundaryNode);
    mvReader.push_back(new CFileReaderBoundaryFace);
    mvReader.push_back(new CFileReaderBoundaryVolume);
    mvReader.push_back(new CFileReaderMaterial);
    mvReader.push_back(new CFileReaderCommMesh);
    mvReader.push_back(new CFileReaderCommNode);
    mvReader.push_back(new CFileReaderCommElement);
    mvReader.push_back(new CFileReaderContactMesh);
    mvReader.push_back(new CFileReaderCommMesh2);
    mvReader.push_back(new CFileReaderCommFace);
    mvReader.push_back(new CFileReaderCommNodeCM2);
}
CFileReaderChunk::~CFileReaderChunk()
{
    for_each(mvReader.begin(), mvReader.end(), pmw::DeleteObject());
}
void CFileReaderChunk::setFactory(pmw::CMeshFactory *pFactory)
{
    uint i;
    for(i=0; i<mvReader.size(); i++){
            mvReader[i]->setFactory(pFactory);
    };
}
void CFileReaderChunk::Read(string filename)
{
    char c_Line[BUFFERLENGTH];
    string s_Line;
    mpLogger->Info(Utility::LoggerMode::MWDebug,"Input Filename => ",filename);
    ifstream ifs(filename.c_str(),ios::in);
    if(ifs){
        uint i;
        while(!ifs.eof()){
            ifs.getline(c_Line,sizeof(c_Line),'\n');
            s_Line = c_Line;
            for(i=0; i < mvReader.size(); i++){
                    mvReader[i]->Read(ifs, s_Line);
            };
        };
    }else{
        mpLogger->Info(Utility::LoggerMode::Error, "File not fount, filename => ", filename);
    }
    ifs.close();
}
void CFileReaderChunk::ReadCnt()
{
    ifstream ifs;
    char   c_Line[BUFFERLENGTH];
    string s_Line;
    ifs.open(msCntFileName.c_str(), ios::in);
    if(ifs){
        while(!ifs.eof()){
            ifs.getline(c_Line,sizeof(c_Line),'\n');
            s_Line = c_Line;
            mpCntReader->Read(ifs, s_Line);
        };
    }else{
        mpLogger->Info(Utility::LoggerMode::Error, "File not fount, filename => ", msCntFileName);
    }
    ifs.close();
}
void CFileReaderChunk::setPath(string& filepath)
{
    string sFullPathName;
    sFullPathName = filepath + msCntFileName;
    msCntFileName = sFullPathName;
}
