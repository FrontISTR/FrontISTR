//
//  FileReaderChank.cpp
//
//			2009.05.01
//			2008.12.09
//			k.Takeda
#include "FileReaderChunk.h"
#include "FileReaderContactMesh.h"
using namespace FileIO;

// Constructor
//
CFileReaderChunk::CFileReaderChunk()
{
    //cntファイル名:とりあえず"mw3.cnt"としてある.<= 2009.09.22
    //--
    msCntFileName= "mw3.cnt";

    mpLogger = Utility::CLogger::Instance();

    mvReader.reserve(22);

    mvReader.push_back(new CFileReaderNode());
    mvReader.push_back(new CFileReaderElement());
    mvReader.push_back(new CFileReaderAssyModel);
    mvReader.push_back(new CFileReaderRefine);

    mvReader.push_back(new CFileReaderBoundaryNode);
    mvReader.push_back(new CFileReaderBoundaryFace);
    mvReader.push_back(new CFileReaderBoundaryVolume);
    mvReader.push_back(new CFileReaderBoundaryEdge);

    mvReader.push_back(new CFileReaderBoundaryNodeMesh);
    mvReader.push_back(new CFileReaderBoundaryFaceMesh);
    mvReader.push_back(new CFileReaderBoundaryVolumeMesh);
    mvReader.push_back(new CFileReaderBoundaryEdgeMesh);

    mvReader.push_back(new CFileReaderMaterial);

    mvReader.push_back(new CFileReaderCommMesh);
    mvReader.push_back(new CFileReaderCommNode);
    mvReader.push_back(new CFileReaderCommElement);

    mvReader.push_back(new CFileReaderContactMesh);

    mvReader.push_back(new CFileReaderCommMesh2);
    mvReader.push_back(new CFileReaderCommFace);
    mvReader.push_back(new CFileReaderCommNodeCM2);

    mvReader.push_back(new CFileReaderElementGroup);
    mvReader.push_back(new CFileReaderElementGroupEntity);
    
}
// Destructor
//
CFileReaderChunk::~CFileReaderChunk()
{
    for_each(mvReader.begin(), mvReader.end(), pmw::DeleteObject());
}

// ------
// Method
// ------
//
//
void CFileReaderChunk::setFactory(pmw::CMeshFactory *pFactory)
{
    uint i;
    for(i=0; i<mvReader.size(); i++){
            mvReader[i]->setFactory(pFactory);
    };
}

// pMWメッシュ書式
// --
void CFileReaderChunk::Read(string filename)
{
    char c_Line[BUFFERLENGTH];
    string s_Line;

    mpLogger->Info(Utility::LoggerMode::MWDebug,"Input Filename => ",filename);//debug

    ifstream ifs(filename.c_str(),ios::in);

    if(ifs){
        while(!ifs.eof()){
            ifs.getline(c_Line,sizeof(c_Line),'\n');
            s_Line = c_Line;
            
            uint i; uint numOfBlock= mvReader.size();
            for(i=0; i < numOfBlock; i++){
                    mvReader[i]->Read(ifs, s_Line);
            };
        };
    }else{
        mpLogger->Info(Utility::LoggerMode::Error, "File not fount, filename => ", filename);
    }

    ifs.close();
}

// cntファイル読み込み
// --
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

// cntファイル名をフルパス名に変更
//
void CFileReaderChunk::setPath(string& filepath)
{
    string sFullPathName;

    sFullPathName = filepath + msCntFileName;
    msCntFileName = sFullPathName;
}









