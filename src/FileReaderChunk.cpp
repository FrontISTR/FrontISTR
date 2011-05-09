
#include <vector>

//
//  FileReaderChank.cpp
//
//			2009.05.01
//			2008.12.09
//			k.Takeda
#include "FileReaderChunk.h"
#include "FileReaderBoundaryNode.h"
#include "FileReaderBoundaryFace.h"
#include "FileReaderBoundaryVolume.h"
using namespace FileIO;

// Constructor
//
CFileReaderChunk::CFileReaderChunk()
{
    mpLogger = Utility::CLogger::Instance();

    mvReader.reserve(8);

    mvReader.push_back(new CFileReaderNode());
    mvReader.push_back(new CFileReaderElement());

    mvReader.push_back(new CFileReaderAssyModel);
    mvReader.push_back(new CFileReaderRefine);

    mvReader.push_back(new CFileReaderBoundaryNode);
    mvReader.push_back(new CFileReaderBoundaryFace);
    mvReader.push_back(new CFileReaderBoundaryVolume);

    mvReader.push_back(new CFileReaderMaterial);

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

// pMW
//
void CFileReaderChunk::Read(string filename)
{
    char c_Line[BUFFERLENGTH];
    string s_Line;

    mpLogger->Info(Utility::LoggerMode::MWDebug,"Input Filename => ",filename);//debug

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
        mpLogger->Info(Utility::LoggerMode::Error, "File not find, filename => ", filename);
    }

    ifs.close();
}
