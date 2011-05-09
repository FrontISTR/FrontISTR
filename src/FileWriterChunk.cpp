//
//  FileWriterChunk.cpp
//
//
//
//                  2009.07.23
//                  2009.07.23
//                  k.Takeda
#include "FileWriterChunk.h"
using namespace FileIO;

CFileWriterChunk::CFileWriterChunk()
{
    mvWriter.reserve(3);

    mvWriter.push_back(new CFileWriterNode);
    mvWriter.push_back(new CFileWriterElement);
    mvWriter.push_back(new CFileWriterBoundaryNode);
}

CFileWriterChunk::~CFileWriterChunk()
{
    for_each(mvWriter.begin(), mvWriter.end(), pmw::DeleteObject());
}

// method
// --
void CFileWriterChunk::Write(string& filename, const uint& numOfLevel)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::MWDebug,"Output Filename => ",filename);//debug

    ofstream ofs(filename.c_str(),ios::out);

    uint i, ilevel;
    for(ilevel=0; ilevel< numOfLevel; ilevel++){

        //debug
        cout << "FileWriterChunk,  ilevel= " << ilevel << endl;

        for(i=0; i < mvWriter.size(); i++){
            mvWriter[i]->Write(ofs, ilevel);
        };
    };

//    //最終レベルだけ出力(debug?)
//    //
//    uint i, ilevel=numOfLevel-1;
//    //debug
//    cout << "FileWriterChunk,  ilevel= " << ilevel << endl;
//
//    for(i=0; i < mvWriter.size(); i++){
//        mvWriter[i]->Write(ofs, ilevel);
//    };


    ofs.close();
}






















