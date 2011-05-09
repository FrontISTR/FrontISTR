/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileWriterChunk.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FileWriterChunk.h"
using namespace FileIO;
CFileWriterChunk::CFileWriterChunk()
{
    mvWriter.reserve(5);
    mvWriter.push_back(new CFileWriterNode);
    mvWriter.push_back(new CFileWriterElement);
    mvWriter.push_back(new CFileWriterBoundaryNode);
    mvWriter.push_back(new CFileWriterContactMesh);
    mvWriter.push_back(new CFileWriterCommMesh2);
}
CFileWriterChunk::~CFileWriterChunk()
{
    for_each(mvWriter.begin(), mvWriter.end(), pmw::DeleteObject());
}
void CFileWriterChunk::Write(string& filename, const uint& numOfLevel)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::MWDebug,"Output Filename => ",filename);
    ofstream ofs(filename.c_str(),ios::out);
    uint i, ilevel;
    for(ilevel=0; ilevel< numOfLevel; ilevel++){
        cout << "FileWriterChunk,  ilevel= " << ilevel << endl;
        for(i=0; i < mvWriter.size(); i++){
            mvWriter[i]->Write(ofs, ilevel);
        };
    };
    ofs.close();
}
