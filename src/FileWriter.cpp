//
//  FileWriter.cpp
//
//
//
//                  2009.07.23
//                  2009.07.23
//                  k.Takeda
#include "FileWriter.h"
using namespace FileIO;

// construct & destruct
//
CFileWriter::CFileWriter()
{
    mpGMGModel= pmw::CGMGModel::Instance();
}
CFileWriter::~CFileWriter()
{
    ;
}

void CFileWriter::setSolutionType(const uint& nSolutionType)
{
    mnSolutionType = nSolutionType;
}



