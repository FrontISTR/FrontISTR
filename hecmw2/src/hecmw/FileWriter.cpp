/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileWriter.cpp
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "HEC_MPI.h"
#include "FileWriter.h"
using namespace FileIO;
CFileWriter::CFileWriter()
{
    mpGMGModel= pmw::CGMGModel::Instance();
}
CFileWriter::~CFileWriter()
{
    ;
}
void CFileWriter::setSolutionType(const uiint& nSolutionType)
{
    mnSolutionType = nSolutionType;
}
