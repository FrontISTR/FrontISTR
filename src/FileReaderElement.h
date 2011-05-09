/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderElement.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef FILE_READER_ELEMENT_HH_F47C3D0E_0594_4815_A4B6_87EC43CFF668
#define FILE_READER_ELEMENT_HH_F47C3D0E_0594_4815_A4B6_87EC43CFF668
#include "CommonStd.h"
#include "TypeDef.h"
#include "FileReader.h"
namespace FileIO{
class CFileReaderElement:public CFileReader{
public:
    CFileReaderElement();
    virtual ~CFileReaderElement();
public:
    virtual bool Read(ifstream& ifs, string& sLine);
};
}
#endif
