/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderElement.h
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "CommonStd.h"
#include "TypeDef.h"
#include "FileReader.h"
#include "FileReaderBinCheck.h"
namespace FileIO
{
#ifndef FILE_READER_ELEMENT_HH_F
#define FILE_READER_ELEMENT_HH_F
class CFileReaderElement:public CFileReader
{
public:
    CFileReaderElement();
    virtual ~CFileReaderElement();
public:
    virtual bool Read(ifstream& ifs, string& sLine);
    virtual bool Read_bin(ifstream& ifs);

    virtual string Name();
};
#endif
}
