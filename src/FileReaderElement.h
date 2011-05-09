//
//	FileReaderElement.h
//
//				2008.12.08
//				2008.12.08
//				k.Takeda
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
//    virtual bool Read(ifstream& ifs, char* line);
    virtual bool Read(ifstream& ifs, string& sLine);
//    virtual bool ParseRead(ifstream& ifs, char* line);
};
}
#endif
