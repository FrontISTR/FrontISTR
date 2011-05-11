//
//	FileReaderElement.h
//
//				2008.12.08
//				2008.12.08
//				k.Takeda
#include "CommonStd.h"
#include "TypeDef.h"

#include "FileReader.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

namespace FileIO{
#ifndef FILE_READER_ELEMENT_HH_F
#define FILE_READER_ELEMENT_HH_F
class CFileReaderElement:public CFileReader{
public:
    CFileReaderElement();
    virtual ~CFileReaderElement();

public:
    virtual bool Read(ifstream& ifs, string& sLine);
    virtual bool Read_bin(ifstream& ifs);

};
#endif
}

