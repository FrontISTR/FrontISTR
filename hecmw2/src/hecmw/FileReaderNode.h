//
//	FileReaderNode.h
//
//				2008.05.26
//				2008.12.08
//				k.Takeda
#ifndef FILE_TAG_NODE_HH_85DBAC1B_EAD2_4ba0_BBD8_2366AD615353
#define FILE_TAG_NODE_HH_85DBAC1B_EAD2_4ba0_BBD8_2366AD615353

#include "CommonStd.h"
#include "FileReader.h"

#include "NodeType.h"

namespace FileIO{
    class CFileReaderNode:public CFileReader{
public:
    CFileReaderNode();
    virtual ~CFileReaderNode();

public:
    virtual bool Read(ifstream& ifs, string& sLine);
};
}

#endif
