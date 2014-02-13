/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderNode.h
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
#include "FileReader.h"
#include "NodeType.h"
#include "FileReaderBinCheck.h"

namespace FileIO
{
#ifndef FILE_TAG_NODE_HH_85DBAC1B_EAD2_4ba0_BBD8_2366AD615353
#define FILE_TAG_NODE_HH_85DBAC1B_EAD2_4ba0_BBD8_2366AD615353
class CFileReaderNode:public CFileReader
{
public:
    CFileReaderNode();
    virtual ~CFileReaderNode();
public:
    virtual bool Read(ifstream& ifs, string& sLine);
    virtual bool Read_bin(ifstream& ifs);

    virtual string Name();
};
#endif //FILE_TAG_NODE_HH_85DBAC1B_EAD2_4ba0_BBD8_2366AD615353
}

