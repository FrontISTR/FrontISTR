/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileReaderNode.h
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
#ifndef FILE_TAG_NODE_HH_85DBAC1B_EAD2_4ba0_BBD8_2366AD615353
#define FILE_TAG_NODE_HH_85DBAC1B_EAD2_4ba0_BBD8_2366AD615353
#include "CommonStd.h"
#include "FileReader.h"
#include "NodeType.h"
#include "FileReaderBinCheck.h" 
namespace FileIO{
class CFileReaderNode:public CFileReader{
public:
    CFileReaderNode();
    virtual ~CFileReaderNode();
public:
    virtual bool Read(ifstream& ifs, string& sLine);
    virtual bool Read_bin(ifstream& ifs);
};
}
#endif
