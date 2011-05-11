/* 
 * File:   FileReaderElementGroup.h
 * Author: ktakeda
 *
 * Created on 2010/10/22, 13:05
 */
#include "FileReader.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

namespace FileIO{
#ifndef FILEREADER_ELEMENTGROUP_H
#define	FILEREADER_ELEMENTGROUP_H
class CFileReaderElementGroup:public CFileReader{
public:
    CFileReaderElementGroup();
    virtual ~CFileReaderElementGroup();
    
public:
    virtual bool Read(ifstream& ifs, string& sline);
    virtual bool Read_bin(ifstream& ifs);
};
#endif	/* FILEREADERELEMENTGROUP_H */
}

