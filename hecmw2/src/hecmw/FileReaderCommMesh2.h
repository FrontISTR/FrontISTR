/* 
 * File:   FileReaderCommMesh2.h
 * Author: ktakeda
 *
 * Created on 2010/03/12, 14:46
 */
#include "FileReader.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

namespace FileIO{
#ifndef _FILEREADERCOMMMESH2_H
#define	_FILEREADERCOMMMESH2_H
class CFileReaderCommMesh2:public CFileReader{
public:
    CFileReaderCommMesh2();
    virtual ~CFileReaderCommMesh2();

public:
    virtual bool Read(ifstream& ifs, string& sLine);
    virtual bool Read_bin(ifstream& ifs);
};
#endif	/* _FILEREADERCOMMMESH2_H */
}




