/* 
 * File:   FileReaderCommMesh.h
 * Author: ktakeda
 *
 * Created on 2009/09/17, 17:13
 */
#include "FileReader.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

namespace FileIO{
#ifndef _FILEREADERCOMMMESH_H
#define	_FILEREADERCOMMMESH_H
class CFileReaderCommMesh:public CFileReader{
public:
    CFileReaderCommMesh();
    virtual ~CFileReaderCommMesh();

public:
    virtual bool Read(ifstream& ifs, string& sline);
    virtual bool Read_bin(ifstream& ifs);
};
#endif	/* _FILEREADERCOMMMESH_H */
}




