/* 
 * File:   FileReaderBoundaryFaceMesh.h
 * Author: ktakeda
 *
 * Created on 2010/04/28, 15:13
 */
#include "FileReader.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

namespace FileIO{
#ifndef _FILEREADERBOUNDARYFACEMESH_H
#define	_FILEREADERBOUNDARYFACEMESH_H
class CFileReaderBoundaryFaceMesh:public CFileReader{
public:
    CFileReaderBoundaryFaceMesh();
    virtual ~CFileReaderBoundaryFaceMesh();

public:
    virtual bool Read(ifstream& ifs, string& sline);
    virtual bool Read_bin(ifstream& ifs);
};
#endif	/* _FILEREADERBOUNDARYFACEMESH_H */
}

