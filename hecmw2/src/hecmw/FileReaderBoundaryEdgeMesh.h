/* 
 * File:   FileReaderBoundaryEdgeMesh.h
 * Author: ktakeda
 *
 * Created on 2010/04/28, 15:15
 */
#include "FileReader.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

namespace FileIO{
#ifndef _FILEREADERBOUNDARYEDGEMESH_H
#define	_FILEREADERBOUNDARYEDGEMESH_H
class CFileReaderBoundaryEdgeMesh:public CFileReader{
public:
    CFileReaderBoundaryEdgeMesh();
    virtual ~CFileReaderBoundaryEdgeMesh();
public:
    virtual bool Read(ifstream& ifs, string& sline);
    virtual bool Read_bin(ifstream& ifs);
};
#endif	/* _FILEREADERBOUNDARYEDGEMESH_H */
}



