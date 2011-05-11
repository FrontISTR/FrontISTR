/* 
 * File:   FileReaderBoundaryVolumeMesh.h
 * Author: ktakeda
 *
 * Created on 2010/04/28, 15:14
 */
#include "FileReader.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

namespace FileIO{
#ifndef _FILEREADERBOUNDARYVOLUMEMESH_H
#define	_FILEREADERBOUNDARYVOLUMEMESH_H
class CFileReaderBoundaryVolumeMesh:public CFileReader{
public:
    CFileReaderBoundaryVolumeMesh();
    virtual ~CFileReaderBoundaryVolumeMesh();

public:
    virtual bool Read(ifstream& ifs, string& sline);
    virtual bool Read_bin(ifstream& ifs);
};
#endif	/* _FILEREADERBOUNDARYVOLUMEMESH_H */
}

