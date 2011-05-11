/* 
 * File:   FileReaderBoundaryNodeMesh.h
 * Author: ktakeda
 *
 * Created on 2010/04/28, 15:12
 */
#include "FileReader.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

namespace FileIO{
#ifndef _FILEREADERBOUNDARY_NODEMESH_H
#define	_FILEREADERBOUNDARY_NODEMESH_H
class CFileReaderBoundaryNodeMesh:public CFileReader{
public:
    CFileReaderBoundaryNodeMesh();
    virtual ~CFileReaderBoundaryNodeMesh();

public:
    virtual bool Read(ifstream& ifs, string& sline);
    virtual bool Read_bin(ifstream& ifs);
};
#endif	/* _FILEREADERBOUNDARY_NODEMESH_H */
}
