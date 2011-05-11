/* 
 * File:   FileReaderBoundaryEdge.h
 * Author: ktakeda
 *
 * Created on 2010/04/28, 16:46
 */
#include "FileReader.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

namespace FileIO{

#ifndef _FILEREADERBOUNDARYEDGE_H
#define	_FILEREADERBOUNDARYEDGE_H
class CFileReaderBoundaryEdge:public CFileReader{
public:
    CFileReaderBoundaryEdge();
    virtual ~CFileReaderBoundaryEdge();
public:
    virtual bool Read(ifstream& ifs, string& sline);
    virtual bool Read_bin(ifstream& ifs);
};
#endif	/* _FILEREADERBOUNDARYEDGE_H */
}

