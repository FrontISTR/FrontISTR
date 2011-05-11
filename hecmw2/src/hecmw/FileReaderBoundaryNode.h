/* 
 * File:   FileReaderBoundaryNode.h
 * Author: ktakeda
 *
 * Created on 2009/05/22, 18:01
 */
#include "FileReader.h"
#include "BoundaryType.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

namespace FileIO{
#ifndef _FILEREADERBOUNDARYNODE_H_
#define	_FILEREADERBOUNDARYNODE_H_
class CFileReaderBoundaryNode:public CFileReader{
public:
    CFileReaderBoundaryNode();
    virtual ~CFileReaderBoundaryNode();

public:
    virtual bool Read(ifstream& ifs, string& sLine);
    virtual bool Read_bin(ifstream& ifs);

};
#endif	/* _FILEREADERBOUNDARYNODE_H */
}



