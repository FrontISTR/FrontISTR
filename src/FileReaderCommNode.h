/* 
 * File:   FileReaderCommNode.h
 *
 * CommMesh内のノード情報読み込み(CommNodeという型は存在しない)
 *
 * Author: ktakeda
 *
 * Created on 2009/09/18, 15:34
 */
#include "FileReader.h"

namespace FileIO{
#ifndef _FILEREADERCOMMNODE_H
#define	_FILEREADERCOMMNODE_H
class CFileReaderCommNode:public CFileReader{
public:
    CFileReaderCommNode();
    virtual ~CFileReaderCommNode();

public:
    virtual bool Read(ifstream& ifs, string& sLine);
};
#endif	/* _FILEREADERCOMMNODE_H */
}





