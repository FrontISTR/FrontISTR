/* 
 * File:   FileReaderBoundaryNode.h
 * Author: ktakeda
 *
 * Created on 2009/05/22, 18:01
 */

#ifndef _FILEREADERBOUNDARYNODE_H_866b7894_b45d_41f8_842c_2516cac806b2
#define	_FILEREADERBOUNDARYNODE_H_866b7894_b45d_41f8_842c_2516cac806b2

#include "FileReader.h"
#include "BoundaryType.h"

namespace FileIO{
class CFileReaderBoundaryNode:public CFileReader{
public:
    CFileReaderBoundaryNode();
    virtual ~CFileReaderBoundaryNode();

public:
    virtual bool Read(ifstream& ifs, string& sLine);

};
}

#endif	/* _FILEREADERBOUNDARYNODE_H */

