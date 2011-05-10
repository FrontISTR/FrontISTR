/* 
 * File:   FileReaderBoundaryEdge.h
 * Author: ktakeda
 *
 * Created on 2010/04/28, 16:46
 */
#include "FileReader.h"

namespace FileIO{
#ifndef _FILEREADERBOUNDARYEDGE_H
#define	_FILEREADERBOUNDARYEDGE_H
class CFileReaderBoundaryEdge:public CFileReader{
public:
    CFileReaderBoundaryEdge();
    virtual ~CFileReaderBoundaryEdge();
public:
    virtual bool Read(ifstream& ifs, string& sline);
};
#endif	/* _FILEREADERBOUNDARYEDGE_H */
}

