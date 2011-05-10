/* 
 * File:   FileReaderBoundaryEdgeMesh.h
 * Author: ktakeda
 *
 * Created on 2010/04/28, 15:15
 */
#include "FileReader.h"

namespace FileIO{
#ifndef _FILEREADERBOUNDARYEDGEMESH_H
#define	_FILEREADERBOUNDARYEDGEMESH_H
class CFileReaderBoundaryEdgeMesh:public CFileReader{
public:
    CFileReaderBoundaryEdgeMesh();
    virtual ~CFileReaderBoundaryEdgeMesh();
public:
    virtual bool Read(ifstream& ifs, string& sline);
};
#endif	/* _FILEREADERBOUNDARYEDGEMESH_H */
}



