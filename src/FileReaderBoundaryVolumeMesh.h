/* 
 * File:   FileReaderBoundaryVolumeMesh.h
 * Author: ktakeda
 *
 * Created on 2010/04/28, 15:14
 */
#include "FileReader.h"

namespace FileIO{
#ifndef _FILEREADERBOUNDARYVOLUMEMESH_H
#define	_FILEREADERBOUNDARYVOLUMEMESH_H
class CFileReaderBoundaryVolumeMesh:public CFileReader{
public:
    CFileReaderBoundaryVolumeMesh();
    virtual ~CFileReaderBoundaryVolumeMesh();

public:
    virtual bool Read(ifstream& ifs, string& sline);
};
#endif	/* _FILEREADERBOUNDARYVOLUMEMESH_H */
}

