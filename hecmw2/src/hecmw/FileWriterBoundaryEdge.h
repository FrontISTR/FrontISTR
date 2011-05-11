/* 
 * File:   FileWriterBoundaryEdge.h
 * Author: ktakeda
 *
 * Created on 2010/05/07, 16:43
 */
#include "FileWriter.h"

namespace FileIO{
#ifndef _FILEWRITERBOUNDARYEDGE_H
#define	_FILEWRITERBOUNDARYEDGE_H
class CFileWriterBoundaryEdge:public CFileWriter{
public:
    CFileWriterBoundaryEdge();
    virtual ~CFileWriterBoundaryEdge();
public:
    virtual void WriteDebug(ofstream& ofs, const uiint& mgLevel);
};
#endif	/* _FILEWRITERBOUNDARYEDGE_H */
}

