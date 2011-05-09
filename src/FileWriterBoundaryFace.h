/* 
 * File:   FileWriterBoundaryFace.h
 * Author: ktakeda
 *
 * Created on 2010/05/07, 16:42
 */
#include "FileWriter.h"

namespace FileIO{
#ifndef _FILEWRITERBOUNDARYFACE_H
#define	_FILEWRITERBOUNDARYFACE_H
class CFileWriterBoundaryFace:public CFileWriter{
public:
    CFileWriterBoundaryFace();
    virtual ~CFileWriterBoundaryFace();
public:
    virtual void Write(ofstream& ofs, const uint& mgLevel);
};
#endif	/* _FILEWRITERBOUNDARYFACE_H */
}

