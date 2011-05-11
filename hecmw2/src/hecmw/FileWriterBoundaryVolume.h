/* 
 * File:   FileWriterBoundaryVolume.h
 * Author: ktakeda
 *
 * Created on 2010/05/07, 16:43
 */
#include "FileWriter.h"

namespace FileIO{
#ifndef _FILEWRITERBOUNDARYVOLUME_H
#define	_FILEWRITERBOUNDARYVOLUME_H
class CFileWriterBoundaryVolume:public CFileWriter{
public:
    CFileWriterBoundaryVolume();
    virtual ~CFileWriterBoundaryVolume();
public:
    virtual void WriteDebug(ofstream& ofs, const uiint& mgLevel);
};
#endif	/* _FILEWRITERBOUNDARYVOLUME_H */
}

