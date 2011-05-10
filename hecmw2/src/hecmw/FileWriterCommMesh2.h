/* 
 * File:   FileWriterCommMesh2.h
 * Author: ktakeda
 *
 * Created on 2010/03/15, 17:28
 */
#include "FileWriter.h"


namespace FileIO{
#ifndef _FILEWRITERCOMMMESH2_H
#define	_FILEWRITERCOMMMESH2_H
class CFileWriterCommMesh2:public CFileWriter{
public:
    CFileWriterCommMesh2();
    virtual ~CFileWriterCommMesh2();

public:
    virtual void Write(ofstream& ofs, const uint& mgLevel);
};
#endif	/* _FILEWRITERCOMMMESH2_H */
}




