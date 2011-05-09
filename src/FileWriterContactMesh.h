/* 
 * File:   FileWriterContactMesh.h
 * Author: ktakeda
 *
 * Created on 2009/11/20, 16:06
 */
#include "FileWriter.h"

namespace FileIO{
#ifndef _FILEWRITERCONTACTMESH_H
#define	_FILEWRITERCONTACTMESH_H
class CFileWriterContactMesh:public CFileWriter{
public:
    CFileWriterContactMesh();
    virtual ~CFileWriterContactMesh();

    virtual void Write(ofstream& ofs, const uint& mgLevel);
};
#endif	/* _FILEWRITERCONTACTMESH_H */
}


