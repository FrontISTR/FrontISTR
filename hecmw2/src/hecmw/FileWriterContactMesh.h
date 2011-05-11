/* 
 * File:   FileWriterContactMesh.h
 * Author: ktakeda
 *
 * Created on 2009/11/20, 16:06
 */
#include "FileWriter.h"

#include "HEC_MPI.h"

namespace FileIO{
#ifndef _FILEWRITERCONTACTMESH_H
#define	_FILEWRITERCONTACTMESH_H
class CFileWriterContactMesh:public CFileWriter{
public:
    CFileWriterContactMesh();
    virtual ~CFileWriterContactMesh();

    virtual void WriteDebug(ofstream& ofs, const uiint& mgLevel);
};
#endif	/* _FILEWRITERCONTACTMESH_H */
}



