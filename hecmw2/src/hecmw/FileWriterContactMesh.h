/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileWriterContactMesh.h
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
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
