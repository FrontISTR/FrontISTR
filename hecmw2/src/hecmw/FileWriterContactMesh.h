/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterContactMesh.h
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "HEC_MPI.h"
#include "FileWriter.h"
namespace FileIO
{
#ifndef _FILEWRITERCONTACTMESH_H
#define	_FILEWRITERCONTACTMESH_H
class CFileWriterContactMesh:public CFileWriter
{
public:
    CFileWriterContactMesh();
    virtual ~CFileWriterContactMesh();
    virtual void WriteDebug(ofstream& ofs, const uiint& mgLevel);
};
#endif	/* _FILEWRITERCONTACTMESH_H */
}
