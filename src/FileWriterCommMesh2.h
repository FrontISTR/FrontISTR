/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileWriterCommMesh2.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
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
