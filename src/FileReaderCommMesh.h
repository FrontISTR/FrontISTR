/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderCommMesh.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FileReader.h"
namespace FileIO{
#ifndef _FILEREADERCOMMMESH_H
#define	_FILEREADERCOMMMESH_H
class CFileReaderCommMesh:public CFileReader{
public:
    CFileReaderCommMesh();
    virtual ~CFileReaderCommMesh();
public:
    virtual bool Read(ifstream& ifs, string& sline);
};
#endif	/* _FILEREADERCOMMMESH_H */
}
