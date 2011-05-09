/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderCommMesh2.h
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
#ifndef _FILEREADERCOMMMESH2_H
#define	_FILEREADERCOMMMESH2_H
class CFileReaderCommMesh2:public CFileReader{
public:
    CFileReaderCommMesh2();
    virtual ~CFileReaderCommMesh2();
public:
    virtual bool Read(ifstream& ifs, string& sLine);
};
#endif	/* _FILEREADERCOMMMESH2_H */
}
