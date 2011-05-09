/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderContactMesh.h
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
#include "boost/lexical_cast.hpp"
#include "ElementType.h"
namespace FileIO{
#ifndef _FILEREADERCONTACTMESH_H
#define	_FILEREADERCONTACTMESH_H
class CFileReaderContactMesh:public CFileReader{
public:
    CFileReaderContactMesh();
    virtual ~CFileReaderContactMesh();
    virtual bool Read(ifstream& ifs, string& sLine);
};
#endif	/* _FILEREADERCONTACTMESH_H */
}
