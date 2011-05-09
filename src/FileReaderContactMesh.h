/* 
 * File:   FileReaderContactMesh.h
 * Author: ktakeda
 *
 * Created on 2009/10/20, 16:59
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




