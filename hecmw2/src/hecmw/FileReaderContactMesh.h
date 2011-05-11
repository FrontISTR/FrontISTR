/* 
 * File:   FileReaderContactMesh.h
 * Author: ktakeda
 *
 * Created on 2009/10/20, 16:59
 */
#include "FileReader.h"
#include "boost/lexical_cast.hpp"

#include "ElementType.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

namespace FileIO{
#ifndef _FILEREADERCONTACTMESH_H
#define	_FILEREADERCONTACTMESH_H
class CFileReaderContactMesh:public CFileReader{
public:
    CFileReaderContactMesh();
    virtual ~CFileReaderContactMesh();

    virtual bool Read(ifstream& ifs, string& sLine);
    virtual bool Read_bin(ifstream& ifs);
};
#endif	/* _FILEREADERCONTACTMESH_H */
}




