/* 
 * File:   FileReaderCommMesh2.h
 * Author: ktakeda
 *
 * Created on 2010/03/12, 14:46
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




