/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderCnt.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "CommonFile.h"
#include "FileReader.h"
namespace FileIO{
#ifndef _FILEREADERCNT_H
#define	_FILEREADERCNT_H
class CFileReaderCnt:public CFileReader{
public:
    CFileReaderCnt();
    virtual ~CFileReaderCnt();
private:
    string msMeshFileBaseName;
public:
    virtual bool Read(ifstream& ifs, string& sLine);
    string& getMeshFileBaseName(){ return msMeshFileBaseName;}
};
#endif	/* _FILEREADERCNT_H */
}
