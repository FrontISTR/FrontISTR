/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderCommFace.h
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
#ifndef _FILEREADERCOMMFACE_H
#define	_FILEREADERCOMMFACE_H
class CFileReaderCommFace:public CFileReader{
public:
    CFileReaderCommFace();
    virtual ~CFileReaderCommFace();
public:
    virtual bool Read(ifstream& ifs, string& sLine);
};
#endif	/* _FILEREADERCOMMFACE_H */
}
