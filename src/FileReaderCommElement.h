/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderCommElement.h
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
#ifndef _FILEREADERCOMMELEMENT_H
#define	_FILEREADERCOMMELEMENT_H
class CFileReaderCommElement:public CFileReader{
public:
    CFileReaderCommElement();
    virtual ~CFileReaderCommElement();
public:
    virtual bool Read(ifstream& ifs, string& sLine);
};
#endif	/* _FILEREADERCOMMELEMENT_H */
}
