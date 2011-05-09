/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderCommNode_CM2.h
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
#ifndef _FILEREADERCOMMNODE_CM2_H
#define	_FILEREADERCOMMNODE_CM2_H
class CFileReaderCommNodeCM2:public CFileReader{
public:
    CFileReaderCommNodeCM2();
    virtual ~CFileReaderCommNodeCM2();
public:
    virtual bool Read(ifstream& ifs, string& sLine);
};
#endif	/* _FILEREADERCOMMNODE_CM2_H */
}
