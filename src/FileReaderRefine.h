/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderRefine.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _FILEREADERREFINE_H_7c9bd15c_4fa3_4a86_9958_b82cf1dd8e0b
#define	_FILEREADERREFINE_H_7c9bd15c_4fa3_4a86_9958_b82cf1dd8e0b
#include "FileReader.h"
using namespace FileIO;
namespace FileIO{
class CFileReaderRefine:public CFileReader{
public:
    CFileReaderRefine();
    ~CFileReaderRefine();
public:
    virtual bool Read(ifstream &ifs, string &sline);
};
}
#endif	/* _FILEREADERREFINE_H_7c9bd15c_4fa3_4a86_9958_b82cf1dd8e0b */
