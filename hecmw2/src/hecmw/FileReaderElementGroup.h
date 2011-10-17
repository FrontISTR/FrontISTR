/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileReaderElementGroup.h
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FileReader.h"
#include "FileReaderBinCheck.h" 
namespace FileIO{
#ifndef FILEREADER_ELEMENTGROUP_H
#define	FILEREADER_ELEMENTGROUP_H
class CFileReaderElementGroup:public CFileReader{
public:
    CFileReaderElementGroup();
    virtual ~CFileReaderElementGroup();
public:
    virtual bool Read(ifstream& ifs, string& sline);
    virtual bool Read_bin(ifstream& ifs);
};
#endif	/* FILEREADERELEMENTGROUP_H */
}
