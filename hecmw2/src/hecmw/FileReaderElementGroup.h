/* 
 * File:   FileReaderElementGroup.h
 * Author: ktakeda
 *
 * Created on 2010/10/22, 13:05
 */
#include "FileReader.h"

namespace FileIO{
#ifndef FILEREADER_ELEMENTGROUP_H
#define	FILEREADER_ELEMENTGROUP_H
class CFileReaderElementGroup:public CFileReader{
public:
    CFileReaderElementGroup();
    virtual ~CFileReaderElementGroup();
    
public:
    virtual bool Read(ifstream& ifs, string& sline);
};
#endif	/* FILEREADERELEMENTGROUP_H */
}

