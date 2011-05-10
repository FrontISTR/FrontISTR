/* 
 * File:   FileReaderCommElement.h
 * Author: ktakeda
 *
 * Created on 2009/09/18, 16:30
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

