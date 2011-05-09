/* 
 * File:   FileReaderCommFace.h
 * Author: ktakeda
 *
 * Created on 2010/03/12, 14:49
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



