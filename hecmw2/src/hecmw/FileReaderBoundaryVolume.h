/* 
 * File:   FileReaderBoundaryVolume.h
 * Author: ktakeda
 *
 * Created on 2009/05/22, 18:13
 */
#include "FileReader.h"

namespace FileIO{
#ifndef _FILEREADERBOUNDARYVOLUME_H_
#define	_FILEREADERBOUNDARYVOLUME_H_   
class CFileReaderBoundaryVolume:public CFileReader{
public:
    CFileReaderBoundaryVolume();
    virtual ~CFileReaderBoundaryVolume();
public:
   virtual bool Read(ifstream& ifs, string& sLine);
};
#endif	/* _FILEREADERBOUNDARYVOLUME_H */
}


