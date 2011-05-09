/* 
 * File:   FileReaderBoundaryVolume.h
 * Author: ktakeda
 *
 * Created on 2009/05/22, 18:13
 */

#ifndef _FILEREADERBOUNDARYVOLUME_H_57126211_b33f_4c7c_9d94_ef379bfa95a8
#define	_FILEREADERBOUNDARYVOLUME_H_57126211_b33f_4c7c_9d94_ef379bfa95a8

#include "FileReader.h"

namespace FileIO{
class CFileReaderBoundaryVolume:public CFileReader{
public:
    CFileReaderBoundaryVolume();
    virtual ~CFileReaderBoundaryVolume();
public:
   virtual bool Read(ifstream& ifs, string& sLine);

};
}

#endif	/* _FILEREADERBOUNDARYVOLUME_H */

