/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderBoundaryVolume.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
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
