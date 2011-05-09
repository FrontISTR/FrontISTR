/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderBoundaryFace.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _FILEREADERBOUNDARYFACE_H_f063dd36_8770_411d_9871_71e18964240e
#define	_FILEREADERBOUNDARYFACE_H_f063dd36_8770_411d_9871_71e18964240e
#include "FileReader.h"
namespace FileIO{
class CFileReaderBoundaryFace:public CFileReader{
public:
    CFileReaderBoundaryFace();
    virtual ~CFileReaderBoundaryFace();
public:
   virtual bool Read(ifstream& ifs, string& sLine);
};
}
#endif	/* _FILEREADERBOUNDARYFACE_H */
