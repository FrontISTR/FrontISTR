/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileReaderBoundaryFace.h
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
#ifndef _FILEREADERBOUNDARYFACE_H_
#define	_FILEREADERBOUNDARYFACE_H_
class CFileReaderBoundaryFace:public CFileReader{
public:
    CFileReaderBoundaryFace();
    virtual ~CFileReaderBoundaryFace();
public:
   virtual bool Read(ifstream& ifs, string& sLine);
   virtual bool Read_bin(ifstream& ifs);
};
#endif	/* _FILEREADERBOUNDARYFACE_H */
}
