/* 
 * File:   FileReaderBoundaryFace.h
 * Author: ktakeda
 *
 * Created on 2009/05/22, 18:09
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

