/* 
 * File:   FileReaderMaterial.h
 * Author: ktakeda
 *
 * Created on 2009/07/27, 16:40
 */

#ifndef _FILEREADERMATERIAL_H_7d6fbf50_4035_4524_b54b_dbf430da3b67
#define	_FILEREADERMATERIAL_H_7d6fbf50_4035_4524_b54b_dbf430da3b67

#include "FileReader.h"

#include "MatrialType.h"

namespace FileIO{
class CFileReaderMaterial:public CFileReader{
public:
    CFileReaderMaterial();
    virtual ~CFileReaderMaterial();

public:
    virtual bool Read(ifstream& ifs, string& sline);
};
}
#endif	/* _FILEREADERMATERIAL_H */

