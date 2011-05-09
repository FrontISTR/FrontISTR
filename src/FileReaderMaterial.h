/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderMaterial.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _FILEREADERMATERIAL_H_7d6fbf50_4035_4524_b54b_dbf430da3b67
#define	_FILEREADERMATERIAL_H_7d6fbf50_4035_4524_b54b_dbf430da3b67
#include "FileReader.h"
#include "MatrialPropType.h"
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
