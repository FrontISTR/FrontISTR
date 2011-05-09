/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileWriterElement.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _FILEWRITERELEMENT_H_c12cc131_bab1_41f8_bf72_bed1a9257326
#define	_FILEWRITERELEMENT_H_c12cc131_bab1_41f8_bf72_bed1a9257326
#include "FileWriter.h"
namespace FileIO{
class CFileWriterElement:public CFileWriter{
public:
    CFileWriterElement();
    virtual ~CFileWriterElement();
private:
    string StrType(const uint& elem_type);
public:
    virtual void Write(ofstream& ofs, const uint& mgLevel);
};
}
#endif	/* _FILEWRITERELEMENT_H */
