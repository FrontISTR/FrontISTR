/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterElement.h
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
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
namespace FileIO
{
class CFileWriterElement:public CFileWriter
{
public:
    CFileWriterElement();
    virtual ~CFileWriterElement();
private:
    string StrType(const uiint& elem_type);
public:
    virtual void WriteDebug(ofstream& ofs, const uiint& mgLevel);
};
}
#endif	/* _FILEWRITERELEMENT_H */
