/* 
 * File:   FileWriterElement.h
 * Author: ktakeda
 *
 * Created on 2009/07/23, 17:43
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
    string StrType(const uiint& elem_type);
public:
    virtual void WriteDebug(ofstream& ofs, const uiint& mgLevel);
};
}

#endif	/* _FILEWRITERELEMENT_H */

