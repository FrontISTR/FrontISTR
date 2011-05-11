/* 
 * File:   FileWriterElementGroup.h
 * Author: ktakeda
 *
 * Created on 2010/10/26, 15:05
 */
#include "FileWriter.h"

namespace FileIO{
#ifndef FILEWRITERELEMENTGROUP_H
#define	FILEWRITERELEMENTGROUP_H
class CFileWriterElementGroup:public CFileWriter{
public:
    CFileWriterElementGroup();
    virtual ~CFileWriterElementGroup();

    virtual void WriteDebug(ofstream& ofs, const uiint& mgLevel);
};
#endif	/* FILEWRITERELEMENTGROUP_H */
}

