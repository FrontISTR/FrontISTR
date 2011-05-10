/* 
 * File:   FileReaderElementGroupEntity.h
 * Author: ktakeda
 *
 * Created on 2010/10/22, 13:07
 */
#include "FileReader.h"


namespace FileIO{
#ifndef FILEREADER_ELEMENTGROUP_ENTITY_H
#define	FILEREADER_ELEMENTGROUP_ENTITY_H
class CFileReaderElementGroupEntity:public CFileReader{
public:
    CFileReaderElementGroupEntity();
    virtual ~CFileReaderElementGroupEntity();

public:
    virtual bool Read(ifstream& ifs, string& sline);
};
#endif	/* FILEREADERELEMENTGROUPENTITY_H */
}



