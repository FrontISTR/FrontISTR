/* 
 * File:   FileReaderElementGroupEntity.h
 * Author: ktakeda
 *
 * Created on 2010/10/22, 13:07
 */
#include "FileReader.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

namespace FileIO{
#ifndef FILEREADER_ELEMENTGROUP_ENTITY_H
#define	FILEREADER_ELEMENTGROUP_ENTITY_H
class CFileReaderElementGroupEntity:public CFileReader{
public:
    CFileReaderElementGroupEntity();
    virtual ~CFileReaderElementGroupEntity();

public:
    virtual bool Read(ifstream& ifs, string& sline);
    virtual bool Read_bin(ifstream& ifs);
};
#endif	/* FILEREADERELEMENTGROUPENTITY_H */
}



