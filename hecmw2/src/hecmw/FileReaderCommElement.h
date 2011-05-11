/* 
 * File:   FileReaderCommElement.h
 * Author: ktakeda
 *
 * Created on 2009/09/18, 16:30
 */
#include "FileReader.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

namespace FileIO{
#ifndef _FILEREADERCOMMELEMENT_H
#define	_FILEREADERCOMMELEMENT_H
class CFileReaderCommElement:public CFileReader{
public:
    CFileReaderCommElement();
    virtual ~CFileReaderCommElement();

public:
    virtual bool Read(ifstream& ifs, string& sLine);
    virtual bool Read_bin(ifstream& ifs);
};
#endif	/* _FILEREADERCOMMELEMENT_H */
}

