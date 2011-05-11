/* 
 * File:   FileReaderCommFace.h
 * Author: ktakeda
 *
 * Created on 2010/03/12, 14:49
 */
#include "FileReader.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

namespace FileIO{
#ifndef _FILEREADERCOMMFACE_H
#define	_FILEREADERCOMMFACE_H
class CFileReaderCommFace:public CFileReader{
public:
    CFileReaderCommFace();
    virtual ~CFileReaderCommFace();
    
public:
    virtual bool Read(ifstream& ifs, string& sLine);
    virtual bool Read_bin(ifstream& ifs);
};
#endif	/* _FILEREADERCOMMFACE_H */
}



