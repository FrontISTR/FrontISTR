/* 
 * File:   FileReaderRefine.h
 * Author: ktakeda
 *
 * Created on 2009/04/08, 16:58
 */
#include "FileReader.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

namespace FileIO{
#ifndef _FILEREADER_REFINE_H_
#define	_FILEREADER_REFINE_H_
class CFileReaderRefine:public CFileReader{
public:
    CFileReaderRefine();
    ~CFileReaderRefine();
public:
    virtual bool Read(ifstream &ifs, string &sline);
    virtual bool Read_bin(ifstream& ifs);
};
#endif	/* _FILEREADERREFINE_H_ */
}


