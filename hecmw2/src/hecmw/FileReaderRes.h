/* 
 * File:   FileReaderRes.h
 * Author: ktakeda
 *
 * リスタート・ファイル(*.res)
 *
 * Created on 2011/02/25, 16:09
 */
#include "FileReader.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

namespace FileIO{
#ifndef FILEREADERRES_H
#define	FILEREADERRES_H
class CFileReaderRes:public CFileReader{
public:
    CFileReaderRes();
    virtual ~CFileReaderRes();

public:
    virtual bool Read(ifstream& ifs, string& sLine);
    virtual bool Read_bin(ifstream& ifs);
};
#endif	/* FILEREADERRES_H */
}



