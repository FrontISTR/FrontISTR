/* 
 * File:   FileReaderAssyModel.h
 * Author: ktakeda
 *
 * Modify     2011/03/18
 * Created on 2009/04/08, 14:17
 */
#include "FileReader.h"

#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定
namespace FileIO{
#ifndef _FILEREADERASSYMODEL_H_
#define	_FILEREADERASSYMODEL_H_
class CFileReaderAssyModel:public CFileReader{
public:
    CFileReaderAssyModel(void);
    ~CFileReaderAssyModel(void);
public:
    virtual bool Read(ifstream& ifs, string& sline);
    virtual bool Read_bin(ifstream& ifs);
};
#endif	/* _FILEREADERASSYMODEL_H_ */
}


