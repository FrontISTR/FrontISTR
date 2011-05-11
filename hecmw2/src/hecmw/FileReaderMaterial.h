/* 
 * File:   FileReaderMaterial.h
 * Author: ktakeda
 *
 * Created on 2009/07/27, 16:40
 */
#include "FileReader.h"
#include "MatrialPropType.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

namespace FileIO{
#ifndef _FILEREADERMATERIAL_H_
#define	_FILEREADERMATERIAL_H_
class CFileReaderMaterial:public CFileReader{
public:
    CFileReaderMaterial();
    virtual ~CFileReaderMaterial();

public:
    virtual bool Read(ifstream& ifs, string& sline);
    virtual bool Read_bin(ifstream& ifs);
};
#endif	/* _FILEREADERMATERIAL_H */
}


