/* 
 * File:   FileReaderBoundaryFace.h
 * Author: ktakeda
 *
 * Created on 2009/05/22, 18:09
 */
#include "FileReader.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

namespace FileIO{
#ifndef _FILEREADERBOUNDARYFACE_H_
#define	_FILEREADERBOUNDARYFACE_H_
class CFileReaderBoundaryFace:public CFileReader{
public:
    CFileReaderBoundaryFace();
    virtual ~CFileReaderBoundaryFace();

public:
   virtual bool Read(ifstream& ifs, string& sLine);
   virtual bool Read_bin(ifstream& ifs);

};
#endif	/* _FILEREADERBOUNDARYFACE_H */
}


