/* 
 * File:   FileReaderBoundaryVolume.h
 * Author: ktakeda
 *
 * Created on 2009/05/22, 18:13
 */
#include "FileReader.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

namespace FileIO{
#ifndef _FILEREADERBOUNDARYVOLUME_H_
#define	_FILEREADERBOUNDARYVOLUME_H_   
class CFileReaderBoundaryVolume:public CFileReader{
public:
    CFileReaderBoundaryVolume();
    virtual ~CFileReaderBoundaryVolume();
public:
   virtual bool Read(ifstream& ifs, string& sLine);
   virtual bool Read_bin(ifstream& ifs);
};
#endif	/* _FILEREADERBOUNDARYVOLUME_H */
}


