/* 
 * File:   FileReaderCnt.h
 * Author: ktakeda
 *
 * Created on 2009/09/22, 16:39
 */
#include "CommonFile.h"
#include "FileReader.h"

namespace FileIO{
#ifndef _FILEREADERCNT_H
#define	_FILEREADERCNT_H
class CFileReaderCnt:public CFileReader{
public:
    CFileReaderCnt();
    virtual ~CFileReaderCnt();

private:
    string msMeshFileBaseName;//メッシュデータのベースネーム

public:
    virtual bool Read(ifstream& ifs, string& sLine);
    string& getMeshFileBaseName(){ return msMeshFileBaseName;}
};
#endif	/* _FILEREADERCNT_H */
}



