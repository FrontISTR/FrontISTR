/* 
 * File:   FileReaderAssyModel.h
 * Author: ktakeda
 *
 * Modify     2009/04/30
 * Created on 2009/04/08, 14:17
 */
#include "FileReader.h"

namespace FileIO{
#ifndef _FILEREADERASSYMODEL_H_9d31e0a0_2a33_4b85_b1bb_2af065f57954
#define	_FILEREADERASSYMODEL_H_9d31e0a0_2a33_4b85_b1bb_2af065f57954    
class CFileReaderAssyModel:public CFileReader{
public:
    CFileReaderAssyModel(void);
    ~CFileReaderAssyModel(void);
public:
    virtual bool Read(ifstream& ifs, string& sline);
};
#endif	/* _FILEREADERASSYMODEL_H_ */
}


