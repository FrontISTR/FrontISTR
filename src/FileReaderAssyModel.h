/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderAssyModel.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _FILEREADERASSYMODEL_H_9d31e0a0_2a33_4b85_b1bb_2af065f57954
#define	_FILEREADERASSYMODEL_H_9d31e0a0_2a33_4b85_b1bb_2af065f57954
#include "FileReader.h"
namespace FileIO{
class CFileReaderAssyModel:public CFileReader{
public:
    CFileReaderAssyModel(void);
    ~CFileReaderAssyModel(void);
public:
    virtual bool Read(ifstream& ifs, string& sline);
};
}
#endif	/* _FILEREADERASSYMODEL_H_9d31e0a0_2a33_4b85_b1bb_2af065f57954 */
