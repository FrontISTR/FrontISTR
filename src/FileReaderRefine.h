/* 
 * File:   FileReaderRefine.h
 * Author: ktakeda
 *
 * Created on 2009/04/08, 16:58
 */

#ifndef _FILEREADERREFINE_H_7c9bd15c_4fa3_4a86_9958_b82cf1dd8e0b
#define	_FILEREADERREFINE_H_7c9bd15c_4fa3_4a86_9958_b82cf1dd8e0b

#include "FileReader.h"
using namespace FileIO;

namespace FileIO{
class CFileReaderRefine:public CFileReader{
public:
    CFileReaderRefine();
    ~CFileReaderRefine();
public:
    virtual bool Read(ifstream &ifs, string &sline);
};
}
#endif	/* _FILEREADERREFINE_H_7c9bd15c_4fa3_4a86_9958_b82cf1dd8e0b */

