/* 
 * File:   FileWriterBoundaryNode.h
 * Author: ktakeda
 *
 * Created on 2009/07/23, 17:57
 */
#include "FileWriter.h"

namespace FileIO{
#ifndef _FILEWRITERBOUNDARYNODE_H_
#define	_FILEWRITERBOUNDARYNODE_H_    
class CFileWriterBoundaryNode:public CFileWriter{
public:
    CFileWriterBoundaryNode();
    virtual ~CFileWriterBoundaryNode();
public:
    virtual void Write(ofstream& ofs, const uint& mgLevel);
};
#endif	/* _FILEWRITERBOUNDARYNODE_H */
}



