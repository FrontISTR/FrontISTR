/* 
 * File:   FileWriterBoundaryNode.h
 * Author: ktakeda
 *
 * Created on 2009/07/23, 17:57
 */

#ifndef _FILEWRITERBOUNDARYNODE_H_127f4064_fe9e_44d9_9ec1_8e65c6a18d73
#define	_FILEWRITERBOUNDARYNODE_H_127f4064_fe9e_44d9_9ec1_8e65c6a18d73

#include "FileWriter.h"

namespace FileIO{
class CFileWriterBoundaryNode:public CFileWriter{
public:
    CFileWriterBoundaryNode();
    virtual ~CFileWriterBoundaryNode();

public:
    virtual void Write(ofstream& ofs, const uint& mgLevel);
};
}

#endif	/* _FILEWRITERBOUNDARYNODE_H */

