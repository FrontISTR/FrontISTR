/* 
 * File:   FileWriterNode.h
 * Author: ktakeda
 *
 * Created on 2009/07/23, 13:28
 */

#ifndef _FILEWRITERNODE_H_f576c116_dcdf_4d04_bc13_e214549933b4
#define	_FILEWRITERNODE_H_f576c116_dcdf_4d04_bc13_e214549933b4

#include "AssyModel.h"
#include "Mesh.h"
#include "Node.h"

#include "FileWriter.h"

namespace FileIO{
class CFileWriterNode:public CFileWriter{
public:
    CFileWriterNode();
    virtual ~CFileWriterNode();

public:
    virtual void Write(ofstream& ofs, const uint& mgLevel);
};
}
#endif	/* _FILEWRITERNODE_H */

