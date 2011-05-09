/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileWriterBoundaryNode.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
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
