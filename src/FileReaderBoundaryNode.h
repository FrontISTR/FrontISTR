/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReaderBoundaryNode.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _FILEREADERBOUNDARYNODE_H_866b7894_b45d_41f8_842c_2516cac806b2
#define	_FILEREADERBOUNDARYNODE_H_866b7894_b45d_41f8_842c_2516cac806b2
#include "FileReader.h"
#include "BoundaryType.h"
namespace FileIO{
class CFileReaderBoundaryNode:public CFileReader{
public:
    CFileReaderBoundaryNode();
    virtual ~CFileReaderBoundaryNode();
public:
    virtual bool Read(ifstream& ifs, string& sLine);
};
}
#endif	/* _FILEREADERBOUNDARYNODE_H */
