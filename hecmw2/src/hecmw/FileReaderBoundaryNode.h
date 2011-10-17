/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileReaderBoundaryNode.h
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FileReader.h"
#include "BoundaryType.h"
#include "FileReaderBinCheck.h" 
namespace FileIO{
#ifndef _FILEREADERBOUNDARYNODE_H_
#define	_FILEREADERBOUNDARYNODE_H_
class CFileReaderBoundaryNode:public CFileReader{
public:
    CFileReaderBoundaryNode();
    virtual ~CFileReaderBoundaryNode();
public:
    virtual bool Read(ifstream& ifs, string& sLine);
    virtual bool Read_bin(ifstream& ifs);
};
#endif	/* _FILEREADERBOUNDARYNODE_H */
}
