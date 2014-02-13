/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderBoundaryNode.h
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
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
namespace FileIO
{
#ifndef _FILEREADERBOUNDARYNODE_H_
#define	_FILEREADERBOUNDARYNODE_H_
class CFileReaderBoundaryNode:public CFileReader
{
public:
    CFileReaderBoundaryNode();
    virtual ~CFileReaderBoundaryNode();
public:
    virtual bool Read(ifstream& ifs, string& sLine);
    virtual bool Read_bin(ifstream& ifs);

    virtual string Name();
};
#endif	/* _FILEREADERBOUNDARYNODE_H */
}
