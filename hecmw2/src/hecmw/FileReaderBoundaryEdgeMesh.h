/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileReaderBoundaryEdgeMesh.h
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
#include "FileReaderBinCheck.h" 
namespace FileIO{
#ifndef _FILEREADERBOUNDARYEDGEMESH_H
#define	_FILEREADERBOUNDARYEDGEMESH_H
class CFileReaderBoundaryEdgeMesh:public CFileReader{
public:
    CFileReaderBoundaryEdgeMesh();
    virtual ~CFileReaderBoundaryEdgeMesh();
public:
    virtual bool Read(ifstream& ifs, string& sline);
    virtual bool Read_bin(ifstream& ifs);
};
#endif	/* _FILEREADERBOUNDARYEDGEMESH_H */
}
