/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderCommMesh.h
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
#include "FileReaderBinCheck.h"
namespace FileIO
{
#ifndef _FILEREADERCOMMMESH_H
#define	_FILEREADERCOMMMESH_H
class CFileReaderCommMesh:public CFileReader
{
public:
    CFileReaderCommMesh();
    virtual ~CFileReaderCommMesh();
public:
    virtual bool Read(ifstream& ifs, string& sline);
    virtual bool Read_bin(ifstream& ifs);

    virtual string Name();
};
#endif	/* _FILEREADERCOMMMESH_H */
}
