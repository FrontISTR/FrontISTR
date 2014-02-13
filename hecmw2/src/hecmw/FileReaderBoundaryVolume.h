/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderBoundaryVolume.h
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
#ifndef _FILEREADERBOUNDARYVOLUME_H_
#define	_FILEREADERBOUNDARYVOLUME_H_
class CFileReaderBoundaryVolume:public CFileReader
{
public:
    CFileReaderBoundaryVolume();
    virtual ~CFileReaderBoundaryVolume();
public:
    virtual bool Read(ifstream& ifs, string& sLine);
    virtual bool Read_bin(ifstream& ifs);

    virtual string Name();
};
#endif	/* _FILEREADERBOUNDARYVOLUME_H */
}
