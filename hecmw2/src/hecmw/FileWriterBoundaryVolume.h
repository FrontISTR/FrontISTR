/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterBoundaryVolume.h
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
#include "FileWriter.h"
namespace FileIO
{
#ifndef _FILEWRITERBOUNDARYVOLUME_H
#define	_FILEWRITERBOUNDARYVOLUME_H
class CFileWriterBoundaryVolume:public CFileWriter
{
public:
    CFileWriterBoundaryVolume();
    virtual ~CFileWriterBoundaryVolume();
public:
    virtual void WriteDebug(ofstream& ofs, const uiint& mgLevel);
};
#endif	/* _FILEWRITERBOUNDARYVOLUME_H */
}
