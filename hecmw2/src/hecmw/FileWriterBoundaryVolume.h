/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileWriterBoundaryVolume.h
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
#include "FileWriter.h"
namespace FileIO{
#ifndef _FILEWRITERBOUNDARYVOLUME_H
#define	_FILEWRITERBOUNDARYVOLUME_H
class CFileWriterBoundaryVolume:public CFileWriter{
public:
    CFileWriterBoundaryVolume();
    virtual ~CFileWriterBoundaryVolume();
public:
    virtual void WriteDebug(ofstream& ofs, const uiint& mgLevel);
};
#endif	/* _FILEWRITERBOUNDARYVOLUME_H */
}
