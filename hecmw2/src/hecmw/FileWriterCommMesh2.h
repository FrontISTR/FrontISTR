/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterCommMesh2.h
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
#ifndef _FILEWRITERCOMMMESH2_H
#define	_FILEWRITERCOMMMESH2_H
class CFileWriterCommMesh2:public CFileWriter
{
public:
    CFileWriterCommMesh2();
    virtual ~CFileWriterCommMesh2();
public:
    virtual void WriteDebug(ofstream& ofs, const uiint& mgLevel);
};
#endif	/* _FILEWRITERCOMMMESH2_H */
}
