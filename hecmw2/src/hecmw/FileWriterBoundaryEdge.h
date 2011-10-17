/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileWriterBoundaryEdge.h
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
#ifndef _FILEWRITERBOUNDARYEDGE_H
#define	_FILEWRITERBOUNDARYEDGE_H
class CFileWriterBoundaryEdge:public CFileWriter{
public:
    CFileWriterBoundaryEdge();
    virtual ~CFileWriterBoundaryEdge();
public:
    virtual void WriteDebug(ofstream& ofs, const uiint& mgLevel);
};
#endif	/* _FILEWRITERBOUNDARYEDGE_H */
}
