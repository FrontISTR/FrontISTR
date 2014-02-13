/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterBoundaryNode.h
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
#ifndef _FILEWRITERBOUNDARYNODE_H_
#define	_FILEWRITERBOUNDARYNODE_H_
class CFileWriterBoundaryNode:public CFileWriter
{
public:
    CFileWriterBoundaryNode();
    virtual ~CFileWriterBoundaryNode();
public:
    virtual void WriteDebug(ofstream& ofs, const uiint& mgLevel);
};
#endif	/* _FILEWRITERBOUNDARYNODE_H */
}
